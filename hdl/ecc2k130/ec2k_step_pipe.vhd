-- ec2k_step_pipe.vhd
-- One ECC2K-130 iteration, R -> R + sigma^j(R), built as a schedule around
-- ONE shared GF(2^131) multiplier.
--
--   j   = 3 + ((HW(x) >> 1) & 7)            HW = normal-basis Hamming weight
--   x2  = sigma^j(x),  y2 = sigma^j(y)       coordinate permutations, free
--   d   = x + x2
--   lam = (y + y2) / d
--   x3  = lam^2 + lam + d                    lam^2 = sigma(lam), free
--   y3  = lam (x + x3) + x3 + y
--
-- The division is Itoh-Tsujii.  Because squaring is wiring, 1/d costs only
-- the multiplies of the addition chain 1,2,4,8,16,32,64,128,130 for the
-- exponent 2^130 - 1: eight of them.  With the two of the addition, a step
-- is ten dependent multiplies -- against 8859 bitsliced instructions per
-- multiply on the GPU client, and 270 multiplies per inversion on a prime
-- field.  That the whole step is ten multiplies and nothing else is the
-- structural reason a binary Koblitz curve suits an FPGA.
--
-- Ten dependent multiplies cannot fill an 11-deep pipeline on their own, so
-- the unit keeps many independent steps in flight -- in the rho application
-- they are different walks -- one per slot.  Scheduling is dataflow rather
-- than a fixed phase rotation: when a product retires, its slot's next
-- operands are computed and the slot is pushed on a ready queue; every clock
-- the head of the queue issues.  New requests issue only when the queue is
-- empty, so work already in flight always has priority and the multiplier
-- never stalls while anything is ready.  With NSLOTS >= MUL_LATENCY + 2 it
-- saturates: one step retires every ten clocks.
--
-- Slot and chain state ride through the multiplier in its tag, so results
-- route back with no matching logic.  Degenerate inputs (d = 0) are not
-- special-cased, matching the client: 1/0 comes out as 0 and the step
-- produces garbage, as it does there; on the challenge curve the case has
-- probability about 2^-131 per step.
--
-- Per-slot state
--   x, y        the input point (needed until the last multiply)
--   j           3 bits, sigma^j is recomputed from x, y when needed
--   acc         running power of d, then lam
--   aux         d^3 (kept for the last chain step), then x3
--   st          which of the ten multiplies is next
--
-- Chain, by state (issue operands -> what the retiring product means)
--   0  d        * sigma^1 (d)     -> d^3        = beta_2
--   1  acc      * sigma^2 (acc)   -> d^15       = beta_4
--   2  acc      * sigma^4 (acc)   -> beta_8
--   3  acc      * sigma^8 (acc)   -> beta_16
--   4  acc      * sigma^16(acc)   -> beta_32
--   5  acc      * sigma^32(acc)   -> beta_64
--   6  acc      * sigma^64(acc)   -> beta_128
--   7  beta_2   * sigma^2 (acc)   -> beta_130 = d^(2^130 - 1)
--   8  (y + y2) * sigma^1 (acc)   -> lam      (sigma(beta_130) = 1/d)
--   9  lam      * (x + x3)        -> y3 = product + x3 + y

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

entity ec2k_step_pipe is
  generic (
    TAG_W     : natural := 8;
    SLOT_W    : natural := 4;            -- NSLOTS = 2**SLOT_W
    DP_WEIGHT : natural := DP_WEIGHT_DEFAULT
  );
  port (
    clk       : in  std_logic;
    rst       : in  std_logic;
    -- request: a point R; accepted when in_ready is high
    in_valid  : in  std_logic;
    in_ready  : out std_logic;
    in_x      : in  gf_t;
    in_y      : in  gf_t;
    in_tag    : in  std_logic_vector(TAG_W - 1 downto 0);
    -- result: R + sigma^j(R), its weight, and the distinguished-point test
    out_valid : out std_logic;
    out_x     : out gf_t;
    out_y     : out gf_t;
    out_hw    : out hw_t;
    out_dp    : out std_logic;
    out_tag   : out std_logic_vector(TAG_W - 1 downto 0)
  );
end entity;

architecture rtl of ec2k_step_pipe is

  constant NSLOTS : natural := 2 ** SLOT_W;
  constant ST_W   : natural := 4;
  constant MTAG_W : natural := SLOT_W + ST_W;

  subtype slot_t is unsigned(SLOT_W - 1 downto 0);
  subtype st_t   is unsigned(ST_W - 1 downto 0);
  subtype tag_t  is std_logic_vector(TAG_W - 1 downto 0);

  type gf_mem_t  is array (0 to NSLOTS - 1) of gf_t;
  type j_mem_t   is array (0 to NSLOTS - 1) of std_logic_vector(2 downto 0);
  type st_mem_t  is array (0 to NSLOTS - 1) of st_t;
  type tag_mem_t is array (0 to NSLOTS - 1) of tag_t;
  type q_t       is array (0 to NSLOTS - 1) of slot_t;

  function q_identity return q_t is
    variable q : q_t;
  begin
    for i in 0 to NSLOTS - 1 loop
      q(i) := to_unsigned(i, SLOT_W);
    end loop;
    return q;
  end function;

  -- per-slot context
  signal m_x, m_y, m_acc, m_aux : gf_mem_t := (others => (others => '0'));
  signal m_j   : j_mem_t   := (others => (others => '0'));
  signal m_st  : st_mem_t  := (others => (others => '0'));
  signal m_tag : tag_mem_t := (others => (others => '0'));

  -- ready queue: slots whose next multiply can issue
  signal rq : q_t := (others => (others => '0'));
  signal rq_wr, rq_rd : unsigned(SLOT_W downto 0) := (others => '0');

  -- free list: slots not holding a step
  signal fl : q_t := q_identity;
  signal fl_wr : unsigned(SLOT_W downto 0) := to_unsigned(NSLOTS, SLOT_W + 1);
  signal fl_rd : unsigned(SLOT_W downto 0) := (others => '0');

  -- pending request: one entry, so the weight is computed a clock before
  -- the operands are formed
  signal p_valid : std_logic := '0';
  signal p_x, p_y : gf_t := (others => '0');
  signal p_hw    : hw_t := (others => '0');
  signal p_tag   : tag_t := (others => '0');
  signal p_take  : std_logic;
  signal in_rdy  : std_logic;

  -- multiplier
  signal mul_valid : std_logic := '0';
  signal mul_a, mul_b : gf_t := (others => '0');
  signal mul_tag   : std_logic_vector(MTAG_W - 1 downto 0) := (others => '0');
  signal res_valid : std_logic;
  signal res_r     : gf_t;
  signal res_tag   : std_logic_vector(MTAG_W - 1 downto 0);

  -- output register stage (weight of x3 is taken from registered data)
  signal o1_valid : std_logic := '0';
  signal o1_x, o1_y : gf_t := (others => '0');
  signal o1_tag   : tag_t := (others => '0');

  signal rq_empty, fl_empty : boolean;

begin

  mul : entity work.gf131_mul
    generic map (TAG_W => MTAG_W)
    port map (
      clk => clk, rst => rst,
      in_valid => mul_valid, in_a => mul_a, in_b => mul_b, in_tag => mul_tag,
      out_valid => res_valid, out_r => res_r, out_tag => res_tag);

  rq_empty <= rq_wr = rq_rd;
  fl_empty <= fl_wr = fl_rd;

  -- the pending request issues when nothing in flight wants the multiplier
  p_take <= '1' when p_valid = '1' and rq_empty and not fl_empty else '0';
  in_rdy <= '1' when rst = '0' and (p_valid = '0' or p_take = '1') else '0';
  in_ready <= in_rdy;

  main : process (clk)
    variable s, rs  : natural range 0 to NSLOTS - 1;
    variable st     : st_t;
    variable rst_st : st_t;
    variable a, b   : gf_t;
    variable acc    : gf_t;
    variable d, lam, x3, y3 : gf_t;
    variable jm3    : std_logic_vector(2 downto 0);
    variable hw     : hw_t;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        rq_wr <= (others => '0'); rq_rd <= (others => '0');
        fl    <= q_identity;
        fl_wr <= to_unsigned(NSLOTS, SLOT_W + 1);
        fl_rd <= (others => '0');
        p_valid   <= '0';
        mul_valid <= '0';
        o1_valid  <= '0';
        out_valid <= '0';
      else
        mul_valid <= '0';
        o1_valid  <= '0';
        out_valid <= '0';

        -- ---------------- accept into the pending register --------------
        if p_take = '1' then
          p_valid <= '0';
        end if;
        if in_valid = '1' and in_rdy = '1' then
          p_valid <= '1';
          p_x     <= in_x;
          p_y     <= in_y;
          p_tag   <= in_tag;
          p_hw    <= gf_weight(in_x);
        end if;

        -- ---------------- issue one multiply ----------------------------
        if not rq_empty then
          s     := to_integer(rq(to_integer(rq_rd(SLOT_W - 1 downto 0))));
          rq_rd <= rq_rd + 1;
          st    := m_st(s);
          acc   := m_acc(s);
          case to_integer(st) is
            when 1      => a := acc;      b := gf_frob(acc, 2);
            when 2      => a := acc;      b := gf_frob(acc, 4);
            when 3      => a := acc;      b := gf_frob(acc, 8);
            when 4      => a := acc;      b := gf_frob(acc, 16);
            when 5      => a := acc;      b := gf_frob(acc, 32);
            when 6      => a := acc;      b := gf_frob(acc, 64);
            when 7      => a := m_aux(s); b := gf_frob(acc, 2);
            when 8      => a := m_y(s) xor gf_sigma_j(m_y(s), m_j(s));
                           b := gf_frob(acc, 1);
            when others => a := acc;      b := m_x(s) xor m_aux(s);
          end case;
          mul_a     <= a;
          mul_b     <= b;
          mul_tag   <= std_logic_vector(to_unsigned(s, SLOT_W)) & std_logic_vector(st);
          mul_valid <= '1';
        elsif p_take = '1' then
          s     := to_integer(fl(to_integer(fl_rd(SLOT_W - 1 downto 0))));
          fl_rd <= fl_rd + 1;
          jm3   := std_logic_vector(p_hw(3 downto 1));
          d     := p_x xor gf_sigma_j(p_x, jm3);
          m_x(s)   <= p_x;
          m_y(s)   <= p_y;
          m_j(s)   <= jm3;
          m_tag(s) <= p_tag;
          mul_a     <= d;
          mul_b     <= gf_frob(d, 1);
          mul_tag   <= std_logic_vector(to_unsigned(s, SLOT_W)) & std_logic_vector(to_unsigned(0, ST_W));
          mul_valid <= '1';
        end if;

        -- ---------------- retire one multiply ----------------------------
        if res_valid = '1' then
          rs     := to_integer(unsigned(res_tag(MTAG_W - 1 downto ST_W)));
          rst_st := unsigned(res_tag(ST_W - 1 downto 0));
          case to_integer(rst_st) is
            when 0 =>
              m_acc(rs) <= res_r;
              m_aux(rs) <= res_r;                       -- beta_2 = d^3
              m_st(rs)  <= to_unsigned(1, ST_W);
              rq(to_integer(rq_wr(SLOT_W - 1 downto 0))) <= to_unsigned(rs, SLOT_W);
              rq_wr <= rq_wr + 1;
            when 1 to 7 =>
              m_acc(rs) <= res_r;
              m_st(rs)  <= rst_st + 1;
              rq(to_integer(rq_wr(SLOT_W - 1 downto 0))) <= to_unsigned(rs, SLOT_W);
              rq_wr <= rq_wr + 1;
            when 8 =>
              lam := res_r;
              d   := m_x(rs) xor gf_sigma_j(m_x(rs), m_j(rs));
              x3  := gf_frob(lam, 1) xor lam xor d;
              m_acc(rs) <= lam;
              m_aux(rs) <= x3;
              m_st(rs)  <= to_unsigned(9, ST_W);
              rq(to_integer(rq_wr(SLOT_W - 1 downto 0))) <= to_unsigned(rs, SLOT_W);
              rq_wr <= rq_wr + 1;
            when others =>
              y3 := res_r xor m_aux(rs) xor m_y(rs);
              o1_x     <= m_aux(rs);
              o1_y     <= y3;
              o1_tag   <= m_tag(rs);
              o1_valid <= '1';
              fl(to_integer(fl_wr(SLOT_W - 1 downto 0))) <= to_unsigned(rs, SLOT_W);
              fl_wr <= fl_wr + 1;
          end case;
        end if;

        -- ---------------- output stage: weight and DP test ---------------
        if o1_valid = '1' then
          hw := gf_weight(o1_x);
          out_x     <= o1_x;
          out_y     <= o1_y;
          out_tag   <= o1_tag;
          out_hw    <= hw;
          if to_integer(hw) <= DP_WEIGHT then
            out_dp <= '1';
          else
            out_dp <= '0';
          end if;
          out_valid <= '1';
        end if;
      end if;
    end if;
  end process;

end architecture;
