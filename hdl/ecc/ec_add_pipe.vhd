-- ec_add_pipe.vhd
-- Affine elliptic-curve point addition over secp256k1, built as a schedule
-- around ONE shared modular multiplier.
--
--   lam = (y2 - y1) * inv        inv = 1/(x2 - x1), supplied by the caller
--   x3  = lam^2 - x1 - x2
--   y3  = lam*(x1 - x3) - y1
--
-- Three dependent multiplies, so a single addition cannot use a pipelined
-- multiplier on its own: with latency L = 12 it would idle 11 cycles out of
-- 12.  The unit therefore keeps many *independent* additions in flight -- in
-- the Pollard-rho application they are different walks -- and issues one
-- multiply every clock, cycling between the three phases.  One addition is
-- accepted every 3 clocks and the multiplier never stalls, so the whole
-- engine costs exactly 3 modular multiplies per point addition.
--
-- Requests carry a caller tag; each is assigned a slot from a rotating pool
-- and the slot travels through the multiplier in its tag field, so results
-- route back to their context with no matching logic.  A slot is reused only
-- after 3*NSLOTS clocks, comfortably longer than the 3L + 3 clocks an
-- addition lives, so no arbitration or back-pressure is needed.
--
-- The inverse comes from outside on purpose: in the rho engine one inversion
-- is shared by a whole batch of walks (Montgomery's trick), which is what
-- makes affine coordinates cheaper than projective ones here.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.fp_pkg.all;

entity ec_add_pipe is
  generic (
    TAG_W  : natural := 8;
    NSLOTS : natural := 16          -- power of two, >= (3*MUL_LATENCY+6)/3
  );
  port (
    clk       : in  std_logic;
    rst       : in  std_logic;
    -- request: accepted when in_ready is high
    in_valid  : in  std_logic;
    in_ready  : out std_logic;
    in_x1     : in  fp_t;
    in_y1     : in  fp_t;
    in_x2     : in  fp_t;
    in_y2     : in  fp_t;
    in_inv    : in  fp_t;
    in_tag    : in  std_logic_vector(TAG_W - 1 downto 0);
    -- result
    out_valid : out std_logic;
    out_x3    : out fp_t;
    out_y3    : out fp_t;
    out_tag   : out std_logic_vector(TAG_W - 1 downto 0)
  );
end entity;

architecture rtl of ec_add_pipe is

  constant SLOT_W : natural := 4;                  -- log2(NSLOTS)
  constant MTAG_W : natural := SLOT_W + 2;         -- slot + phase

  subtype slot_t is unsigned(SLOT_W - 1 downto 0);

  type fp_mem_t  is array (0 to NSLOTS - 1) of fp_t;
  type tag_mem_t is array (0 to NSLOTS - 1) of std_logic_vector(TAG_W - 1 downto 0);

  -- per-slot context
  signal m_x1, m_y1, m_x2 : fp_mem_t := (others => (others => '0'));
  signal m_lam, m_x3, m_d : fp_mem_t := (others => (others => '0'));
  signal m_tag            : tag_mem_t := (others => (others => '0'));

  -- multiplier interface
  signal mul_valid : std_logic := '0';
  signal mul_a, mul_b : fp_t := (others => '0');
  signal mul_tag   : std_logic_vector(MTAG_W - 1 downto 0) := (others => '0');
  signal res_valid : std_logic;
  signal res_r     : fp_t;
  signal res_tag   : std_logic_vector(MTAG_W - 1 downto 0);

  -- phase-2 and phase-3 work queues (slots ready to issue)
  type q_t is array (0 to NSLOTS - 1) of slot_t;
  signal q2, q3 : q_t := (others => (others => '0'));
  signal q2_wr, q2_rd, q3_wr, q3_rd : unsigned(SLOT_W downto 0) := (others => '0');

  signal phase    : unsigned(1 downto 0) := (others => '0');
  signal next_slot : slot_t := (others => '0');

  function q_empty (wr, rd : unsigned) return boolean is
  begin
    return wr = rd;
  end function;

begin

  mul : entity work.fp_mul_secp256k1
    generic map (TAG_W => MTAG_W)
    port map (
      clk => clk, rst => rst,
      in_valid => mul_valid, in_a => mul_a, in_b => mul_b, in_tag => mul_tag,
      out_valid => res_valid, out_r => res_r, out_tag => res_tag);

  -- a new request is taken only on phase 0, and only if no phase-2 or
  -- phase-3 work is waiting (those have priority: they are already in
  -- flight and holding a slot)
  in_ready <= '1' when phase = 0 and rst = '0' else '0';

  main : process (clk)
    variable rslot  : integer;
    variable rphase : unsigned(1 downto 0);
    variable x3v, dv, y3v : fp_t;
    variable s : integer;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        phase     <= (others => '0');
        next_slot <= (others => '0');
        q2_wr <= (others => '0'); q2_rd <= (others => '0');
        q3_wr <= (others => '0'); q3_rd <= (others => '0');
        mul_valid <= '0';
        out_valid <= '0';
      else
        mul_valid <= '0';
        out_valid <= '0';

        -- ---------------- issue one multiply, rotating phases ----------
        case to_integer(phase) is
          when 0 =>
            if in_valid = '1' then
              s := to_integer(next_slot);
              m_x1(s)  <= in_x1;
              m_y1(s)  <= in_y1;
              m_x2(s)  <= in_x2;
              m_tag(s) <= in_tag;
              mul_a     <= fp_sub(in_y2, in_y1);
              mul_b     <= in_inv;
              mul_tag   <= std_logic_vector(next_slot) & "00";
              mul_valid <= '1';
              next_slot <= next_slot + 1;
            end if;
          when 1 =>
            if not q_empty(q2_wr, q2_rd) then
              s := to_integer(q2(to_integer(q2_rd(SLOT_W - 1 downto 0))));
              mul_a     <= m_lam(s);
              mul_b     <= m_lam(s);
              mul_tag   <= std_logic_vector(to_unsigned(s, SLOT_W)) & "01";
              mul_valid <= '1';
              q2_rd <= q2_rd + 1;
            end if;
          when others =>
            if not q_empty(q3_wr, q3_rd) then
              s := to_integer(q3(to_integer(q3_rd(SLOT_W - 1 downto 0))));
              mul_a     <= m_lam(s);
              mul_b     <= m_d(s);
              mul_tag   <= std_logic_vector(to_unsigned(s, SLOT_W)) & "10";
              mul_valid <= '1';
              q3_rd <= q3_rd + 1;
            end if;
        end case;

        if phase = 2 then
          phase <= (others => '0');
        else
          phase <= phase + 1;
        end if;

        -- ---------------- retire one multiply --------------------------
        if res_valid = '1' then
          rslot  := to_integer(unsigned(res_tag(MTAG_W - 1 downto 2)));
          rphase := unsigned(res_tag(1 downto 0));
          case to_integer(rphase) is
            when 0 =>
              -- lam
              m_lam(rslot) <= res_r;
              q2(to_integer(q2_wr(SLOT_W - 1 downto 0))) <= to_unsigned(rslot, SLOT_W);
              q2_wr <= q2_wr + 1;
            when 1 =>
              -- lam^2 -> x3 = lam^2 - x1 - x2, d = x1 - x3
              x3v := fp_sub(fp_sub(res_r, m_x1(rslot)), m_x2(rslot));
              dv  := fp_sub(m_x1(rslot), x3v);
              m_x3(rslot) <= x3v;
              m_d(rslot)  <= dv;
              q3(to_integer(q3_wr(SLOT_W - 1 downto 0))) <= to_unsigned(rslot, SLOT_W);
              q3_wr <= q3_wr + 1;
            when others =>
              -- lam*(x1-x3) -> y3
              y3v       := fp_sub(res_r, m_y1(rslot));
              out_x3    <= m_x3(rslot);
              out_y3    <= y3v;
              out_tag   <= m_tag(rslot);
              out_valid <= '1';
          end case;
        end if;
      end if;
    end if;
  end process;

end architecture;
