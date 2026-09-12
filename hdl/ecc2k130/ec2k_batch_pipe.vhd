-- ec2k_batch_pipe.vhd
-- ECC2K-130 iterations, R -> R + sigma^j(R), for W walks at a time around
-- ONE shared GF(2^131) multiplier, with the W inversions batched through a
-- product tree (Montgomery's trick).
--
-- Per step the arithmetic is
--
--   d   = x + sigma^j(x)                 j = 3 + ((HW(x) >> 1) & 7)
--   lam = (y + sigma^j(y)) / d
--   x3  = lam^2 + lam + d                lam^2 is a permutation, free
--   y3  = lam (x + x3) + x3 + y
--
-- and the one expensive thing is 1/d.  ec2k_step_pipe inverts every d on
-- its own, eight multiplies each, ten per step.  Here W walks form a batch
-- and their d's the leaves of a binary tree:
--
--   forward    t(n) = t(2n) t(2n+1)            W-1 multiplies, log W levels
--   invert     t(1) = 1/t(1)                    8 multiplies (Itoh-Tsujii)
--   backward   t(2n) = t(n) t(2n+1),            2(W-1) multiplies, log W levels
--              t(2n+1) = t(n) t(2n)             (leaves end holding 1/d_i)
--   lam        W multiplies
--   y3         W multiplies
--
-- 5W + 5 multiplies per batch: 5 + 5/W per step, 5.3 at W = 16, against 10.
-- That is the same trick the GPU client plays with 32 walks per word, and
-- it is the single largest lever on throughput in this design.
--
-- Scheduling.  Each tree level, each inversion step, and each of the lam
-- and y3 passes is a *burst* of independent multiplies that can issue on
-- consecutive clocks.  A batch is a sequence of 2 log W + 10 bursts; when
-- the last product of a burst retires the batch is pushed on a ready queue
-- with its next burst, and the issue engine streams bursts from the queue
-- head with the following batch prefetched so no clock is lost between
-- them.  Because the multiplier is a fixed-latency pipe, products retire in
-- issue order, so "last of the burst retired" means the whole burst has.
-- Bursts from different batches interleave freely; a handful of batches in
-- flight keeps the multiplier saturated (measured: see README).
--
-- Overwriting the tree in place is safe because every read of a value that
-- a burst will overwrite happens at issue, and the overwrite happens at
-- retire MUL_LATENCY clocks later; this needs MUL_LATENCY >= 2.
--
-- Memories.  Each array has exactly one writer so it maps to distributed
-- RAM rather than flip-flops: x, y, j, tag, valid are written only when a
-- batch fills; the tree only when a product retires; z (which holds d, then
-- x3) only from the issue side.  Reads happen one clock before operands
-- are formed, so no path runs RAM -> sigma -> XOR in one clock.
--
-- Batches shorter than W -- the tail of a run, or a testbench -- would wait
-- forever for leaves that never come, so a batch that has been partly
-- filled for FLUSH_CLK clocks with nothing arriving is padded with dummy
-- leaves (weight-1 x, d /= 0) that produce no output.
--
-- Slot memory per batch of W: W (x, y, z) + 2W tree words, all 131 bits.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

entity ec2k_batch_pipe is
  generic (
    TAG_W     : natural := 8;
    LOG_W     : natural := 4;                -- W = 2**LOG_W walks per batch, >= 1
    LOG_NB    : natural := 3;                -- 2**LOG_NB batches in flight
    FLUSH_CLK : natural := 32;
    DP_WEIGHT : natural := DP_WEIGHT_DEFAULT
  );
  port (
    clk       : in  std_logic;
    rst       : in  std_logic;
    in_valid  : in  std_logic;
    in_ready  : out std_logic;
    in_x      : in  gf_t;
    in_y      : in  gf_t;
    in_tag    : in  std_logic_vector(TAG_W - 1 downto 0);
    out_valid : out std_logic;
    out_x     : out gf_t;
    out_y     : out gf_t;
    out_hw    : out hw_t;
    out_dp    : out std_logic;
    out_tag   : out std_logic_vector(TAG_W - 1 downto 0)
  );
end entity;

architecture rtl of ec2k_batch_pipe is

  constant W  : natural := 2 ** LOG_W;
  constant NB : natural := 2 ** LOG_NB;

  function maxn (a, b : natural) return natural is
  begin
    if a > b then return a; else return b; end if;
  end function;

  -- node/leaf index in the tag: a tree node (LOG_W+1 bits) or an inversion
  -- step (3 bits)
  constant IDX_W  : natural := maxn(LOG_W + 1, 3);
  constant LVL_W  : natural := maxn(LOG_W, 3);
  constant KIND_W : natural := 3;
  constant MTAG_W : natural := LOG_NB + KIND_W + IDX_W + 1;

  constant PH_FWD : natural := 0;
  constant PH_INV : natural := 1;
  constant PH_BWD : natural := 2;
  constant PH_LAM : natural := 3;
  constant PH_FIN : natural := 4;

  subtype bid_t  is unsigned(LOG_NB - 1 downto 0);
  subtype leaf_t is unsigned(LOG_W - 1 downto 0);
  subtype node_t is unsigned(LOG_W downto 0);
  subtype idx_t  is unsigned(IDX_W - 1 downto 0);
  subtype lvl_t  is unsigned(LVL_W - 1 downto 0);
  subtype ph_t   is unsigned(KIND_W - 1 downto 0);
  subtype tag_t  is std_logic_vector(TAG_W - 1 downto 0);
  subtype mtag_t is std_logic_vector(MTAG_W - 1 downto 0);

  type gf_leaf_mem_t is array (0 to NB * W - 1) of gf_t;
  type gf_tree_mem_t is array (0 to NB * 2 * W - 1) of gf_t;
  type j_mem_t   is array (0 to NB * W - 1) of std_logic_vector(2 downto 0);
  type tag_mem_t is array (0 to NB * W - 1) of tag_t;
  type v_mem_t   is array (0 to NB * W - 1) of std_logic;
  type ph_mem_t  is array (0 to NB - 1) of ph_t;
  type lvl_mem_t is array (0 to NB - 1) of lvl_t;
  type q_t       is array (0 to NB - 1) of bid_t;

  function q_identity return q_t is
    variable q : q_t;
  begin
    for i in 0 to NB - 1 loop
      q(i) := to_unsigned(i, LOG_NB);
    end loop;
    return q;
  end function;

  function leaf_addr (b : bid_t; i : unsigned) return natural is
  begin
    return to_integer(b & i(LOG_W - 1 downto 0));
  end function;

  function tree_addr (b : bid_t; n : node_t) return natural is
  begin
    return to_integer(b & n);
  end function;

  -- last index of a burst
  function burst_end (ph : ph_t; lvl : lvl_t) return idx_t is
  begin
    case to_integer(ph) is
      when PH_FWD => return to_unsigned(2 ** to_integer(lvl) - 1, IDX_W);
      when PH_INV => return to_unsigned(0, IDX_W);
      when PH_BWD => return to_unsigned(2 ** (to_integer(lvl) + 1) - 1, IDX_W);
      when others => return to_unsigned(W - 1, IDX_W);
    end case;
  end function;

  -- Frobenius power for inversion step k: beta_{2k} = beta_k sigma^k(beta_k)
  function inv_frob (v : gf_t; k : unsigned(2 downto 0)) return gf_t is
  begin
    case to_integer(k) is
      when 0      => return gf_frob(v, 1);
      when 1      => return gf_frob(v, 2);
      when 2      => return gf_frob(v, 4);
      when 3      => return gf_frob(v, 8);
      when 4      => return gf_frob(v, 16);
      when 5      => return gf_frob(v, 32);
      when 6      => return gf_frob(v, 64);
      when others => return gf_frob(v, 2);
    end case;
  end function;

  -- ------------------------------------------------------------------ --
  -- memories
  -- ------------------------------------------------------------------ --
  signal lx, ly, lz : gf_leaf_mem_t := (others => (others => '0'));
  signal lj   : j_mem_t   := (others => (others => '0'));
  signal ltag : tag_mem_t := (others => (others => '0'));
  signal lval : v_mem_t   := (others => '0');
  signal tr   : gf_tree_mem_t := (others => (others => '0'));

  signal b_ph  : ph_mem_t  := (others => (others => '0'));
  signal b_lvl : lvl_mem_t := (others => (others => '0'));

  -- ------------------------------------------------------------------ --
  -- input: two pipeline registers so the weight has two clocks
  -- ------------------------------------------------------------------ --
  signal p0_valid, p1_valid : std_logic := '0';
  signal p0_x, p0_y, p1_x, p1_y : gf_t := (others => '0');
  signal p0_tag, p1_tag : tag_t := (others => '0');
  signal p0_parts : hw_parts_t := (others => (others => '0'));
  signal p1_hw    : hw_t := (others => '0');
  signal p0_adv, p1_take : std_logic;
  signal in_rdy   : std_logic;

  -- ------------------------------------------------------------------ --
  -- fill
  -- ------------------------------------------------------------------ --
  signal fb        : bid_t := (others => '0');
  signal fb_valid  : std_logic := '0';
  signal fill_cnt  : unsigned(LOG_W downto 0) := (others => '0');
  signal fill_pend : std_logic := '0';       -- a full batch waiting for the queue
  signal pend_b    : bid_t := (others => '0'); -- which one (fb is reallocated meanwhile)
  signal idle_cnt  : unsigned(15 downto 0) := (others => '0');
  signal flushing  : std_logic := '0';
  signal fill_ok, dummy_fill : std_logic;

  -- ready queue and free list of batch ids
  signal rq : q_t := (others => (others => '0'));
  signal rq_wr, rq_rd : unsigned(LOG_NB downto 0) := (others => '0');
  signal fl : q_t := q_identity;
  signal fl_wr : unsigned(LOG_NB downto 0) := to_unsigned(NB, LOG_NB + 1);
  signal fl_rd : unsigned(LOG_NB downto 0) := (others => '0');
  signal rq_empty, fl_empty : boolean;

  -- ------------------------------------------------------------------ --
  -- burst engine: cur issues, nxt is prefetched
  -- ------------------------------------------------------------------ --
  signal cur_valid, nxt_valid : std_logic := '0';
  signal cur_b, nxt_b     : bid_t := (others => '0');
  signal cur_ph, nxt_ph   : ph_t := (others => '0');
  signal cur_lvl, nxt_lvl : lvl_t := (others => '0');
  signal cur_idx, cur_end : idx_t := (others => '0');

  -- stage A: raw operands read from memory
  signal ra_valid : std_logic := '0';
  signal ra_a, ra_b, ra_c : gf_t := (others => '0');
  signal ra_ph    : ph_t := (others => '0');
  signal ra_leafy : std_logic := '0';        -- this level touches leaves
  signal ra_ja, ra_jb : std_logic_vector(2 downto 0) := (others => '0');
  signal ra_k     : unsigned(2 downto 0) := (others => '0');
  signal ra_tag   : mtag_t := (others => '0');
  signal ra_zaddr : natural range 0 to NB * W - 1 := 0;

  -- multiplier
  signal mul_valid : std_logic := '0';
  signal mul_a, mul_b : gf_t := (others => '0');
  signal mul_tag   : mtag_t := (others => '0');
  signal res_valid : std_logic;
  signal res_r     : gf_t;
  signal res_tag   : mtag_t;

  -- output: two registers so the weight of x3 has two clocks
  signal o1_valid, o2_valid : std_logic := '0';
  signal o1_x, o1_y, o2_x, o2_y : gf_t := (others => '0');
  signal o1_tag, o2_tag : tag_t := (others => '0');
  signal o2_parts : hw_parts_t := (others => (others => '0'));

begin

  assert LOG_W >= 1 report "ec2k_batch_pipe needs at least two walks per batch" severity failure;
  assert MUL_LATENCY >= 2 report "in-place tree update needs MUL_LATENCY >= 2" severity failure;

  mul : entity work.gf131_mul
    generic map (TAG_W => MTAG_W)
    port map (
      clk => clk, rst => rst,
      in_valid => mul_valid, in_a => mul_a, in_b => mul_b, in_tag => mul_tag,
      out_valid => res_valid, out_r => res_r, out_tag => res_tag);

  rq_empty <= rq_wr = rq_rd;
  fl_empty <= fl_wr = fl_rd;

  -- ---------------------------------------------------------------- --
  -- input handshake
  -- ---------------------------------------------------------------- --
  fill_ok    <= fb_valid and not fill_pend;
  p1_take    <= p1_valid and fill_ok;
  p0_adv     <= p0_valid and (not p1_valid or p1_take);
  in_rdy     <= '1' when rst = '0' and (p0_valid = '0' or p0_adv = '1') else '0';
  in_ready   <= in_rdy;
  dummy_fill <= flushing and fill_ok and not p1_valid;

  main : process (clk)
    variable b       : bid_t;
    variable ph      : ph_t;
    variable lvl     : lvl_t;
    variable idx     : idx_t;
    variable n, c, s : node_t;
    variable i       : leaf_t;
    variable last    : std_logic;
    variable leafy   : boolean;
    variable oa, ob  : gf_t;
    variable x3, y3  : gf_t;
    variable rb      : bid_t;
    variable rkind   : ph_t;
    variable ridx    : idx_t;
    variable rlast   : std_logic;
    variable rn      : node_t;
    variable ri      : leaf_t;
    variable nxt_consumed : boolean;
    variable rq_pushed    : boolean;
    variable la      : natural range 0 to NB * W - 1;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        p0_valid <= '0'; p1_valid <= '0';
        fb_valid <= '0'; fill_cnt <= (others => '0'); fill_pend <= '0';
        idle_cnt <= (others => '0'); flushing <= '0';
        rq_wr <= (others => '0'); rq_rd <= (others => '0');
        fl <= q_identity;
        fl_wr <= to_unsigned(NB, LOG_NB + 1); fl_rd <= (others => '0');
        cur_valid <= '0'; nxt_valid <= '0';
        ra_valid <= '0'; mul_valid <= '0';
        o1_valid <= '0'; o2_valid <= '0'; out_valid <= '0';
      else
        rq_pushed := false;

        -- ============ input pipeline ============
        if p1_take = '1' then
          p1_valid <= '0';
        end if;
        if p0_adv = '1' then
          p1_valid <= '1';
          p1_x     <= p0_x;
          p1_y     <= p0_y;
          p1_tag   <= p0_tag;
          p1_hw    <= hw_sum(p0_parts);
          p0_valid <= '0';
        end if;
        if in_valid = '1' and in_rdy = '1' then
          p0_valid <= '1';
          p0_x     <= in_x;
          p0_y     <= in_y;
          p0_tag   <= in_tag;
          p0_parts <= gf_weight_parts(in_x);
        end if;

        -- ============ fill ============
        if fb_valid = '0' then
          if not fl_empty then
            fb       <= fl(to_integer(fl_rd(LOG_NB - 1 downto 0)));
            fl_rd    <= fl_rd + 1;
            fb_valid <= '1';
            fill_cnt <= (others => '0');
            idle_cnt <= (others => '0');
            flushing <= '0';
          end if;
        elsif p1_take = '1' or dummy_fill = '1' then
          la := leaf_addr(fb, fill_cnt(LOG_W - 1 downto 0));
          if p1_take = '1' then
            lx(la)   <= p1_x;
            ly(la)   <= p1_y;
            lj(la)   <= std_logic_vector(p1_hw(3 downto 1));
            ltag(la) <= p1_tag;
            lval(la) <= '1';
          else
            lx(la)   <= DUMMY_X;
            ly(la)   <= (others => '0');
            lj(la)   <= (others => '0');
            lval(la) <= '0';
          end if;
          idle_cnt <= (others => '0');
          if fill_cnt = W - 1 then
            -- batch complete: first forward level is the parents of the leaves
            b_ph(to_integer(fb))  <= to_unsigned(PH_FWD, KIND_W);
            b_lvl(to_integer(fb)) <= to_unsigned(LOG_W - 1, LVL_W);
            fill_pend <= '1';
            pend_b    <= fb;
            fb_valid  <= '0';
            flushing  <= '0';
          else
            fill_cnt <= fill_cnt + 1;
          end if;
        elsif fill_cnt /= 0 and fill_pend = '0' then
          if idle_cnt = FLUSH_CLK then
            flushing <= '1';
          else
            idle_cnt <= idle_cnt + 1;
          end if;
        end if;

        -- ============ retire ============
        o1_valid <= '0';
        if res_valid = '1' then
          rb    := unsigned(res_tag(MTAG_W - 1 downto KIND_W + IDX_W + 1));
          rkind := unsigned(res_tag(KIND_W + IDX_W downto IDX_W + 1));
          ridx  := unsigned(res_tag(IDX_W downto 1));
          rlast := res_tag(0);
          rn    := ridx(LOG_W downto 0);
          ri    := ridx(LOG_W - 1 downto 0);
          case to_integer(rkind) is
            when PH_FWD =>
              tr(tree_addr(rb, rn)) <= res_r;
              if rlast = '1' then
                if b_lvl(to_integer(rb)) = 0 then
                  b_ph(to_integer(rb))  <= to_unsigned(PH_INV, KIND_W);
                  b_lvl(to_integer(rb)) <= (others => '0');
                else
                  b_lvl(to_integer(rb)) <= b_lvl(to_integer(rb)) - 1;
                end if;
              end if;
            when PH_INV =>
              -- t(1) holds the root, then beta_2, then 1/root; t(0) the accumulator
              case to_integer(ridx(2 downto 0)) is
                when 0      => tr(tree_addr(rb, to_unsigned(1, LOG_W + 1))) <= res_r;
                when 7      => tr(tree_addr(rb, to_unsigned(1, LOG_W + 1))) <= gf_frob(res_r, 1);
                when others => tr(tree_addr(rb, to_unsigned(0, LOG_W + 1))) <= res_r;
              end case;
              if ridx(2 downto 0) = 7 then
                b_ph(to_integer(rb))  <= to_unsigned(PH_BWD, KIND_W);
                b_lvl(to_integer(rb)) <= (others => '0');
              else
                b_lvl(to_integer(rb)) <= resize(ridx(2 downto 0) + 1, LVL_W);
              end if;
            when PH_BWD =>
              tr(tree_addr(rb, rn)) <= res_r;
              if rlast = '1' then
                if b_lvl(to_integer(rb)) = LOG_W - 1 then
                  b_ph(to_integer(rb)) <= to_unsigned(PH_LAM, KIND_W);
                else
                  b_lvl(to_integer(rb)) <= b_lvl(to_integer(rb)) + 1;
                end if;
              end if;
            when PH_LAM =>
              tr(tree_addr(rb, ('1' & ri))) <= res_r;           -- leaf W+i := lam
              if rlast = '1' then
                b_ph(to_integer(rb)) <= to_unsigned(PH_FIN, KIND_W);
              end if;
            when others =>
              la := leaf_addr(rb, ri);
              y3 := res_r xor lz(la) xor ly(la);
              o1_valid <= lval(la);
              o1_x     <= lz(la);
              o1_y     <= y3;
              o1_tag   <= ltag(la);
              if rlast = '1' then
                fl(to_integer(fl_wr(LOG_NB - 1 downto 0))) <= rb;
                fl_wr <= fl_wr + 1;
              end if;
          end case;
          if rlast = '1' and to_integer(rkind) /= PH_FIN then
            rq(to_integer(rq_wr(LOG_NB - 1 downto 0))) <= rb;
            rq_wr <= rq_wr + 1;
            rq_pushed := true;
          end if;
        end if;

        -- a completed fill enters the queue on a clock no retire is using it
        if fill_pend = '1' and not rq_pushed then
          rq(to_integer(rq_wr(LOG_NB - 1 downto 0))) <= pend_b;
          rq_wr <= rq_wr + 1;
          fill_pend <= '0';
        end if;

        -- ============ burst engine ============
        nxt_consumed := false;
        if cur_valid = '1' then
          if cur_idx = cur_end then
            if nxt_valid = '1' then
              cur_b   <= nxt_b;  cur_ph <= nxt_ph;  cur_lvl <= nxt_lvl;
              cur_idx <= (others => '0');
              cur_end <= burst_end(nxt_ph, nxt_lvl);
              nxt_consumed := true;
            else
              cur_valid <= '0';
            end if;
          else
            cur_idx <= cur_idx + 1;
          end if;
        elsif nxt_valid = '1' then
          cur_valid <= '1';
          cur_b   <= nxt_b;  cur_ph <= nxt_ph;  cur_lvl <= nxt_lvl;
          cur_idx <= (others => '0');
          cur_end <= burst_end(nxt_ph, nxt_lvl);
          nxt_consumed := true;
        end if;

        if (nxt_valid = '0' or nxt_consumed) then
          if not rq_empty then
            b := rq(to_integer(rq_rd(LOG_NB - 1 downto 0)));
            rq_rd     <= rq_rd + 1;
            nxt_valid <= '1';
            nxt_b     <= b;
            nxt_ph    <= b_ph(to_integer(b));
            nxt_lvl   <= b_lvl(to_integer(b));
          else
            nxt_valid <= '0';
          end if;
        end if;

        -- ============ stage A: address and read ============
        ra_valid <= cur_valid;
        if cur_valid = '1' then
          b   := cur_b;  ph := cur_ph;  lvl := cur_lvl;  idx := cur_idx;
          last := '0';
          if idx = cur_end then last := '1'; end if;
          leafy := false;
          ra_k  <= lvl(2 downto 0);
          case to_integer(ph) is
            when PH_FWD =>
              n := to_unsigned(2 ** to_integer(lvl), LOG_W + 1) + resize(idx, LOG_W + 1);
              if to_integer(lvl) = LOG_W - 1 then
                -- children are leaves 2 idx and 2 idx + 1: d is formed in stage B
                leafy := true;
                i := shift_left(resize(idx, LOG_W), 1);
                ra_a  <= lx(leaf_addr(b, i));
                ra_ja <= lj(leaf_addr(b, i));
                ra_b  <= lx(leaf_addr(b, i or to_unsigned(1, LOG_W)));
                ra_jb <= lj(leaf_addr(b, i or to_unsigned(1, LOG_W)));
              else
                ra_a <= tr(tree_addr(b, n(LOG_W - 1 downto 0) & '0'));
                ra_b <= tr(tree_addr(b, n(LOG_W - 1 downto 0) & '1'));
              end if;
              ra_tag <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(resize(n, IDX_W)) & last;
            when PH_INV =>
              -- one multiply per burst; the chain step k is the batch level
              case to_integer(lvl(2 downto 0)) is
                when 0 | 1 =>
                  ra_a <= tr(tree_addr(b, to_unsigned(1, LOG_W + 1)));
                  ra_b <= tr(tree_addr(b, to_unsigned(1, LOG_W + 1)));
                when 7 =>
                  ra_a <= tr(tree_addr(b, to_unsigned(1, LOG_W + 1)));
                  ra_b <= tr(tree_addr(b, to_unsigned(0, LOG_W + 1)));
                when others =>
                  ra_a <= tr(tree_addr(b, to_unsigned(0, LOG_W + 1)));
                  ra_b <= tr(tree_addr(b, to_unsigned(0, LOG_W + 1)));
              end case;
              ra_tag <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(resize(lvl(2 downto 0), IDX_W)) & '1';
            when PH_BWD =>
              n := to_unsigned(2 ** to_integer(lvl), LOG_W + 1) + resize(idx(IDX_W - 1 downto 1), LOG_W + 1);
              c := n(LOG_W - 1 downto 0) & idx(0);
              s := n(LOG_W - 1 downto 0) & (not idx(0));
              ra_a <= tr(tree_addr(b, n));
              if to_integer(lvl) = LOG_W - 1 then
                leafy := true;
                ra_b  <= lx(leaf_addr(b, s(LOG_W - 1 downto 0)));
                ra_jb <= lj(leaf_addr(b, s(LOG_W - 1 downto 0)));
              else
                ra_b <= tr(tree_addr(b, s));
              end if;
              ra_tag <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(resize(c, IDX_W)) & last;
            when PH_LAM =>
              i := resize(idx, LOG_W);
              ra_a  <= ly(leaf_addr(b, i));
              ra_ja <= lj(leaf_addr(b, i));
              ra_b  <= tr(tree_addr(b, '1' & i));
              ra_c  <= lx(leaf_addr(b, i));
              ra_zaddr <= leaf_addr(b, i);
              ra_tag <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(resize(i, IDX_W)) & last;
            when others =>
              i := resize(idx, LOG_W);
              ra_a  <= tr(tree_addr(b, '1' & i));
              ra_b  <= lx(leaf_addr(b, i));
              ra_c  <= lz(leaf_addr(b, i));
              ra_zaddr <= leaf_addr(b, i);
              ra_tag <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(resize(i, IDX_W)) & last;
          end case;
          ra_ph <= ph;
          if leafy then ra_leafy <= '1'; else ra_leafy <= '0'; end if;
        end if;

        -- ============ stage B: form operands ============
        mul_valid <= ra_valid;
        mul_tag   <= ra_tag;
        if ra_valid = '1' then
          oa := ra_a;
          ob := ra_b;
          case to_integer(ra_ph) is
            when PH_FWD =>
              if ra_leafy = '1' then
                oa := ra_a xor gf_sigma_j(ra_a, ra_ja);
                ob := ra_b xor gf_sigma_j(ra_b, ra_jb);
              end if;
            when PH_INV =>
              ob := inv_frob(ra_b, ra_k);
            when PH_BWD =>
              if ra_leafy = '1' then
                ob := ra_b xor gf_sigma_j(ra_b, ra_jb);
              end if;
            when PH_LAM =>
              oa := ra_a xor gf_sigma_j(ra_a, ra_ja);
              lz(ra_zaddr) <= ra_c xor gf_sigma_j(ra_c, ra_ja);   -- d, for x3 later
            when others =>
              -- ra_a = lam, ra_b = x, ra_c = d
              x3 := gf_frob(ra_a, 1) xor ra_a xor ra_c;
              ob := ra_b xor x3;
              lz(ra_zaddr) <= x3;
          end case;
          mul_a <= oa;
          mul_b <= ob;
        end if;

        -- ============ output: weight and DP test over two clocks ============
        o2_valid <= o1_valid;
        o2_x     <= o1_x;
        o2_y     <= o1_y;
        o2_tag   <= o1_tag;
        o2_parts <= gf_weight_parts(o1_x);

        out_valid <= o2_valid;
        out_x     <= o2_x;
        out_y     <= o2_y;
        out_tag   <= o2_tag;
        out_hw    <= hw_sum(o2_parts);
        if to_integer(hw_sum(o2_parts)) <= DP_WEIGHT then
          out_dp <= '1';
        else
          out_dp <= '0';
        end if;
      end if;
    end if;
  end process;

end architecture;
