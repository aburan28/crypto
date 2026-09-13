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
-- retire MUL_LATENCY clocks later; this needs MUL_LATENCY >= 2.  A burst
-- that reads what the previous one wrote issues at least two clocks after
-- that write has landed, so no word is read on the clock it is written and
-- the memories' collision behaviour never matters (checked in simulation).
--
-- Memories.  The field-element arrays are memory blocks: the address
-- register, the synchronous read, the output register, three clocks from
-- the burst engine to the data; one write port and one read port each, so
-- an array read through two addresses is kept twice.  The first version had
-- them in distributed RAM, 3k LUTs per engine; that fitted, but in a full
-- device the write address and enable of a 131 x 256 LUTRAM fan out to a
-- thousand LUTs spread over the SLICEMs of a whole region, and those nets
-- were what failed timing in every engine of a 48-engine build.  A memory
-- block's address pins fan out to two.
--
--   leaf    (x, y, j, valid)  block RAM, written at fill, read at issue
--   d_a, d_b                  block RAM, written at fill, read at issue
--   t_a, t_b                  UltraRAM: the tree, written at retire, read
--                             at issue
--
-- The tree is the largest array (2W words per batch) and it goes in
-- UltraRAM, of which the device has 960 blocks the design otherwise leaves
-- empty, while block RAM is what bounds the number of engines: 17 tiles per
-- engine with the tree in block RAM, 13 with it in UltraRAM.  Two 72-bit
-- UltraRAMs hold a 131 x 512 copy; their depth (4096) is mostly unused, but
-- the ports are the resource, not the bits.  Both memories' write and read
-- addresses come from registers of their own so the placer can put them
-- beside the blocks, which for the UltraRAM columns are further from an
-- engine's logic than its block RAMs.
--
-- Issue reads are addressed from the burst engine's registers and land
-- three clocks later in stage A3.  Nothing is read at retire: the final multiply
-- of a walk (lam times x + x3) has y, the tag and x3 in hand when it issues,
-- and they ride a shift register beside the multiplier (SRLs, ~300 LUTs)
-- and meet the product on the way out.  An earlier version kept a second
-- leaf table (y, tag, valid) and an x3 table for the retire side to read
-- through the multiplier's look-ahead tag; that was 4 RAMB36 + 1 RAMB18 of
-- the engine's 20 + 2, and block RAM is the resource that bounds the number
-- of engines in a device.  The retire path itself is a memory write plus
-- one XOR.
--
-- Batches shorter than W -- the tail of a run, or a testbench -- would wait
-- forever for leaves that never come, so a batch that has been partly
-- filled for FLUSH_CLK clocks with nothing arriving is padded with dummy
-- leaves (weight-1 x, d /= 0) that produce no output.
--
-- Slot memory per batch of W: W (x, y, d) + 2W tree words, all 131 bits.

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

  -- clocks from the burst engine's registers to the memories' data: the
  -- address register, the synchronous read, the output register
  constant RD_LAT : natural := 3;

  function maxn (a, b : natural) return natural is
  begin
    if a > b then return a; else return b; end if;
  end function;

  -- node/leaf index in the tag: a tree node (LOG_W+1 bits) or an inversion
  -- step (3 bits)
  constant IDX_W  : natural := maxn(LOG_W + 1, 3);
  constant LVL_W  : natural := maxn(LOG_W, 3);
  constant KIND_W : natural := 3;
  -- batch, phase, level, index, last: the level rides along so the
  -- retire side computes the batch's next state from the tag alone and
  -- never reads the level array it writes (a 16:1 mux into a decrement
  -- and a write decode, one of the engine's longest paths)
  constant MTAG_W : natural := LOG_NB + KIND_W + LVL_W + IDX_W + 1;

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

  -- memory words
  -- leaf word: x, y, j, tag, valid (fields from the top)
  constant LA_W   : natural := 2 * M + 3 + TAG_W + 1;
  constant LA_X   : natural := LA_W - 1;             -- x  = (LA_X downto LA_Y + 1)
  constant LA_Y   : natural := M + 3 + TAG_W;        -- y  = (LA_Y downto LA_J + 1)
  constant LA_J   : natural := TAG_W + 3;            -- j  = (LA_J downto LA_J - 2)
  constant LA_TAG : natural := TAG_W;                -- tag = (LA_TAG downto 1), valid = (0)
  subtype la_word_t is std_logic_vector(LA_W - 1 downto 0);
  subtype laddr_t is natural range 0 to NB * W - 1;
  subtype taddr_t is natural range 0 to NB * 2 * W - 1;
  -- what the final multiply carries to its retire: x3, y, tag, valid
  constant FP_W : natural := 2 * M + TAG_W + 1;
  subtype fp_word_t is std_logic_vector(FP_W - 1 downto 0);
  -- loaded in stage A3, read when the product retires: A3 -> B -> the
  -- multiplier's MUL_LATENCY -> the retire clock
  constant FP_DEPTH : natural := MUL_LATENCY + 2;
  type fp_pipe_t is array (0 to FP_DEPTH - 1) of fp_word_t;

  type la_mem_t is array (0 to NB * W - 1) of la_word_t;
  type gf_leaf_mem_t is array (0 to NB * W - 1) of gf_t;
  type gf_tree_mem_t is array (0 to NB * 2 * W - 1) of gf_t;
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

  function leaf_addr (b : bid_t; i : unsigned) return laddr_t is
  begin
    return to_integer(b & i(LOG_W - 1 downto 0));
  end function;

  function tree_addr (b : bid_t; n : node_t) return taddr_t is
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
  -- memories (the leaf tables block RAM, the tree UltraRAM; see the header)
  -- ------------------------------------------------------------------ --
  signal m_la       : la_mem_t;
  signal m_da, m_db : gf_leaf_mem_t;
  signal m_ta, m_tb : gf_tree_mem_t;

  attribute ram_style : string;
  attribute ram_style of m_la : signal is "block";
  attribute ram_style of m_da : signal is "block";
  attribute ram_style of m_db : signal is "block";
  attribute ram_style of m_ta : signal is "ultra";
  attribute ram_style of m_tb : signal is "ultra";

  -- read side: address (combinational, then registered), latch output,
  -- output register
  signal ra_qa, ra_qb : laddr_t;
  signal ra_ta, ra_tb : taddr_t;
  signal ad_qa, ad_qb : laddr_t := 0;
  signal ad_ta, ad_tb : taddr_t := 0;
  signal rd_la : la_word_t;
  signal rd_da, rd_db, rd_ta, rd_tb : gf_t;
  signal rr_la : la_word_t;
  signal rr_da, rr_db, rr_ta, rr_tb : gf_t;

  -- the final multiply's companions, shifted every clock like the
  -- multiplier's own tag so bubbles keep them aligned
  signal fp : fp_pipe_t := (others => (others => '0'));

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
  -- the registered leaf write (address once per table)
  signal w_en      : std_logic := '0';
  signal w_la_a, w_da_a, w_db_a : laddr_t := 0;
  signal w_la_word : la_word_t := (others => '0');
  signal w_d       : gf_t := (others => '0');
  -- the registered tree write, likewise
  signal tw_en     : std_logic := '0';
  signal tw_a, tw_b : taddr_t := 0;
  signal tw_d      : gf_t := (others => '0');

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

  -- the dummy leaf's d, a constant
  constant DUMMY_D : gf_t := DUMMY_X xor gf_sigma_j(DUMMY_X, "000");

  -- stages A0 .. A2 ride beside the memory read; A3 has the operands
  type ph_pipe_t   is array (0 to RD_LAT - 1) of ph_t;
  type k_pipe_t    is array (0 to RD_LAT - 1) of unsigned(2 downto 0);
  type mtag_pipe_t is array (0 to RD_LAT - 1) of mtag_t;
  signal a_valid : std_logic_vector(0 to RD_LAT - 1) := (others => '0');
  signal a_leafy : std_logic_vector(0 to RD_LAT - 1) := (others => '0');
  signal a_ph    : ph_pipe_t := (others => (others => '0'));
  signal a_k     : k_pipe_t := (others => (others => '0'));
  signal a_tag   : mtag_pipe_t := (others => (others => '0'));

  signal ra_valid : std_logic := '0';
  signal ra_a, ra_b : gf_t := (others => '0');
  signal ra_x3    : gf_t := (others => '0');    -- FIN: lam^2 + lam + d
  signal ra_ph    : ph_t := (others => '0');
  signal ra_ja    : std_logic_vector(2 downto 0) := (others => '0');
  signal ra_k     : unsigned(2 downto 0) := (others => '0');
  signal ra_tag   : mtag_t := (others => '0');

  -- multiplier
  signal mul_valid : std_logic := '0';
  signal mul_a, mul_b : gf_t := (others => '0');
  signal mul_tag   : mtag_t := (others => '0');
  signal res_valid : std_logic;
  signal res_r     : gf_t;
  signal res_tag   : mtag_t;

  -- output: four registers so the weight of x3 has four clocks (group
  -- popcounts, sums of four groups, their sum, the compare); the 22-way
  -- sum in one clock was the engine's worst path at 3 ns
  signal o1_valid, o2_valid, o3_valid, o4_valid : std_logic := '0';
  signal o1_x, o1_y, o2_x, o2_y, o3_x, o3_y, o4_x, o4_y : gf_t := (others => '0');
  signal o1_tag, o2_tag, o3_tag, o4_tag : tag_t := (others => '0');
  signal o2_parts : hw_parts_t := (others => (others => '0'));
  signal o3_quads : hw_quads_t := (others => (others => '0'));
  signal o4_hw    : hw_t := (others => '0');

begin

  assert LOG_W >= 1 report "ec2k_batch_pipe needs at least two walks per batch" severity failure;
  assert MUL_LATENCY >= 2 report "in-place tree update needs MUL_LATENCY >= 2" severity failure;

  mul : entity work.gf131_mul
    generic map (TAG_W => MTAG_W)
    port map (
      clk => clk, rst => rst,
      in_valid => mul_valid, in_a => mul_a, in_b => mul_b, in_tag => mul_tag,
      out_valid => res_valid, out_r => res_r, out_tag => res_tag,
      ahead_valid => open, ahead_tag => open);

  rq_empty <= rq_wr = rq_rd;
  fl_empty <= fl_wr = fl_rd;

  -- ---------------------------------------------------------------- --
  -- memory read ports: the address register (one per port, so each sits
  -- beside its memory), the synchronous read, then the output register.
  -- The writes are in the main process (one writer each).
  -- ---------------------------------------------------------------- --
  mem_rd : process (clk)
  begin
    if rising_edge(clk) then
      ad_qa <= ra_qa;
      ad_qb <= ra_qb;
      ad_ta <= ra_ta;
      ad_tb <= ra_tb;
      rd_la <= m_la(ad_qa);
      rd_da <= m_da(ad_qa);
      rd_db <= m_db(ad_qb);
      rd_ta <= m_ta(ad_ta);
      rd_tb <= m_tb(ad_tb);
    end if;
  end process;

  -- in simulation: a tree word is never read, by a multiply that will use
  -- it, on the clock it is written (the design does not depend on which
  -- value a colliding read returns; the address registers do follow the
  -- idle burst engine, so only reads with a_valid count)
  -- pragma translate_off
  tree_collision : process (clk)
  begin
    if rising_edge(clk) then
      assert not (tw_en = '1' and a_valid(0) = '1' and (tw_a = ad_ta or tw_b = ad_tb))
        report "ec2k_batch_pipe: tree read/write collision" severity failure;
    end if;
  end process;
  -- pragma translate_on

  mem_oreg : process (clk)
  begin
    if rising_edge(clk) then
      rr_la <= rd_la;
      rr_da <= rd_da;
      rr_db <= rd_db;
      rr_ta <= rd_ta;
      rr_tb <= rd_tb;
    end if;
  end process;

  -- ---------------------------------------------------------------- --
  -- issue-side read addresses: one per memory port, from the burst
  -- engine's registers; the phase later picks which port feeds an operand
  -- ---------------------------------------------------------------- --
  addr : process (cur_b, cur_ph, cur_lvl, cur_idx)
    variable b       : bid_t;
    variable ph      : ph_t;
    variable lvl     : lvl_t;
    variable idx     : idx_t;
    variable n, s    : node_t;
    variable pa, pb  : node_t;
    variable i       : leaf_t;
    variable qa, qb  : leaf_t;
  begin
    b := cur_b;  ph := cur_ph;  lvl := cur_lvl;  idx := cur_idx;
    i  := resize(idx, LOG_W);
    qa := i;  qb := i;
    pa := '1' & i;  pb := '1' & i;
    case to_integer(ph) is
      when PH_FWD =>
        n  := to_unsigned(2 ** to_integer(lvl), LOG_W + 1) + resize(idx, LOG_W + 1);
        pa := n(LOG_W - 1 downto 0) & '0';
        pb := n(LOG_W - 1 downto 0) & '1';
        qa := pa(LOG_W - 1 downto 0);
        qb := pb(LOG_W - 1 downto 0);
      when PH_INV =>
        case to_integer(lvl(2 downto 0)) is
          when 0 | 1 =>
            pa := to_unsigned(1, LOG_W + 1);  pb := to_unsigned(1, LOG_W + 1);
          when 7 =>
            pa := to_unsigned(1, LOG_W + 1);  pb := to_unsigned(0, LOG_W + 1);
          when others =>
            pa := to_unsigned(0, LOG_W + 1);  pb := to_unsigned(0, LOG_W + 1);
        end case;
      when PH_BWD =>
        n  := to_unsigned(2 ** to_integer(lvl), LOG_W + 1) + resize(idx(IDX_W - 1 downto 1), LOG_W + 1);
        s  := n(LOG_W - 1 downto 0) & (not idx(0));
        pa := n;
        pb := s;
        qb := s(LOG_W - 1 downto 0);
      when others =>
        null;
    end case;
    ra_ta <= tree_addr(b, pa);
    ra_tb <= tree_addr(b, pb);
    ra_qa <= leaf_addr(b, qa);
    ra_qb <= leaf_addr(b, qb);
  end process;

  -- ---------------------------------------------------------------- --
  -- input handshake.  A stage takes a walk only when it is empty, so each
  -- stage's data enable is its own valid bit and nothing else: one walk
  -- per two clocks through here, against the one per 5.3 the step unit
  -- consumes.  Letting a stage take and pass on the same clock chained
  -- the fill's state through three stages into 300 clock enables, the
  -- worst path of the routed 64-engine image (fill_pend -> p0_x, 0.03 ns).
  -- ---------------------------------------------------------------- --
  fill_ok    <= fb_valid and not fill_pend;
  p1_take    <= p1_valid and fill_ok;
  p0_adv     <= p0_valid and not p1_valid;
  in_rdy     <= not p0_valid and not rst;
  in_ready   <= in_rdy;
  dummy_fill <= flushing and fill_ok and not p1_valid;

  main : process (clk)
    variable b       : bid_t;
    variable ph      : ph_t;
    variable lvl     : lvl_t;
    variable idx     : idx_t;
    variable n, c    : node_t;
    variable i       : leaf_t;
    variable last    : std_logic;
    variable oa, ob  : gf_t;
    variable x3, y3  : gf_t;
    variable fp_in   : fp_word_t;
    variable rb      : bid_t;
    variable rkind   : ph_t;
    variable rlvl    : lvl_t;
    variable ridx    : idx_t;
    variable rlast   : std_logic;
    variable rn      : node_t;
    variable ri      : leaf_t;
    variable nxt_consumed : boolean;
    variable rq_pushed    : boolean;
    variable j       : std_logic_vector(2 downto 0);
    variable la      : laddr_t;
  begin
    if rising_edge(clk) then
      -- The body runs every clock; the reset, applied last, overrides the
      -- control registers only.  Data registers never see it, so no reset
      -- net reaches thousands of clock enables (see ec2k_axil's stage).
      rq_pushed := false;

      -- ============ input pipeline ============
      if p1_take = '1' then
        p1_valid <= '0';
      end if;
      if p1_valid = '0' then
        p1_x     <= p0_x;
        p1_y     <= p0_y;
        p1_tag   <= p0_tag;
        p1_hw    <= hw_sum(p0_parts);
      end if;
      if p0_adv = '1' then
        p1_valid <= '1';
        p0_valid <= '0';
      end if;
      if p0_valid = '0' then
        p0_x     <= in_x;
        p0_y     <= in_y;
        p0_tag   <= in_tag;
        p0_parts <= gf_weight_parts(in_x);
        if in_valid = '1' and rst = '0' then
          p0_valid <= '1';
        end if;
      end if;

      -- ============ fill ============
      -- The leaf writes go through one register stage: the write address
      -- of each table is its own copy, so the placer can put it beside the
      -- table's block RAMs (in the routed 64-engine image the fill's batch
      -- id into a leaf table's write address was the worst path, 3.0 ns).
      w_en <= '0';
      if w_en = '1' then
        m_la(w_la_a) <= w_la_word;
        m_da(w_da_a) <= w_d;
        m_db(w_db_a) <= w_d;
      end if;
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
        w_en   <= '1';
        w_la_a <= la;  w_da_a <= la;  w_db_a <= la;
        if p1_take = '1' then
          j         := std_logic_vector(p1_hw(3 downto 1));
          w_la_word <= p1_x & p1_y & j & p1_tag & '1';
          w_d       <= p1_x xor gf_sigma_j(p1_x, j);
        else
          w_la_word <= DUMMY_X & GF_ZERO & "000" & tag_t'(others => '0') & '0';
          w_d       <= DUMMY_D;
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
      -- The tree writes go through one register stage too: the product
      -- register feeds a mux (the inversion's Frobenius) and the address
      -- decode, and the UltraRAM columns are further from an engine's
      -- logic than its block RAMs.
      tw_en <= '0';
      if tw_en = '1' then
        m_ta(tw_a) <= tw_d;
        m_tb(tw_b) <= tw_d;
      end if;
      o1_valid <= '0';
      if res_valid = '1' then
        rb    := unsigned(res_tag(MTAG_W - 1 downto KIND_W + LVL_W + IDX_W + 1));
        rkind := unsigned(res_tag(KIND_W + LVL_W + IDX_W downto LVL_W + IDX_W + 1));
        rlvl  := unsigned(res_tag(LVL_W + IDX_W downto IDX_W + 1));
        ridx  := unsigned(res_tag(IDX_W downto 1));
        rlast := res_tag(0);
        rn    := ridx(LOG_W downto 0);
        ri    := ridx(LOG_W - 1 downto 0);
        case to_integer(rkind) is
          when PH_FWD =>
            tw_en <= '1';
            tw_a  <= tree_addr(rb, rn);  tw_b <= tree_addr(rb, rn);
            tw_d  <= res_r;
            if rlast = '1' then
              if rlvl = 0 then
                b_ph(to_integer(rb))  <= to_unsigned(PH_INV, KIND_W);
                b_lvl(to_integer(rb)) <= (others => '0');
              else
                b_lvl(to_integer(rb)) <= rlvl - 1;
              end if;
            end if;
          when PH_INV =>
            -- t(1) holds the root, then beta_2, then 1/root; t(0) the accumulator
            tw_en <= '1';
            case to_integer(ridx(2 downto 0)) is
              when 0 =>
                tw_a <= tree_addr(rb, to_unsigned(1, LOG_W + 1));
                tw_b <= tree_addr(rb, to_unsigned(1, LOG_W + 1));
                tw_d <= res_r;
              when 7 =>
                tw_a <= tree_addr(rb, to_unsigned(1, LOG_W + 1));
                tw_b <= tree_addr(rb, to_unsigned(1, LOG_W + 1));
                tw_d <= gf_frob(res_r, 1);
              when others =>
                tw_a <= tree_addr(rb, to_unsigned(0, LOG_W + 1));
                tw_b <= tree_addr(rb, to_unsigned(0, LOG_W + 1));
                tw_d <= res_r;
            end case;
            if ridx(2 downto 0) = 7 then
              b_ph(to_integer(rb))  <= to_unsigned(PH_BWD, KIND_W);
              b_lvl(to_integer(rb)) <= (others => '0');
            else
              b_lvl(to_integer(rb)) <= resize(ridx(2 downto 0) + 1, LVL_W);
            end if;
          when PH_BWD =>
            tw_en <= '1';
            tw_a  <= tree_addr(rb, rn);  tw_b <= tree_addr(rb, rn);
            tw_d  <= res_r;
            if rlast = '1' then
              if rlvl = LOG_W - 1 then
                b_ph(to_integer(rb)) <= to_unsigned(PH_LAM, KIND_W);
              else
                b_lvl(to_integer(rb)) <= rlvl + 1;
              end if;
            end if;
          when PH_LAM =>
            tw_en <= '1';                                         -- leaf W+i := lam
            tw_a  <= tree_addr(rb, ('1' & ri));  tw_b <= tree_addr(rb, ('1' & ri));
            tw_d  <= res_r;
            if rlast = '1' then
              b_ph(to_integer(rb)) <= to_unsigned(PH_FIN, KIND_W);
            end if;
          when others =>
            -- x3, y, tag and valid rode beside the multiply in fp
            y3 := res_r xor fp(FP_DEPTH - 1)(FP_W - 1 downto M + TAG_W + 1)
                        xor fp(FP_DEPTH - 1)(M + TAG_W downto TAG_W + 1);
            o1_valid <= fp(FP_DEPTH - 1)(0);
            o1_x     <= fp(FP_DEPTH - 1)(FP_W - 1 downto M + TAG_W + 1);
            o1_y     <= y3;
            o1_tag   <= fp(FP_DEPTH - 1)(TAG_W downto 1);
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

      -- ============ stage A0: tag, beside the memory address ============
      a_valid(0) <= cur_valid;
      if cur_valid = '1' then
        b   := cur_b;  ph := cur_ph;  lvl := cur_lvl;  idx := cur_idx;
        last := '0';
        if idx = cur_end then last := '1'; end if;
        a_k(0) <= lvl(2 downto 0);
        i := resize(idx, LOG_W);
        case to_integer(ph) is
          when PH_FWD =>
            n := to_unsigned(2 ** to_integer(lvl), LOG_W + 1) + resize(idx, LOG_W + 1);
            a_tag(0) <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(lvl)
                        & std_logic_vector(resize(n, IDX_W)) & last;
          when PH_INV =>
            a_tag(0) <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(lvl)
                        & std_logic_vector(resize(lvl(2 downto 0), IDX_W)) & '1';
          when PH_BWD =>
            n := to_unsigned(2 ** to_integer(lvl), LOG_W + 1) + resize(idx(IDX_W - 1 downto 1), LOG_W + 1);
            c := n(LOG_W - 1 downto 0) & idx(0);
            a_tag(0) <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(lvl)
                        & std_logic_vector(resize(c, IDX_W)) & last;
          when others =>
            -- LAM: y_i (+ sigma^j in stage B) times 1/d_i from leaf W+i
            -- FIN: lam_i from leaf W+i times x_i + x3_i, with d_i for x3
            a_tag(0) <= std_logic_vector(b) & std_logic_vector(ph)
                        & std_logic_vector(lvl)
                        & std_logic_vector(resize(i, IDX_W)) & last;
        end case;
        a_ph(0)    <= ph;
        if to_integer(lvl) = LOG_W - 1 then a_leafy(0) <= '1'; else a_leafy(0) <= '0'; end if;
      end if;

      -- ============ stages A1, A2: the memories' read and output register ============
      for s in 1 to RD_LAT - 1 loop
        a_valid(s) <= a_valid(s - 1);
        a_ph(s)    <= a_ph(s - 1);
        a_k(s)     <= a_k(s - 1);
        a_tag(s)   <= a_tag(s - 1);
        a_leafy(s) <= a_leafy(s - 1);
      end loop;

      -- ============ stage A3: operands off the memories ============
      -- The FIN leaf's x3 = lam^2 + lam + d is formed here (lam from the
      -- tree, d from its table) and enters fp with y, the tag and valid;
      -- fp shifts every clock, whatever is in A3.
      x3 := gf_frob(rr_ta, 1) xor rr_ta xor rr_da;
      fp_in := x3 & rr_la(LA_Y downto LA_J + 1) & rr_la(LA_TAG downto 0);
      fp <= fp_in & fp(0 to FP_DEPTH - 2);
      ra_valid <= a_valid(RD_LAT - 1);
      if a_valid(RD_LAT - 1) = '1' then
        ra_ja    <= rr_la(LA_J downto LA_J - 2);
        ra_x3    <= x3;
        ra_k     <= a_k(RD_LAT - 1);
        ra_tag   <= a_tag(RD_LAT - 1);
        case to_integer(a_ph(RD_LAT - 1)) is
          when PH_FWD =>
            if a_leafy(RD_LAT - 1) = '1' then ra_a <= rr_da; else ra_a <= rr_ta; end if;
            if a_leafy(RD_LAT - 1) = '1' then ra_b <= rr_db; else ra_b <= rr_tb; end if;
          when PH_INV =>
            ra_a <= rr_ta;
            ra_b <= rr_tb;
          when PH_BWD =>
            ra_a <= rr_ta;
            if a_leafy(RD_LAT - 1) = '1' then ra_b <= rr_db; else ra_b <= rr_tb; end if;
          when PH_LAM =>
            ra_a <= rr_la(LA_Y downto LA_J + 1);                  -- y
            ra_b <= rr_tb;
          when others =>
            ra_a <= rr_ta;
            ra_b <= rr_la(LA_X downto LA_Y + 1);                  -- x
        end case;
        ra_ph <= a_ph(RD_LAT - 1);
      end if;

      -- ============ stage B: form operands ============
      mul_valid <= ra_valid;
      mul_tag   <= ra_tag;
      if ra_valid = '1' then
        oa := ra_a;
        ob := ra_b;
        case to_integer(ra_ph) is
          when PH_INV =>
            ob := inv_frob(ra_b, ra_k);
          when PH_LAM =>
            oa := ra_a xor gf_sigma_j(ra_a, ra_ja);
          when PH_FIN =>
            -- ra_a = lam, ra_b = x
            ob := ra_b xor ra_x3;
          when others =>
            null;
        end case;
        mul_a <= oa;
        mul_b <= ob;
      end if;

      -- ============ output: weight and DP test over four clocks ============
      o2_valid <= o1_valid;
      o2_x     <= o1_x;
      o2_y     <= o1_y;
      o2_tag   <= o1_tag;
      o2_parts <= gf_weight_parts(o1_x);

      o3_valid <= o2_valid;
      o3_x     <= o2_x;
      o3_y     <= o2_y;
      o3_tag   <= o2_tag;
      o3_quads <= hw_quads(o2_parts);

      o4_valid <= o3_valid;
      o4_x     <= o3_x;
      o4_y     <= o3_y;
      o4_tag   <= o3_tag;
      o4_hw    <= hw_sum(o3_quads);

      out_valid <= o4_valid;
      out_x     <= o4_x;
      out_y     <= o4_y;
      out_tag   <= o4_tag;
      out_hw    <= o4_hw;
      if to_integer(o4_hw) <= DP_WEIGHT then
        out_dp <= '1';
      else
        out_dp <= '0';
      end if;

      if rst = '1' then
        p0_valid <= '0'; p1_valid <= '0';
        fb_valid <= '0'; fill_cnt <= (others => '0'); fill_pend <= '0'; w_en <= '0'; tw_en <= '0';
        idle_cnt <= (others => '0'); flushing <= '0';
        rq_wr <= (others => '0'); rq_rd <= (others => '0');
        fl <= q_identity;
        fl_wr <= to_unsigned(NB, LOG_NB + 1); fl_rd <= (others => '0');
        cur_valid <= '0'; nxt_valid <= '0';
        a_valid <= (others => '0'); ra_valid <= '0'; mul_valid <= '0';
        o1_valid <= '0'; o2_valid <= '0'; o3_valid <= '0'; o4_valid <= '0'; out_valid <= '0';
      end if;
    end if;
  end process;

end architecture;
