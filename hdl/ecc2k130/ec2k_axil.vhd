-- ec2k_axil.vhd
-- AXI4-Lite register block around NENG ec2k_walker engines: the whole
-- host-visible surface of the FPGA.  On AWS F2 it hangs off the shell's OCL
-- port (BAR0 of the application PF); anywhere else it is an ordinary 32-bit
-- AXI-Lite slave.
--
-- The host sees NENG * 2**ID_W walks with flat ids gid = eng * 2**ID_W + id.
-- It loads a start point by writing LD_ID, LD_X0..4, LD_Y0..4 and LD_GO,
-- then polls STATUS.LD_BUSY clear before the next load.  Distinguished
-- points are queued in a FIFO; the host polls STATUS.DP_AVAIL, reads DP_*,
-- and writes DP_POP.  The engines produce a distinguished point every
-- 2**26 steps or so per walk, so a dozen register accesses per event costs
-- nothing next to the walk itself.
--
-- Register map (byte offsets, 32-bit, little-endian words of the 131-bit
-- field elements: word k = bits 32k+31..32k, word 4 = bits 130..128):
--
--   0x000  MAGIC     RO  0x2C13_0001
--   0x004  CTRL      RW  bit0 RUN (engines held in reset while 0)
--                        bit1 CLEAR (self-clearing: counters and DP FIFO)
--   0x008  STATUS    RO  bit0 LD_BUSY, bit1 DP_AVAIL, bit2 always 0 (was
--                        DP_OVERFLOW), bits 23:16 DP FIFO occupancy
--   0x00C  GEOM      RO  [7:0] ID_W, [11:8] LOG_W, [15:12] LOG_NB,
--                        [23:16] DP_WEIGHT, [31:24] NENG
--   0x010  STEPS_LO  RO  completed steps, all engines; reading LO latches HI
--   0x014  STEPS_HI  RO
--   0x018  DPS       RO  distinguished points queued
--   0x01C  DROPPED   RO  always 0: the engines hold their reports until the
--                        queue takes them, nothing is dropped
--   0x020  LD_ID     RW  gid of the walk to (re)start
--   0x024  LD_X0..4  RW  0x024 0x028 0x02C 0x030 0x034
--   0x038  LD_Y0..4  RW  0x038 0x03C 0x040 0x044 0x048
--   0x04C  LD_GO     WO  any write: hand (LD_ID, LD_X, LD_Y) to its engine
--   0x050  CLOCK     RO  engine clock in kHz (the CLK_KHZ generic; 0 if the
--                        image did not say)
--   0x080  DP_ID     RO  gid of the walk at the head of the queue
--   0x084  DP_STEPS_LO RO  steps that walk took from its start point
--   0x088  DP_STEPS_HI RO
--   0x090  DP_X0..4  RO  0x090 0x094 0x098 0x09C 0x0A0
--   0x0A4  DP_Y0..4  RO  0x0A4 0x0A8 0x0AC 0x0B0 0x0B4
--   0x0B8  DP_POP    WO  any write: drop the head of the queue
--
-- The engines hang off a spine: one register stage per engine, chained, so
-- that nothing crosses the die in a clock -- with 48 or 64 engines over
-- three SLRs a flat bus (load registers fanned out to every engine, a 48:1
-- mux of reports back) is what breaks timing, not the arithmetic.  Every
-- wire between stages, and between a stage and its engine, runs register
-- to register.
--
--   down   a load (gid, x, y) shifts outward one stage per clock; the stage
--          whose engine it names keeps a copy and offers it to the engine
--          until taken.  Credits (below) travel the same way.
--   up     a report slot (gid, steps, x, y) shifts inward; an engine drops
--          its report into an empty slot passing by.  Beside it ride the
--          load acknowledgement, unused credits, and a running sum of step
--          pulses (so STEPS needs no popcount over the die either).
--
-- Reports are never dropped and the chain never stalls, by credit: this
-- block issues one credit per free queue slot not already spoken for, an
-- engine with a report waiting takes the first credit that passes and
-- only then inserts, and a credit that reaches the far end unused comes
-- back up so it can be reissued.  Queue occupancy plus outstanding credits
-- never exceeds the queue depth, so every report that enters the chain has
-- a slot waiting.  An engine holds its report in its own register (and its
-- FIFO behind that) until acknowledged, so a slow host or a burst of
-- distinguished points only delays reports.
--
-- Clearing RUN drops every walk in flight; the reset sweeps down the spine
-- a stage per clock and this block ignores the chain for 2 NENG + 8 clocks
-- after RUN rises, long enough for the sweep to have cleared everything.
-- The host client treats a restart like the GPU client treats a checkpoint
-- restore of a fresh run: new seeds for every walk, the corpus already
-- uploaded is unaffected.  A load issued while RUN is clear is dropped.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

entity ec2k_axil is
  generic (
    NENG      : natural := 1;                -- walker engines
    ID_W      : natural := 8;                -- walks per engine = 2**ID_W
    LOG_W     : natural := 4;
    LOG_NB    : natural := 3;
    FLUSH_CLK : natural := 32;
    CNT_W     : natural := 32;
    DP_WEIGHT : natural := DP_WEIGHT_DEFAULT;
    DP_FIFO_W : natural := 6;                -- queue depth = 2**DP_FIFO_W
    CLK_KHZ   : natural := 0                 -- reported in CLOCK, nothing else
  );
  port (
    clk       : in  std_logic;
    rst       : in  std_logic;               -- synchronous, active high
    -- AXI4-Lite slave
    s_awaddr  : in  std_logic_vector(31 downto 0);
    s_awvalid : in  std_logic;
    s_awready : out std_logic;
    s_wdata   : in  std_logic_vector(31 downto 0);
    s_wstrb   : in  std_logic_vector(3 downto 0);
    s_wvalid  : in  std_logic;
    s_wready  : out std_logic;
    s_bresp   : out std_logic_vector(1 downto 0);
    s_bvalid  : out std_logic;
    s_bready  : in  std_logic;
    s_araddr  : in  std_logic_vector(31 downto 0);
    s_arvalid : in  std_logic;
    s_arready : out std_logic;
    s_rdata   : out std_logic_vector(31 downto 0);
    s_rresp   : out std_logic_vector(1 downto 0);
    s_rvalid  : out std_logic;
    s_rready  : in  std_logic
  );
end entity;

architecture rtl of ec2k_axil is

  constant MAGIC : std_logic_vector(31 downto 0) := x"2C130001";

  function clog2 (n : natural) return natural is
    variable r : natural := 0;
  begin
    while 2 ** r < n loop
      r := r + 1;
    end loop;
    return r;
  end function;

  constant ENG_W : natural := clog2(NENG);
  constant GID_W : natural := ID_W + ENG_W;

  subtype word_t is std_logic_vector(31 downto 0);
  subtype id_t   is unsigned(ID_W - 1 downto 0);
  subtype gid_t  is unsigned(GID_W - 1 downto 0);
  subtype cnt_t  is unsigned(CNT_W - 1 downto 0);

  type gf_arr_t  is array (natural range <>) of gf_t;
  type id_arr_t  is array (natural range <>) of id_t;
  type cnt_arr_t is array (natural range <>) of cnt_t;

  function eng_of (g : gid_t) return natural is
  begin
    if ENG_W = 0 then
      return 0;
    else
      return to_integer(g(GID_W - 1 downto ID_W));
    end if;
  end function;

  function word_of (v : gf_t; k : natural) return word_t is
    variable w : word_t := (others => '0');
  begin
    for b in 0 to 31 loop
      if 32 * k + b < M then
        w(b) := v(32 * k + b);
      end if;
    end loop;
    return w;
  end function;

  procedure put_word (signal v : inout gf_t; k : natural; w : word_t) is
  begin
    for b in 0 to 31 loop
      if 32 * k + b < M then
        v(32 * k + b) <= w(b);
      end if;
    end loop;
  end procedure;

  function mk_gid (e : natural; id : id_t) return gid_t is
    variable g : gid_t := (others => '0');
  begin
    g(ID_W - 1 downto 0) := id;
    if ENG_W > 0 then
      g(GID_W - 1 downto ID_W) := to_unsigned(e, ENG_W);
    end if;
    return g;
  end function;

  -- engines, each behind its spine stage; the reset is a shift register
  -- down the spine, one flop per engine
  signal rst_eng    : std_logic;
  signal rst_src    : std_logic_vector(0 to NENG - 1);
  signal rst_eng_r  : std_logic_vector(0 to NENG - 1) := (others => '1');
  signal e_ld_ready : std_logic_vector(0 to NENG - 1);
  signal e_dp_valid : std_logic_vector(0 to NENG - 1);
  signal e_dp_ack   : std_logic_vector(0 to NENG - 1) := (others => '0');
  signal e_dp_id    : id_arr_t(0 to NENG - 1);
  signal e_dp_steps : cnt_arr_t(0 to NENG - 1);
  signal e_dp_x, e_dp_y : gf_arr_t(0 to NENG - 1);
  signal e_step     : std_logic_vector(0 to NENG - 1);

  -- spine.  Index i of a down signal is what stage i sees (stage i drives
  -- i+1; index 0 is this block); index i of an up signal is what stage i
  -- drives (index NENG is the far end, constants).
  constant NS_W : natural := clog2(NENG + 1);
  subtype nstep_t is unsigned(NS_W - 1 downto 0);
  type gid_arr_t   is array (natural range <>) of gid_t;
  type nstep_arr_t is array (natural range <>) of nstep_t;

  signal dn_ldv, dn_cr   : std_logic_vector(0 to NENG) := (others => '0');
  signal dn_gid          : gid_arr_t(0 to NENG);
  signal dn_x, dn_y      : gf_arr_t(0 to NENG);
  signal up_valid, up_cret, up_lddone : std_logic_vector(0 to NENG) := (others => '0');
  -- up_valid_n(i) is what up_valid(i) registers next; stage i - 1 keeps
  -- its own copy of "the slot above is free" from it, so the select of
  -- its 300-bit insert mux is a register beside the mux, not the stage
  -- above's valid bit routed across the die into 300 LUTs (a 0.06 ns
  -- path in the routed 64-engine image)
  signal up_valid_n      : std_logic_vector(0 to NENG) := (others => '0');
  signal up_slot_free    : std_logic_vector(0 to NENG - 1) := (others => '1');
  signal pend, ins       : std_logic_vector(0 to NENG - 1);
  -- the select of the 300-bit insert mux: from one LUT it was the worst
  -- path of the routed 80-engine image (2.8 ns of route to 310 loads);
  -- replicated, each copy sits among its loads
  attribute MAX_FANOUT : string;
  attribute MAX_FANOUT of ins : signal is "64";
  signal up_gid          : gid_arr_t(0 to NENG);
  signal up_steps        : cnt_arr_t(0 to NENG);
  signal up_x, up_y      : gf_arr_t(0 to NENG);
  signal up_nstep        : nstep_arr_t(0 to NENG) := (others => (others => '0'));

  -- per stage: the load kept for this engine, and a credit held for a
  -- report waiting to go
  signal l_valid, hold   : std_logic_vector(0 to NENG - 1) := (others => '0');
  signal l_id            : id_arr_t(0 to NENG - 1);
  signal l_x, l_y        : gf_arr_t(0 to NENG - 1);

  -- control / status
  signal run        : std_logic := '0';
  signal steps      : unsigned(63 downto 0) := (others => '0');
  signal steps_hi   : word_t := (others => '0');
  signal dps        : unsigned(31 downto 0) := (others => '0');
  -- the per-clock step count is registered before it is added into the
  -- 64-bit counter
  signal nstep_r    : nstep_t := (others => '0');

  -- spine head: credits out, warm-up after RUN rises
  constant WARM_CLK : natural := 2 * NENG + 8;
  signal warm       : natural range 0 to WARM_CLK := 0;

  -- load
  signal ld_pend    : std_logic := '0';                 -- LD_BUSY
  signal ld_sent    : std_logic := '0';                 -- and it is on the spine
  signal ld_gid     : gid_t := (others => '0');
  signal ld_x, ld_y : gf_t := (others => '0');
  signal ld_id_reg  : word_t := (others => '0');

  -- distinguished-point queue
  constant QD : natural := 2 ** DP_FIFO_W;
  signal credits    : natural range 0 to QD := 0;      -- issued, not yet back
  -- credits + queued, kept as a counter so the issue decision is a
  -- register against a constant and not pointer subtraction, addition and
  -- compare in one clock (the block's worst path at 3 ns otherwise)
  signal occ        : natural range 0 to QD := 0;
  type gid_mem_t is array (0 to QD - 1) of gid_t;
  type cnt_mem_t is array (0 to QD - 1) of cnt_t;
  type gf_mem_t  is array (0 to QD - 1) of gf_t;
  signal q_gid   : gid_mem_t;
  signal q_steps : cnt_mem_t;
  signal q_x, q_y : gf_mem_t;
  signal q_wr, q_rd : unsigned(DP_FIFO_W downto 0) := (others => '0');
  -- the head entry, read out of the LUTRAM every clock at the read
  -- pointer, so an AXI read is a mux of registers and the RAM's address
  -- is a register with nothing in front of it.  It lags a pop by one
  -- clock, and a read decode waits out that clock (head_wait) so the
  -- head is exact whatever the master's timing.
  signal q_head_gid   : gid_t := (others => '0');
  signal q_head_steps : cnt_t := (others => '0');
  signal q_head_x, q_head_y : gf_t := (others => '0');
  signal head_wait : std_logic := '0';
  signal q_empty, q_full : boolean;
  signal q_count : unsigned(DP_FIFO_W downto 0);

  -- AXI
  signal aw_got, w_got, ar_got : std_logic := '0';
  signal aw_addr, ar_addr : std_logic_vector(11 downto 2) := (others => '0');
  -- left alone, Vivado replicates the read address once per LUT of the
  -- read mux, 1 200 flip-flops for ten bits; a fanout limit gets a
  -- handful of copies instead
  attribute MAX_FANOUT of ar_addr : signal is "64";
  signal w_data   : word_t := (others => '0');
  signal bvalid   : std_logic := '0';
  signal rvalid   : std_logic := '0';
  signal rdata    : word_t := (others => '0');

begin

  rst_eng <= rst or not run;

  -- the head of the down chain is this block's registers; the far end of
  -- the up chain is empty, except that unused credits turn round there
  dn_gid(0) <= ld_gid;
  dn_x(0)   <= ld_x;
  dn_y(0)   <= ld_y;
  up_valid(NENG)  <= '0';
  up_valid_n(NENG) <= '0';
  up_cret(NENG)   <= dn_cr(NENG);
  up_lddone(NENG) <= '0';
  up_gid(NENG)    <= (others => '0');
  up_steps(NENG)  <= (others => '0');
  up_x(NENG)      <= (others => '0');
  up_y(NENG)      <= (others => '0');
  up_nstep(NENG)  <= (others => '0');

  engines : for i in 0 to NENG - 1 generate
    r0 : if i = 0 generate
      rst_src(0) <= rst_eng;
    end generate;
    rn : if i > 0 generate
      rst_src(i) <= rst_eng_r(i - 1);
    end generate;

    eng : entity work.ec2k_walker
      generic map (ID_W => ID_W, LOG_W => LOG_W, LOG_NB => LOG_NB,
                   FLUSH_CLK => FLUSH_CLK, CNT_W => CNT_W, DP_WEIGHT => DP_WEIGHT)
      port map (
        clk => clk, rst => rst_eng_r(i),
        ld_valid => l_valid(i), ld_ready => e_ld_ready(i),
        ld_id => l_id(i), ld_x => l_x(i), ld_y => l_y(i),
        dp_valid => e_dp_valid(i), dp_ack => e_dp_ack(i),
        dp_id => e_dp_id(i), dp_steps => e_dp_steps(i),
        dp_x => e_dp_x(i), dp_y => e_dp_y(i), step_pulse => e_step(i));

    -- a report is pending from the clock the walker raises dp_valid until
    -- the clock after our ack (the walker updates on seeing it); it is
    -- inserted into the passing slot once a credit is held and the slot
    -- is free
    pend(i) <= e_dp_valid(i) and not e_dp_ack(i);
    ins(i)  <= hold(i) and pend(i) and up_slot_free(i);
    up_valid_n(i) <= '0' when rst_eng_r(i) = '1' else
                     '1' when ins(i) = '1' else up_valid(i + 1);

    stage : process (clk)
      variable pending, accept, insert : boolean;
    begin
      if rising_edge(clk) then
        rst_eng_r(i) <= rst_src(i);
        up_slot_free(i) <= not up_valid_n(i + 1);

        -- Data moves unconditionally and is never reset: a reset that held
        -- 1 300 data flip-flops per stage still was a net from the stage's
        -- reset register into every one of their clock enables, and the
        -- first routed 48-engine image closed by 0.02 ns on exactly that
        -- net.  Only the valid bits, the credit and the ack see the reset.
        dn_gid(i + 1) <= dn_gid(i);
        dn_x(i + 1)   <= dn_x(i);
        dn_y(i + 1)   <= dn_y(i);
        if dn_ldv(i) = '1' and eng_of(dn_gid(i)) = i then
          l_id(i) <= dn_gid(i)(ID_W - 1 downto 0);
          l_x(i)  <= dn_x(i);
          l_y(i)  <= dn_y(i);
        end if;

        pending := pend(i) = '1';
        insert  := ins(i) = '1';
        if insert then
          up_gid(i)   <= mk_gid(i, e_dp_id(i));
          up_steps(i) <= e_dp_steps(i);
          up_x(i)     <= e_dp_x(i);
          up_y(i)     <= e_dp_y(i);
        else
          up_gid(i)   <= up_gid(i + 1);
          up_steps(i) <= up_steps(i + 1);
          up_x(i)     <= up_x(i + 1);
          up_y(i)     <= up_y(i + 1);
        end if;

        -- down: loads shift through; the one for this engine is kept and
        -- offered until the walker takes it (a walker takes a load only on
        -- clocks with no step retiring, so this can be a while)
        dn_ldv(i + 1) <= dn_ldv(i);
        accept := l_valid(i) = '1' and e_ld_ready(i) = '1';
        if accept then
          l_valid(i) <= '0';
        end if;
        if dn_ldv(i) = '1' and eng_of(dn_gid(i)) = i then
          l_valid(i) <= '1';
        end if;

        -- credits: keep the first one that passes while a report waits
        if dn_cr(i) = '1' and hold(i) = '0' and pending then
          hold(i)      <= '1';
          dn_cr(i + 1) <= '0';
        else
          dn_cr(i + 1) <= dn_cr(i);
        end if;

        -- up: insert into an empty slot once a credit is held, else pass
        -- the slot on
        e_dp_ack(i) <= '0';
        up_valid(i) <= up_valid_n(i);
        if insert then
          e_dp_ack(i) <= '1';
          hold(i)     <= '0';
        end if;
        up_cret(i) <= up_cret(i + 1);
        if accept then
          up_lddone(i) <= '1';
        else
          up_lddone(i) <= up_lddone(i + 1);
        end if;
        if e_step(i) = '1' then
          up_nstep(i) <= up_nstep(i + 1) + 1;
        else
          up_nstep(i) <= up_nstep(i + 1);
        end if;

        if rst_eng_r(i) = '1' then
          dn_ldv(i + 1)  <= '0';
          dn_cr(i + 1)   <= '0';
          up_valid(i)    <= '0';
          up_cret(i)     <= '0';
          up_lddone(i)   <= '0';
          up_nstep(i)    <= (others => '0');
          l_valid(i)     <= '0';
          hold(i)        <= '0';
          e_dp_ack(i)    <= '0';
        end if;
      end if;
    end process;
  end generate;

  q_empty <= q_wr = q_rd;
  q_full  <= q_wr(DP_FIFO_W) /= q_rd(DP_FIFO_W)
             and q_wr(DP_FIFO_W - 1 downto 0) = q_rd(DP_FIFO_W - 1 downto 0);
  q_count <= q_wr - q_rd;

  s_awready <= not aw_got;
  s_wready  <= not w_got;
  s_bvalid  <= bvalid;
  s_bresp   <= "00";
  s_arready <= not (rvalid or ar_got);
  s_rvalid  <= rvalid;
  s_rdata   <= rdata;
  s_rresp   <= "00";

  main : process (clk)
    variable do_write, do_read : boolean;
    variable waddr, raddr : natural range 0 to 1023;
    variable clear, pop, push : boolean;
    variable cr, oc : natural range 0 to 2 * QD + 1;
  begin
    if rising_edge(clk) then
      clear := false;
      pop   := false;
      push  := false;

      -- spine head.  Nothing moves until the reset sweep after RUN rose
      -- has passed the far end and everything it dropped has drained.
      dn_ldv(0) <= '0';
      dn_cr(0)  <= '0';
      if run = '0' or warm /= WARM_CLK then
        -- no credits out, so the occupancy is the queue's
        cr := 0;
        oc := to_integer(q_count);
        ld_pend <= '0';
        ld_sent <= '0';
        if run = '0' then
          warm <= 0;
        else
          warm <= warm + 1;
        end if;
      else
        -- occ = credits + queued: a report moves one to the other, an
        -- issue adds one, an unused credit or a pop frees one.  The issue
        -- decision reads the counter as it stood at the clock edge; a
        -- credit returning this clock is reissued next clock, not this.
        cr := credits;
        oc := occ;
        if up_valid(0) = '1' then
          push := true;
          cr   := cr - 1;
        end if;
        if up_cret(0) = '1' then
          cr := cr - 1;
          oc := oc - 1;
        end if;
        if up_lddone(0) = '1' then
          ld_pend <= '0';
          ld_sent <= '0';
        end if;
        if ld_pend = '1' and ld_sent = '0' then
          dn_ldv(0) <= '1';
          ld_sent   <= '1';
        end if;
        if occ < QD then
          dn_cr(0) <= '1';
          cr := cr + 1;
          oc := oc + 1;
        end if;
      end if;

      -- AXI write channel: address and data arrive independently
      if s_awvalid = '1' and aw_got = '0' then
        aw_addr <= s_awaddr(11 downto 2);
        aw_got  <= '1';
      end if;
      if s_wvalid = '1' and w_got = '0' then
        w_data <= s_wdata;
        w_got  <= '1';
      end if;
      if bvalid = '1' and s_bready = '1' then
        bvalid <= '0';
      end if;
      -- a write waits for the previous response to be taken rather than
      -- overlapping it, so nothing here depends on the master's bready in
      -- the same clock (that path ran from the bridge into the queue)
      do_write := aw_got = '1' and w_got = '1' and bvalid = '0';
      if do_write then
        waddr := to_integer(unsigned(aw_addr));
        aw_got <= '0';
        w_got  <= '0';
        bvalid <= '1';
        case waddr is
          when 1 =>                                -- CTRL
            run   <= w_data(0);
            clear := w_data(1) = '1';
          when 8 =>                                -- LD_ID
            ld_id_reg <= w_data;
            ld_gid    <= resize(unsigned(w_data), GID_W);
          when 9 to 13  => put_word(ld_x, waddr - 9, w_data);
          when 14 to 18 => put_word(ld_y, waddr - 14, w_data);
          when 19 =>                               -- LD_GO
            if run = '1' and ld_pend = '0' and eng_of(ld_gid) < NENG then
              ld_pend <= '1';
            end if;
          when 46 =>                               -- DP_POP
            pop := not q_empty;
          when others => null;
        end case;
      end if;

      -- AXI read channel
      if rvalid = '1' and s_rready = '1' then
        rvalid <= '0';
      end if;
      -- the address is registered a clock before the decode, so the read
      -- mux starts from a register in this block, not from the bridge
      if s_arvalid = '1' and ar_got = '0' and rvalid = '0' then
        ar_addr <= s_araddr(11 downto 2);
        ar_got  <= '1';
      end if;
      do_read := ar_got = '1' and head_wait = '0';
      if do_read then
        raddr := to_integer(unsigned(ar_addr));
        ar_got <= '0';
        rvalid <= '1';
        rdata  <= (others => '0');
        case raddr is
          when 0 => rdata <= MAGIC;
          when 1 => rdata <= (0 => run, others => '0');
          when 2 =>
            rdata(0) <= ld_pend;
            rdata(1) <= '0' when q_empty else '1';
            rdata(2) <= '0';
            rdata(16 + DP_FIFO_W downto 16) <= std_logic_vector(q_count);
          when 3 =>
            rdata(7 downto 0)   <= std_logic_vector(to_unsigned(ID_W, 8));
            rdata(11 downto 8)  <= std_logic_vector(to_unsigned(LOG_W, 4));
            rdata(15 downto 12) <= std_logic_vector(to_unsigned(LOG_NB, 4));
            rdata(23 downto 16) <= std_logic_vector(to_unsigned(DP_WEIGHT, 8));
            rdata(31 downto 24) <= std_logic_vector(to_unsigned(NENG, 8));
          when 4 =>
            rdata    <= std_logic_vector(steps(31 downto 0));
            steps_hi <= std_logic_vector(steps(63 downto 32));
          when 5 => rdata <= steps_hi;
          when 6 => rdata <= std_logic_vector(dps);
          when 7 => rdata <= (others => '0');
          when 8 => rdata <= ld_id_reg;
          when 9 to 13  => rdata <= word_of(ld_x, raddr - 9);
          when 14 to 18 => rdata <= word_of(ld_y, raddr - 14);
          when 20 => rdata <= std_logic_vector(to_unsigned(CLK_KHZ, 32));
          when 32 => rdata <= std_logic_vector(resize(q_head_gid, 32));
          when 33 => rdata <= std_logic_vector(resize(q_head_steps, 64)(31 downto 0));
          when 34 => rdata <= std_logic_vector(resize(q_head_steps, 64)(63 downto 32));
          when 36 to 40 => rdata <= word_of(q_head_x, raddr - 36);
          when 41 to 45 => rdata <= word_of(q_head_y, raddr - 41);
          when others => null;
        end case;
      end if;

      -- a report off the spine into the queue; it brought a credit, so
      -- there is room
      if push then
        assert not q_full report "ec2k_axil: report arrived at a full queue" severity failure;
        q_gid(to_integer(q_wr(DP_FIFO_W - 1 downto 0)))   <= up_gid(0);
        q_steps(to_integer(q_wr(DP_FIFO_W - 1 downto 0))) <= up_steps(0);
        q_x(to_integer(q_wr(DP_FIFO_W - 1 downto 0)))     <= up_x(0);
        q_y(to_integer(q_wr(DP_FIFO_W - 1 downto 0)))     <= up_y(0);
        q_wr <= q_wr + 1;
        dps  <= dps + 1;
      end if;
      if pop then
        q_rd <= q_rd + 1;
        oc   := oc - 1;
      end if;

      -- steps: the spine's running sum, NENG clocks late
      nstep_r <= up_nstep(0);
      steps   <= steps + nstep_r;

      if clear then
        steps    <= (others => '0');
        dps      <= (others => '0');
        nstep_r  <= (others => '0');
        q_wr     <= (others => '0');
        q_rd     <= (others => '0');
        oc       := cr;                        -- the queue is empty, credits stay out
      end if;
      credits <= cr;
      occ     <= oc;

      head_wait    <= '0';
      if pop or push or clear then
        head_wait <= '1';
      end if;
      q_head_gid   <= q_gid(to_integer(q_rd(DP_FIFO_W - 1 downto 0)));
      q_head_steps <= q_steps(to_integer(q_rd(DP_FIFO_W - 1 downto 0)));
      q_head_x     <= q_x(to_integer(q_rd(DP_FIFO_W - 1 downto 0)));
      q_head_y     <= q_y(to_integer(q_rd(DP_FIFO_W - 1 downto 0)));

      if rst = '1' then
        run      <= '0';
        credits  <= 0;
        occ      <= 0;
        steps    <= (others => '0');
        dps      <= (others => '0');
        nstep_r  <= (others => '0');
        q_wr     <= (others => '0');
        q_rd     <= (others => '0');
        aw_got   <= '0';
        w_got    <= '0';
        ar_got   <= '0';
        head_wait <= '0';
        bvalid   <= '0';
        rvalid   <= '0';
      end if;
    end if;
  end process;

end architecture;
