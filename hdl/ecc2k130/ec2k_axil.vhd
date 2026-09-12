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
-- Reports flow engine -> holding register -> queue.  Each engine holds a
-- report in its own register (and its FIFO behind that) until this block
-- acknowledges it, so a slow host or a burst of distinguished points only
-- delays reports, never loses them.  The holding registers are drained
-- round-robin, one per clock, into the queue the host reads.
--
-- Clearing RUN drops every walk in flight.  The host client treats a
-- restart like the GPU client treats a checkpoint restore of a fresh run:
-- new seeds for every walk, the corpus already uploaded is unaffected.

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

  -- engines; the reset is registered once per engine so that no single
  -- flop fans out to every engine on the die
  signal rst_eng    : std_logic;
  signal rst_eng_r  : std_logic_vector(0 to NENG - 1) := (others => '1');
  signal e_ld_valid : std_logic_vector(0 to NENG - 1);
  signal e_ld_ready : std_logic_vector(0 to NENG - 1);
  signal e_dp_valid : std_logic_vector(0 to NENG - 1);
  signal e_dp_ack   : std_logic_vector(0 to NENG - 1) := (others => '0');
  signal e_dp_id    : id_arr_t(0 to NENG - 1);
  signal e_dp_steps : cnt_arr_t(0 to NENG - 1);
  signal e_dp_x, e_dp_y : gf_arr_t(0 to NENG - 1);
  signal e_step     : std_logic_vector(0 to NENG - 1);

  -- control / status
  signal run        : std_logic := '0';
  signal steps      : unsigned(63 downto 0) := (others => '0');
  signal steps_hi   : word_t := (others => '0');
  signal dps        : unsigned(31 downto 0) := (others => '0');
  -- the per-clock step count is registered before it is added into the
  -- 64-bit counter, so no path runs popcount -> 64-bit add
  signal nstep_r    : natural range 0 to NENG := 0;

  -- load
  signal ld_pend    : std_logic := '0';
  signal ld_gid     : gid_t := (others => '0');
  signal ld_x, ld_y : gf_t := (others => '0');
  signal ld_id_reg  : word_t := (others => '0');

  -- per-engine holding register between the walker and the queue; a report
  -- is taken when the register is free (or being drained this clock) and
  -- the previous ack has been seen, then acked for one clock
  signal h_valid    : std_logic_vector(0 to NENG - 1) := (others => '0');
  signal h_id       : id_arr_t(0 to NENG - 1);
  signal h_steps    : cnt_arr_t(0 to NENG - 1);
  signal h_x, h_y   : gf_arr_t(0 to NENG - 1);
  signal h_ptr      : natural range 0 to NENG - 1 := 0;

  -- distinguished-point queue
  constant QD : natural := 2 ** DP_FIFO_W;
  type gid_mem_t is array (0 to QD - 1) of gid_t;
  type cnt_mem_t is array (0 to QD - 1) of cnt_t;
  type gf_mem_t  is array (0 to QD - 1) of gf_t;
  signal q_gid   : gid_mem_t;
  signal q_steps : cnt_mem_t;
  signal q_x, q_y : gf_mem_t;
  signal q_wr, q_rd : unsigned(DP_FIFO_W downto 0) := (others => '0');
  signal q_empty, q_full : boolean;
  signal q_count : unsigned(DP_FIFO_W downto 0);

  -- AXI
  signal aw_got, w_got : std_logic := '0';
  signal aw_addr  : std_logic_vector(11 downto 2) := (others => '0');
  signal w_data   : word_t := (others => '0');
  signal bvalid   : std_logic := '0';
  signal rvalid   : std_logic := '0';
  signal rdata    : word_t := (others => '0');

begin

  rst_eng <= rst or not run;

  engines : for i in 0 to NENG - 1 generate
    rst_reg : process (clk)
    begin
      if rising_edge(clk) then
        rst_eng_r(i) <= rst_eng;
      end if;
    end process;

    eng : entity work.ec2k_walker
      generic map (ID_W => ID_W, LOG_W => LOG_W, LOG_NB => LOG_NB,
                   FLUSH_CLK => FLUSH_CLK, CNT_W => CNT_W, DP_WEIGHT => DP_WEIGHT)
      port map (
        clk => clk, rst => rst_eng_r(i),
        ld_valid => e_ld_valid(i), ld_ready => e_ld_ready(i),
        ld_id => ld_gid(ID_W - 1 downto 0), ld_x => ld_x, ld_y => ld_y,
        dp_valid => e_dp_valid(i), dp_ack => e_dp_ack(i),
        dp_id => e_dp_id(i), dp_steps => e_dp_steps(i),
        dp_x => e_dp_x(i), dp_y => e_dp_y(i), step_pulse => e_step(i));

    e_ld_valid(i) <= '1' when ld_pend = '1' and run = '1' and eng_of(ld_gid) = i else '0';
  end generate;

  q_empty <= q_wr = q_rd;
  q_full  <= q_wr(DP_FIFO_W) /= q_rd(DP_FIFO_W)
             and q_wr(DP_FIFO_W - 1 downto 0) = q_rd(DP_FIFO_W - 1 downto 0);
  q_count <= q_wr - q_rd;

  s_awready <= not aw_got;
  s_wready  <= not w_got;
  s_bvalid  <= bvalid;
  s_bresp   <= "00";
  s_arready <= not rvalid;
  s_rvalid  <= rvalid;
  s_rdata   <= rdata;
  s_rresp   <= "00";

  main : process (clk)
    variable do_write, do_read : boolean;
    variable waddr, raddr : natural range 0 to 1023;
    variable clear, pop, push, drain : boolean;
    variable nstep : natural range 0 to NENG;
    variable g     : gid_t;
  begin
    if rising_edge(clk) then
      clear := false;
      pop   := false;

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
      do_write := aw_got = '1' and w_got = '1' and (bvalid = '0' or s_bready = '1');
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
            ld_pend <= '1';
          when 46 =>                               -- DP_POP
            pop := not q_empty;
          when others => null;
        end case;
      end if;

      -- AXI read channel
      if rvalid = '1' and s_rready = '1' then
        rvalid <= '0';
      end if;
      do_read := s_arvalid = '1' and rvalid = '0';
      if do_read then
        raddr := to_integer(unsigned(s_araddr(11 downto 2)));
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
          when 32 => rdata <= std_logic_vector(resize(q_gid(to_integer(q_rd(DP_FIFO_W - 1 downto 0))), 32));
          when 33 => rdata <= std_logic_vector(resize(q_steps(to_integer(q_rd(DP_FIFO_W - 1 downto 0))), 64)(31 downto 0));
          when 34 => rdata <= std_logic_vector(resize(q_steps(to_integer(q_rd(DP_FIFO_W - 1 downto 0))), 64)(63 downto 32));
          when 36 to 40 => rdata <= word_of(q_x(to_integer(q_rd(DP_FIFO_W - 1 downto 0))), raddr - 36);
          when 41 to 45 => rdata <= word_of(q_y(to_integer(q_rd(DP_FIFO_W - 1 downto 0))), raddr - 41);
          when others => null;
        end case;
      end if;

      -- load handshake with the selected engine
      if ld_pend = '1' and run = '1' and e_ld_ready(eng_of(ld_gid)) = '1' then
        ld_pend <= '0';
      end if;

      -- drain one holding register per clock into the queue
      drain := h_valid(h_ptr) = '1' and not q_full;
      push  := drain;
      if drain then
        g := (others => '0');
        g(ID_W - 1 downto 0) := h_id(h_ptr);
        if ENG_W > 0 then
          g(GID_W - 1 downto ID_W) := to_unsigned(h_ptr, ENG_W);
        end if;
        q_gid(to_integer(q_wr(DP_FIFO_W - 1 downto 0)))   <= g;
        q_steps(to_integer(q_wr(DP_FIFO_W - 1 downto 0))) <= h_steps(h_ptr);
        q_x(to_integer(q_wr(DP_FIFO_W - 1 downto 0)))     <= h_x(h_ptr);
        q_y(to_integer(q_wr(DP_FIFO_W - 1 downto 0)))     <= h_y(h_ptr);
        h_valid(h_ptr) <= '0';
      end if;
      if h_ptr = NENG - 1 then
        h_ptr <= 0;
      else
        h_ptr <= h_ptr + 1;
      end if;

      -- take reports into free holding registers and ack them
      for i in 0 to NENG - 1 loop
        e_dp_ack(i) <= '0';
        if e_dp_valid(i) = '1' and e_dp_ack(i) = '0'
           and (h_valid(i) = '0' or (drain and h_ptr = i)) then
          h_valid(i)  <= '1';
          h_id(i)     <= e_dp_id(i);
          h_steps(i)  <= e_dp_steps(i);
          h_x(i)      <= e_dp_x(i);
          h_y(i)      <= e_dp_y(i);
          e_dp_ack(i) <= '1';
        end if;
      end loop;

      -- queue pointers and counters (steps lags by two clocks)
      if push then
        q_wr <= q_wr + 1;
        dps  <= dps + 1;
      end if;
      if pop then
        q_rd <= q_rd + 1;
      end if;
      nstep := 0;
      for i in 0 to NENG - 1 loop
        if e_step(i) = '1' then
          nstep := nstep + 1;
        end if;
      end loop;
      nstep_r <= nstep;
      steps   <= steps + nstep_r;

      if clear then
        steps    <= (others => '0');
        dps      <= (others => '0');
        nstep_r  <= 0;
        q_wr     <= (others => '0');
        q_rd     <= (others => '0');
        h_valid  <= (others => '0');
      end if;

      if rst = '1' then
        run      <= '0';
        ld_pend  <= '0';
        steps    <= (others => '0');
        dps      <= (others => '0');
        nstep_r  <= 0;
        e_dp_ack <= (others => '0');
        q_wr     <= (others => '0');
        q_rd     <= (others => '0');
        h_valid  <= (others => '0');
        h_ptr    <= 0;
        aw_got   <= '0';
        w_got    <= '0';
        bvalid   <= '0';
        rvalid   <= '0';
      end if;
    end if;
  end process;

end architecture;
