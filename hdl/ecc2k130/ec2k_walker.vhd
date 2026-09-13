-- ec2k_walker.vhd
-- The rho sequencer around one ec2k_batch_pipe: NWALK independent walks,
-- each stepped until it lands on a distinguished point, which is reported
-- to the host with its walk id and step count.  The host then loads a fresh
-- start point into that id.  Nothing else crosses the boundary, which is the
-- same division of labour as the GPU client: the walk kernel reports
-- (seed, endpoint) pairs and the host owns restarts, corpus and collisions.
--
-- Walks live in a ready FIFO of (id, steps, x, y, dp) rather than a
-- memory: a walk is either inside the step unit, waiting in the FIFO, or
-- sitting in the report register.  Every completed step goes back into the
-- FIFO; the distinguished ones carry the dp flag, and when one reaches the
-- head it is moved into the report register instead of the step unit,
-- where it waits for the host's acknowledgement.  So a report is never
-- lost however slowly the host drains, and however many distinguished
-- points retire on consecutive clocks: the FIFO is the backlog and it can
-- hold every walk.  A load from the host enters the same FIFO; it is
-- accepted only on clocks when no step completes, so the FIFO has one
-- write port.
--
-- The FIFO is block RAM (304 x NWALK), read synchronously two clocks ahead
-- into a four-entry register buffer whose head is what the step unit and
-- the report register see; reads are issued while buffer plus reads in
-- flight are under four, so the head is always the oldest entry and the
-- buffer never overflows.  Distributed RAM would put the write pointer on
-- a net of 1400 LUTs, which is what failed timing in a full device.
--
-- The step count travels with the walk, through the step unit's tag and
-- back into the FIFO plus one, rather than in a per-walk counter memory:
-- there is then no read-modify-write of a RAM on the retire path, and no
-- state indexed by walk id at all.
--
-- Report handshake: dp_valid holds with stable data until the clock after
-- dp_ack is seen high; the host side pulses dp_ack for one clock when it
-- has taken the report, and must not take it again while its own ack is
-- still in flight (ec2k_axil does exactly that).  Both ends are registers,
-- so the wires between an engine and the register block have a full clock.
--
-- Throughput is that of the step unit: 5 + 5/W clocks per step per
-- multiplier, provided NWALK comfortably covers the 2**(LOG_W + LOG_NB)
-- walks the step unit holds, so a full batch is always forming.  When
-- fewer walks than that are live (start-up, or a host that is slow to
-- reload) the step unit's flush completes partial batches with dummy leaves
-- after FLUSH_CLK idle clocks, at a cost in efficiency but never a stall.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

entity ec2k_walker is
  generic (
    ID_W      : natural := 8;                -- NWALK = 2**ID_W
    LOG_W     : natural := 4;                -- step unit: walks per batch
    LOG_NB    : natural := 3;                -- step unit: batches in flight
    FLUSH_CLK : natural := 32;
    CNT_W     : natural := 32;
    DP_WEIGHT : natural := DP_WEIGHT_DEFAULT
  );
  port (
    clk        : in  std_logic;
    rst        : in  std_logic;
    -- host: start (or restart) walk ld_id at (ld_x, ld_y)
    ld_valid   : in  std_logic;
    ld_ready   : out std_logic;
    ld_id      : in  unsigned(ID_W - 1 downto 0);
    ld_x       : in  gf_t;
    ld_y       : in  gf_t;
    -- distinguished point: walk dp_id reached (dp_x, dp_y) after dp_steps
    -- steps; held until the clock after dp_ack
    dp_valid   : out std_logic;
    dp_ack     : in  std_logic;
    dp_id      : out unsigned(ID_W - 1 downto 0);
    dp_steps   : out unsigned(CNT_W - 1 downto 0);
    dp_x       : out gf_t;
    dp_y       : out gf_t;
    -- one pulse per completed step, for rate monitoring
    step_pulse : out std_logic
  );
end entity;

architecture rtl of ec2k_walker is

  constant NWALK : natural := 2 ** ID_W;

  subtype id_t  is unsigned(ID_W - 1 downto 0);
  subtype cnt_t is unsigned(CNT_W - 1 downto 0);

  -- FIFO word: x, y, id, steps, dp
  constant FW : natural := 2 * M + ID_W + CNT_W + 1;
  subtype fword_t is std_logic_vector(FW - 1 downto 0);
  constant X_LO   : natural := FW - M;
  constant Y_LO   : natural := X_LO - M;
  constant ID_LO  : natural := Y_LO - ID_W;
  constant CNT_LO : natural := 1;

  type f_mem_t is array (0 to NWALK - 1) of fword_t;
  signal f_mem : f_mem_t;
  attribute ram_style : string;
  attribute ram_style of f_mem : signal is "block";

  signal f_wr, f_rd : unsigned(ID_W downto 0) := (others => '0');
  signal f_empty, f_full : boolean;
  signal rd_w, rr_w : fword_t;                -- read latch, output register
  signal rv         : std_logic_vector(0 to 1) := (others => '0');  -- reads in flight

  -- output buffer: the FIFO's head lives here
  constant OB_LOG : natural := 2;
  constant OB_N   : natural := 2 ** OB_LOG;
  type ob_mem_t is array (0 to OB_N - 1) of fword_t;
  signal ob : ob_mem_t := (others => (others => '0'));
  signal ob_wr, ob_rd : unsigned(OB_LOG downto 0) := (others => '0');
  signal ob_empty : boolean;
  signal head_w   : fword_t;
  signal head_dp  : std_logic;
  -- the buffer's write enables are one-hot registers decoded a clock
  -- ahead (from the write pointer and the reads in flight), its read
  -- pointer selects a 304-bit mux, the report register's load is 304
  -- more: each such net from one driver is a long route in a full
  -- device, so the drivers are replicated
  signal ob_we    : std_logic_vector(0 to OB_N - 1) := (others => '0');
  attribute MAX_FANOUT : string;
  attribute MAX_FANOUT of ob_we : signal is "100";
  attribute MAX_FANOUT of ob_rd : signal is "100";

  -- step unit; the tag is (id, steps so far)
  constant TAG_W : natural := ID_W + CNT_W;
  signal s_in_valid, s_in_ready : std_logic;
  signal s_in_tag       : std_logic_vector(TAG_W - 1 downto 0);
  signal s_out_valid, s_out_dp : std_logic;
  signal s_out_x, s_out_y : gf_t;
  signal s_out_hw       : hw_t;
  signal s_out_tag      : std_logic_vector(TAG_W - 1 downto 0);

  -- report register
  signal dpv     : std_logic := '0';
  signal take_dp : boolean;
  attribute MAX_FANOUT of take_dp : signal is "100";

  signal ld_rdy  : std_logic;

begin

  step : entity work.ec2k_batch_pipe
    generic map (TAG_W => TAG_W, LOG_W => LOG_W, LOG_NB => LOG_NB,
                 FLUSH_CLK => FLUSH_CLK, DP_WEIGHT => DP_WEIGHT)
    port map (
      clk => clk, rst => rst,
      in_valid => s_in_valid, in_ready => s_in_ready,
      in_x => head_w(FW - 1 downto X_LO), in_y => head_w(X_LO - 1 downto Y_LO),
      in_tag => s_in_tag,
      out_valid => s_out_valid, out_x => s_out_x, out_y => s_out_y,
      out_hw => s_out_hw, out_dp => s_out_dp, out_tag => s_out_tag);

  f_empty  <= f_wr = f_rd;
  f_full   <= f_wr(ID_W) /= f_rd(ID_W) and f_wr(ID_W - 1 downto 0) = f_rd(ID_W - 1 downto 0);
  ob_empty <= ob_wr = ob_rd;
  head_w   <= ob(to_integer(ob_rd(OB_LOG - 1 downto 0)));
  head_dp  <= head_w(0);

  -- every completed step takes the FIFO's write port, so a load is
  -- accepted only on clocks with no retirement
  ld_rdy   <= '1' when rst = '0' and s_out_valid = '0' and not f_full else '0';
  ld_ready <= ld_rdy;

  -- head of the FIFO: a walk is offered to the step unit, a distinguished
  -- point to the report register once that is free (or being freed)
  s_in_valid <= '0' when ob_empty or rst = '1' or head_dp = '1' else '1';
  s_in_tag   <= head_w(Y_LO - 1 downto ID_LO) & head_w(ID_LO - 1 downto CNT_LO);
  take_dp    <= not ob_empty and rst = '0' and head_dp = '1' and (dpv = '0' or dp_ack = '1');

  dp_valid <= dpv;

  -- the FIFO memory's read port and output register
  mem_rd : process (clk)
  begin
    if rising_edge(clk) then
      rd_w <= f_mem(to_integer(f_rd(ID_W - 1 downto 0)));
      rr_w <= rd_w;
    end if;
  end process;

  main : process (clk)
    variable wr    : natural range 0 to NWALK - 1;
    variable issue : boolean;
    variable held  : natural range 0 to OB_N + 2;
  begin
    if rising_edge(clk) then
      -- the reset, applied last, touches the pointers and flags only; the
      -- report register and the buffer are data (see ec2k_axil's stage)
      step_pulse <= '0';

      -- prefetch: keep the buffer fed, never more than it can hold
      held := to_integer(ob_wr - ob_rd);
      if rv(0) = '1' then held := held + 1; end if;
      if rv(1) = '1' then held := held + 1; end if;
      issue := not f_empty and held < OB_N;
      if issue then
        f_rd <= f_rd + 1;
      end if;
      if issue then rv(0) <= '1'; else rv(0) <= '0'; end if;
      rv(1) <= rv(0);
      -- the word read for rv(0) lands next clock at the pointer as it
      -- will then be (plus one if a word lands this clock)
      for k in 0 to OB_N - 1 loop
        ob_we(k) <= '0';
        if rv(0) = '1' then
          if rv(1) = '1' then
            if to_integer(ob_wr(OB_LOG - 1 downto 0) + 1) = k then ob_we(k) <= '1'; end if;
          else
            if to_integer(ob_wr(OB_LOG - 1 downto 0)) = k then ob_we(k) <= '1'; end if;
          end if;
        end if;
        if ob_we(k) = '1' then
          ob(k) <= rr_w;
        end if;
      end loop;
      if rv(1) = '1' then
        ob_wr <= ob_wr + 1;
      end if;

      -- pop: the step unit took the head, or the head is a report (the
      -- two exclude each other: a report is never offered to the step
      -- unit)
      if (s_in_valid = '1' and s_in_ready = '1') or take_dp then
        ob_rd <= ob_rd + 1;
      end if;
      if take_dp then
        dpv      <= '1';
        dp_id    <= unsigned(head_w(Y_LO - 1 downto ID_LO));
        dp_x     <= head_w(FW - 1 downto X_LO);
        dp_y     <= head_w(X_LO - 1 downto Y_LO);
        dp_steps <= unsigned(head_w(ID_LO - 1 downto CNT_LO));
      end if;
      if dp_ack = '1' and not take_dp then
        dpv <= '0';
      end if;

      -- a completed step: re-queue it with its count plus one, flagged
      -- if distinguished; else a host load, which starts at zero.  One
      -- if/elsif chain so the memory has exactly one write port.
      wr := to_integer(f_wr(ID_W - 1 downto 0));
      if s_out_valid = '1' then
        step_pulse <= '1';
        f_mem(wr) <= s_out_x & s_out_y & s_out_tag(TAG_W - 1 downto CNT_W)
                     & std_logic_vector(unsigned(s_out_tag(CNT_W - 1 downto 0)) + 1)
                     & s_out_dp;
        f_wr <= f_wr + 1;
      elsif ld_valid = '1' and ld_rdy = '1' then
        f_mem(wr) <= ld_x & ld_y & std_logic_vector(ld_id)
                     & std_logic_vector(to_unsigned(0, CNT_W)) & '0';
        f_wr <= f_wr + 1;
      end if;

      if rst = '1' then
        f_wr  <= (others => '0');
        f_rd  <= (others => '0');
        ob_wr <= (others => '0');
        ob_rd <= (others => '0');
        rv    <= (others => '0');
        ob_we <= (others => '0');
        dpv        <= '0';
        step_pulse <= '0';
      end if;
    end if;
  end process;

end architecture;
