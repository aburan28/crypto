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
-- write port (all distributed RAM).
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

  subtype id_t is unsigned(ID_W - 1 downto 0);

  subtype cnt_t is unsigned(CNT_W - 1 downto 0);

  type gf_mem_t  is array (0 to NWALK - 1) of gf_t;
  type id_mem_t  is array (0 to NWALK - 1) of id_t;
  type v_mem_t   is array (0 to NWALK - 1) of std_logic;
  type cnt_mem_t is array (0 to NWALK - 1) of cnt_t;

  -- ready FIFO
  signal f_x, f_y : gf_mem_t := (others => (others => '0'));
  signal f_id     : id_mem_t := (others => (others => '0'));
  signal f_cnt    : cnt_mem_t := (others => (others => '0'));
  signal f_dp     : v_mem_t  := (others => '0');
  signal f_wr, f_rd : unsigned(ID_W downto 0) := (others => '0');
  signal f_empty, f_full : boolean;
  signal head     : natural range 0 to NWALK - 1;
  signal head_dp  : std_logic;

  -- step unit; the tag is (id, steps so far)
  constant TAG_W : natural := ID_W + CNT_W;
  signal s_in_valid, s_in_ready : std_logic;
  signal s_in_x, s_in_y : gf_t;
  signal s_in_tag       : std_logic_vector(TAG_W - 1 downto 0);
  signal s_out_valid, s_out_dp : std_logic;
  signal s_out_x, s_out_y : gf_t;
  signal s_out_hw       : hw_t;
  signal s_out_tag      : std_logic_vector(TAG_W - 1 downto 0);

  -- report register
  signal dpv     : std_logic := '0';
  signal take_dp : boolean;

  signal ld_rdy  : std_logic;

begin

  step : entity work.ec2k_batch_pipe
    generic map (TAG_W => TAG_W, LOG_W => LOG_W, LOG_NB => LOG_NB,
                 FLUSH_CLK => FLUSH_CLK, DP_WEIGHT => DP_WEIGHT)
    port map (
      clk => clk, rst => rst,
      in_valid => s_in_valid, in_ready => s_in_ready,
      in_x => s_in_x, in_y => s_in_y, in_tag => s_in_tag,
      out_valid => s_out_valid, out_x => s_out_x, out_y => s_out_y,
      out_hw => s_out_hw, out_dp => s_out_dp, out_tag => s_out_tag);

  f_empty <= f_wr = f_rd;
  f_full  <= f_wr(ID_W) /= f_rd(ID_W) and f_wr(ID_W - 1 downto 0) = f_rd(ID_W - 1 downto 0);
  head    <= to_integer(f_rd(ID_W - 1 downto 0));
  head_dp <= f_dp(head);

  -- every completed step takes the FIFO's write port, so a load is
  -- accepted only on clocks with no retirement (keeps it single-writer,
  -- hence distributed RAM)
  ld_rdy   <= '1' when rst = '0' and s_out_valid = '0' and not f_full else '0';
  ld_ready <= ld_rdy;

  -- head of the FIFO: a walk is offered to the step unit, a distinguished
  -- point to the report register once that is free (or being freed)
  s_in_valid <= '0' when f_empty or rst = '1' or head_dp = '1' else '1';
  s_in_x     <= f_x(head);
  s_in_y     <= f_y(head);
  s_in_tag   <= std_logic_vector(f_id(head)) & std_logic_vector(f_cnt(head));
  take_dp    <= not f_empty and rst = '0' and head_dp = '1' and (dpv = '0' or dp_ack = '1');

  dp_valid <= dpv;

  main : process (clk)
    variable wr : natural range 0 to NWALK - 1;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        f_wr <= (others => '0');
        f_rd <= (others => '0');
        dpv        <= '0';
        step_pulse <= '0';
      else
        step_pulse <= '0';

        -- pop: the step unit took the head, or the head is a report
        if s_in_valid = '1' and s_in_ready = '1' then
          f_rd <= f_rd + 1;
        elsif take_dp then
          f_rd     <= f_rd + 1;
          dpv      <= '1';
          dp_id    <= f_id(head);
          dp_x     <= f_x(head);
          dp_y     <= f_y(head);
          dp_steps <= f_cnt(head);
        end if;
        if dp_ack = '1' and not take_dp then
          dpv <= '0';
        end if;

        -- a completed step: re-queue it with its count plus one, flagged
        -- if distinguished; else a host load, which starts at zero.  One
        -- if/elsif chain so every array below has exactly one write port.
        wr := to_integer(f_wr(ID_W - 1 downto 0));
        if s_out_valid = '1' then
          step_pulse <= '1';
          f_x(wr)   <= s_out_x;
          f_y(wr)   <= s_out_y;
          f_id(wr)  <= unsigned(s_out_tag(TAG_W - 1 downto CNT_W));
          f_cnt(wr) <= unsigned(s_out_tag(CNT_W - 1 downto 0)) + 1;
          f_dp(wr)  <= s_out_dp;
          f_wr <= f_wr + 1;
        elsif ld_valid = '1' and ld_rdy = '1' then
          f_x(wr)   <= ld_x;
          f_y(wr)   <= ld_y;
          f_id(wr)  <= ld_id;
          f_cnt(wr) <= (others => '0');
          f_dp(wr)  <= '0';
          f_wr <= f_wr + 1;
        end if;
      end if;
    end if;
  end process;

end architecture;
