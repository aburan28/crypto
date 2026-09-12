-- ec2k_walker.vhd
-- The rho sequencer around one ec2k_batch_pipe: NWALK independent walks,
-- each stepped until it lands on a distinguished point, which is reported
-- to the host with its walk id and step count.  The host then loads a fresh
-- start point into that id.  Nothing else crosses the boundary, which is the
-- same division of labour as the GPU client: the walk kernel reports
-- (seed, endpoint) pairs and the host owns restarts, corpus and collisions.
--
-- Walks live in a ready FIFO of (id, x, y) rather than a memory: a walk is
-- either inside the step unit or waiting in the FIFO, and a completed step
-- goes straight back into the FIFO unless it is distinguished.  A load from
-- the host enters the same FIFO; it is accepted only on clocks when no step
-- completes, so the FIFO and the per-walk step-count memory each have one
-- write port (both are distributed RAM).
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
    -- distinguished point: walk dp_id reached (dp_x, dp_y) after dp_steps steps
    dp_valid   : out std_logic;
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

  type gf_mem_t  is array (0 to NWALK - 1) of gf_t;
  type id_mem_t  is array (0 to NWALK - 1) of id_t;
  type cnt_mem_t is array (0 to NWALK - 1) of unsigned(CNT_W - 1 downto 0);

  -- ready FIFO
  signal f_x, f_y : gf_mem_t := (others => (others => '0'));
  signal f_id     : id_mem_t := (others => (others => '0'));
  signal f_wr, f_rd : unsigned(ID_W downto 0) := (others => '0');
  signal f_empty, f_full : boolean;

  signal cnt : cnt_mem_t := (others => (others => '0'));

  -- step unit
  signal s_in_valid, s_in_ready : std_logic;
  signal s_in_x, s_in_y : gf_t;
  signal s_in_tag       : std_logic_vector(ID_W - 1 downto 0);
  signal s_out_valid, s_out_dp : std_logic;
  signal s_out_x, s_out_y : gf_t;
  signal s_out_hw       : hw_t;
  signal s_out_tag      : std_logic_vector(ID_W - 1 downto 0);

  signal requeue : boolean;
  signal ld_rdy  : std_logic;

begin

  step : entity work.ec2k_batch_pipe
    generic map (TAG_W => ID_W, LOG_W => LOG_W, LOG_NB => LOG_NB,
                 FLUSH_CLK => FLUSH_CLK, DP_WEIGHT => DP_WEIGHT)
    port map (
      clk => clk, rst => rst,
      in_valid => s_in_valid, in_ready => s_in_ready,
      in_x => s_in_x, in_y => s_in_y, in_tag => s_in_tag,
      out_valid => s_out_valid, out_x => s_out_x, out_y => s_out_y,
      out_hw => s_out_hw, out_dp => s_out_dp, out_tag => s_out_tag);

  f_empty <= f_wr = f_rd;
  f_full  <= f_wr(ID_W) /= f_rd(ID_W) and f_wr(ID_W - 1 downto 0) = f_rd(ID_W - 1 downto 0);

  -- a completed, non-distinguished step takes the FIFO's write port; any
  -- completed step takes the counter array's write port, so a load is
  -- accepted only on clocks with no retirement (keeps cnt single-writer,
  -- hence distributed RAM instead of 2**ID_W * CNT_W flip-flops)
  requeue <= s_out_valid = '1' and s_out_dp = '0';
  ld_rdy  <= '1' when rst = '0' and s_out_valid = '0' and not f_full else '0';
  ld_ready <= ld_rdy;

  -- head of the FIFO is offered to the step unit every clock
  s_in_valid <= '0' when f_empty or rst = '1' else '1';
  s_in_x     <= f_x(to_integer(f_rd(ID_W - 1 downto 0)));
  s_in_y     <= f_y(to_integer(f_rd(ID_W - 1 downto 0)));
  s_in_tag   <= std_logic_vector(f_id(to_integer(f_rd(ID_W - 1 downto 0))));

  main : process (clk)
    variable w  : natural range 0 to NWALK - 1;
    variable n  : unsigned(CNT_W - 1 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        f_wr <= (others => '0');
        f_rd <= (others => '0');
        dp_valid   <= '0';
        step_pulse <= '0';
      else
        dp_valid   <= '0';
        step_pulse <= '0';

        -- pop: the step unit took the head
        if s_in_valid = '1' and s_in_ready = '1' then
          f_rd <= f_rd + 1;
        end if;

        -- a completed step: count it, then re-queue or report; else a host
        -- load, which starts the walk's counter.  One if/elsif chain so
        -- every array below has exactly one write port.
        if s_out_valid = '1' then
          w := to_integer(unsigned(s_out_tag));
          n := cnt(w) + 1;
          cnt(w)     <= n;
          step_pulse <= '1';
          if s_out_dp = '1' then
            dp_valid <= '1';
            dp_id    <= unsigned(s_out_tag);
            dp_steps <= n;
            dp_x     <= s_out_x;
            dp_y     <= s_out_y;
          else
            f_x(to_integer(f_wr(ID_W - 1 downto 0)))  <= s_out_x;
            f_y(to_integer(f_wr(ID_W - 1 downto 0)))  <= s_out_y;
            f_id(to_integer(f_wr(ID_W - 1 downto 0))) <= unsigned(s_out_tag);
            f_wr <= f_wr + 1;
          end if;
        elsif ld_valid = '1' and ld_rdy = '1' then
          w := to_integer(ld_id);
          cnt(w) <= (others => '0');
          f_x(to_integer(f_wr(ID_W - 1 downto 0)))  <= ld_x;
          f_y(to_integer(f_wr(ID_W - 1 downto 0)))  <= ld_y;
          f_id(to_integer(f_wr(ID_W - 1 downto 0))) <= ld_id;
          f_wr <= f_wr + 1;
        end if;
      end if;
    end if;
  end process;

end architecture;
