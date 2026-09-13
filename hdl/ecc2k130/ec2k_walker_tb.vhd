-- ec2k_walker_tb.vhd
-- End-to-end test of the rho engine: the testbench plays host to
-- ec2k_walker, loading start points and consuming distinguished points.
--
-- Each WALK record of vectors_ecc2k130.txt is a start point of order l on
-- the challenge curve, the number of steps the client's reference model
-- took to reach a distinguished point under a deliberately loose cutoff
-- (so walks report within a few dozen steps instead of 2^30), and that
-- point.  The host loads the first NWALK records into the walker, and on
-- every report checks the walk id, step count and point against the
-- record, then loads the next record into that id, until all are done.
-- Walks interleave arbitrarily inside the engine; the per-id bookkeeping is
-- what makes the check independent of that order.
--
-- The clock count here is not a throughput figure: the population of live
-- walks dwindles as records run out, and the tail is spent in partial
-- batches waiting for the flush.  ec2k_batch_tb measures the saturated
-- rate.
--
--   ghdl -a --std=08 gf131_pkg.vhd gf131_tb_pkg.vhd gf2_dsp_leaf.vhd gf2_kmul.vhd gf131_mul.vhd ec2k_batch_pipe.vhd ec2k_walker.vhd ec2k_walker_tb.vhd
--   ghdl -e --std=08 ec2k_walker_tb
--   ghdl -r --std=08 ec2k_walker_tb

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

use work.gf131_pkg.all;
use work.gf131_tb_pkg.all;

entity ec2k_walker_tb is
  generic (
    VECTORS : string := "vectors_ecc2k130.txt";
    ID_W    : natural := 5;
    LOG_W   : natural := 3;
    LOG_NB  : natural := 2;
    -- > 0: no checking; load NWALK walks, reload every report with its own
    -- point, and print the steady-state clocks per step over RATE_CLK clocks
    RATE_CLK : natural := 0
  );
end entity;

architecture sim of ec2k_walker_tb is

  constant NWALK  : natural := 2 ** ID_W;
  constant CNT_W  : natural := 32;
  constant MAXBLK : natural := 256;

  signal clk : std_logic := '0';
  signal rst : std_logic := '1';
  signal running : boolean := true;

  signal ld_valid, ld_ready : std_logic := '0';
  signal ld_id     : unsigned(ID_W - 1 downto 0) := (others => '0');
  signal ld_x, ld_y : gf_t := (others => '0');
  signal dp_valid  : std_logic;
  signal dp_ack    : std_logic := '0';
  signal dp_id     : unsigned(ID_W - 1 downto 0);
  signal dp_steps  : unsigned(CNT_W - 1 downto 0);
  signal dp_x, dp_y : gf_t;
  signal step_pulse : std_logic;

  type gf_array_t  is array (0 to MAXBLK - 1) of gf_t;
  type nat_array_t is array (0 to MAXBLK - 1) of natural;
  signal bx0, by0, bdx, bdy : gf_array_t := (others => (others => '0'));
  signal bk     : nat_array_t := (others => 0);
  signal nblk   : natural := 0;
  signal dpw    : natural := 0;
  signal loaded : boolean := false;

  signal steps, cycles, reports, errors : natural := 0;

begin

  clk <= not clk after 5 ns when running else '0';

  dut : entity work.ec2k_walker
    generic map (ID_W => ID_W, LOG_W => LOG_W, LOG_NB => LOG_NB, FLUSH_CLK => 16,
                 CNT_W => CNT_W, DP_WEIGHT => 56)
    port map (
      clk => clk, rst => rst,
      ld_valid => ld_valid, ld_ready => ld_ready,
      ld_id => ld_id, ld_x => ld_x, ld_y => ld_y,
      dp_valid => dp_valid, dp_ack => dp_ack, dp_id => dp_id, dp_steps => dp_steps,
      dp_x => dp_x, dp_y => dp_y, step_pulse => step_pulse);

  loader : process
    file f        : text;
    variable ln   : line;
    variable st   : file_open_status;
    variable kind : string(1 to 5);
    variable h    : string(1 to HEXW);
    variable ch   : character;
    variable ok   : boolean;
    variable iv   : integer;
    variable n    : natural := 0;
  begin
    file_open(st, f, VECTORS, read_mode);
    assert st = open_ok report "cannot open " & VECTORS severity failure;
    while not endfile(f) and n < MAXBLK loop
      readline(f, ln);
      if ln'length >= 5 then
        read(ln, kind, ok);
        if ok and kind = "WALK " then
          read(ln, iv, ok); dpw <= iv;                read(ln, ch, ok);
          read(ln, h, ok);  bx0(n) <= hex_to_gf(h);  read(ln, ch, ok);
          read(ln, h, ok);  by0(n) <= hex_to_gf(h);
          read(ln, iv, ok); bk(n) <= iv;              read(ln, ch, ok);
          read(ln, h, ok);  bdx(n) <= hex_to_gf(h);  read(ln, ch, ok);
          read(ln, h, ok);  bdy(n) <= hex_to_gf(h);
          n := n + 1;
        end if;
      end if;
    end loop;
    file_close(f);
    nblk <= n;
    wait for 1 ns;
    loaded <= true;
    report "loaded " & integer'image(n) & " WALK records";
    wait;
  end process;

  host : process
    type blk_of_t is array (0 to NWALK - 1) of integer;
    variable cur      : blk_of_t := (others => -1);   -- record each id is walking
    variable pending  : blk_of_t := (others => -1);   -- record waiting to be loaded
    variable next_blk : natural := 0;
    variable ld_cur   : integer := -1;
    variable w        : natural;
    variable found    : boolean;
    variable n_loaded, steps0, cyc0 : natural := 0;
  begin
    wait until loaded;
    assert dpw = 56 report "vector cutoff " & integer'image(dpw)
      & " does not match the DP_WEIGHT the walker was built with" severity failure;

    if RATE_CLK > 0 then
      wait until rising_edge(clk);
      rst <= '0';
      -- load every id from some record, then keep the population constant
      -- by restarting each reported walk at its own point
      while n_loaded < NWALK loop
        wait until rising_edge(clk);
        if ld_valid = '1' and ld_ready = '1' then
          n_loaded := n_loaded + 1;
          ld_valid <= '0';
        end if;
        if ld_valid = '0' or ld_ready = '1' then
          if n_loaded + 1 <= NWALK then
            ld_valid <= '1';
            ld_id    <= to_unsigned(n_loaded mod NWALK, ID_W);
            ld_x     <= bx0(n_loaded mod nblk);
            ld_y     <= by0(n_loaded mod nblk);
          end if;
        end if;
        dp_ack <= '0';
        if dp_valid = '1' and dp_ack = '0' then
          dp_ack <= '1';
        end if;
      end loop;
      ld_valid <= '0';
      wait until rising_edge(clk);
      steps0 := steps;  cyc0 := cycles;
      while cycles < cyc0 + RATE_CLK loop
        wait until rising_edge(clk);
        if ld_valid = '1' and ld_ready = '1' then
          ld_valid <= '0';
        end if;
        dp_ack <= '0';
        if dp_valid = '1' and dp_ack = '0' and ld_valid = '0' then
          dp_ack   <= '1';
          ld_valid <= '1';
          ld_id    <= dp_id;
          ld_x     <= dp_x;
          ld_y     <= dp_y;
        end if;
      end loop;
      report "ec2k_walker_tb: rate " & integer'image(steps - steps0) & " steps in "
             & integer'image(cycles - cyc0) & " clk = "
             & integer'image((100 * (cycles - cyc0)) / (steps - steps0)) & "/100 clk per step ("
             & integer'image(NWALK) & " walks, batches of " & integer'image(2 ** LOG_W)
             & ", " & integer'image(2 ** LOG_NB) & " in flight)";
      running <= false;
      wait;
    end if;

    assert nblk >= NWALK report "need at least " & integer'image(NWALK) & " WALK records"
      severity failure;
    for i in 0 to NWALK - 1 loop
      pending(i) := i;
    end loop;
    next_blk := NWALK;

    wait until rising_edge(clk);
    rst <= '0';

    while reports < nblk loop
      wait until rising_edge(clk);
      -- values read here are those the walker saw at this edge

      if ld_valid = '1' and ld_ready = '1' then
        cur(to_integer(ld_id)) := ld_cur;
        ld_valid <= '0';
        ld_cur   := -1;
      end if;

      -- take a report the way ec2k_axil does: ack for one clock, and skip
      -- the clock after, when the walker is still showing the acked one
      dp_ack <= '0';
      if dp_valid = '1' and dp_ack = '0' then
        dp_ack <= '1';
        w := to_integer(dp_id);
        if cur(w) < 0 then
          report "report from idle walk " & integer'image(w) severity error;
          errors <= errors + 1;
        else
          if to_integer(dp_steps) /= bk(cur(w)) then
            report "walk " & integer'image(w) & " record " & integer'image(cur(w))
                   & ": reported after " & integer'image(to_integer(dp_steps))
                   & " steps, want " & integer'image(bk(cur(w)))
              severity error;
            errors <= errors + 1;
          end if;
          if dp_x /= bdx(cur(w)) or dp_y /= bdy(cur(w)) then
            report "walk " & integer'image(w) & " record " & integer'image(cur(w))
                   & ": point mismatch, got (" & gf_to_hex(dp_x) & ", " & gf_to_hex(dp_y)
                   & ") want (" & gf_to_hex(bdx(cur(w))) & ", " & gf_to_hex(bdy(cur(w))) & ")"
              severity error;
            errors <= errors + 1;
          end if;
          cur(w) := -1;
          if next_blk < nblk then
            pending(w) := next_blk;
            next_blk := next_blk + 1;
          end if;
        end if;
        reports <= reports + 1;
      end if;

      if ld_cur < 0 then
        found := false;
        for i in 0 to NWALK - 1 loop
          if not found and pending(i) >= 0 then
            found    := true;
            ld_cur   := pending(i);
            pending(i) := -1;
            ld_valid <= '1';
            ld_id    <= to_unsigned(i, ID_W);
            ld_x     <= bx0(ld_cur);
            ld_y     <= by0(ld_cur);
          end if;
        end loop;
      end if;

      assert cycles < 100000 report "timeout after " & integer'image(reports)
        & " reports" severity failure;
    end loop;

    wait until rising_edge(clk);
    if errors = 0 then
      report "ec2k_walker_tb: " & integer'image(reports) & " distinguished points from "
             & integer'image(steps) & " steps in " & integer'image(cycles) & " clk ("
             & integer'image(NWALK) & " walks, batches of " & integer'image(2 ** LOG_W)
             & ", " & integer'image(2 ** LOG_NB) & " in flight)";
      report "ec2k_walker_tb: PASS";
    else
      report "ec2k_walker_tb: FAIL -- " & integer'image(errors) & " errors" severity failure;
    end if;
    running <= false;
    wait;
  end process;

  monitor : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '0' then
        cycles <= cycles + 1;
        if step_pulse = '1' then
          steps <= steps + 1;
        end if;
      end if;
    end if;
  end process;

end architecture;
