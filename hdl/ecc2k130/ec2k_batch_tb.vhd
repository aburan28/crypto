-- ec2k_batch_tb.vhd
-- Self-checking testbench for the batched ECC2K-130 step unit.
--
-- Same protocol and vectors as ec2k_step_tb: every STEP record of
-- vectors_ecc2k130.txt goes in as fast as the unit takes it, and x3, y3 and
-- the weight of x3 come back by tag and are checked against the file.  The
-- vector count is chosen not to be a multiple of the batch size, so the
-- last batch is completed by the flush path with dummy leaves.
--
-- The clock count is the throughput figure: with the multiplier saturated
-- it approaches 5 + 5/W clocks per step (5W + 5 multiplies per batch of W).
-- ROUNDS > 1 streams the first NFEED vectors (a multiple of W, so every
-- batch is full) that many times over, which amortises fill and drain and
-- gives the steady-state rate; the check is by tag so repeats are fine.
--
--   ghdl -a --std=08 gf131_pkg.vhd gf131_tb_pkg.vhd gf2_dsp_leaf.vhd gf2_kmul.vhd gf131_mul.vhd ec2k_batch_pipe.vhd ec2k_batch_tb.vhd
--   ghdl -e --std=08 ec2k_batch_tb
--   ghdl -r --std=08 ec2k_batch_tb
--   ghdl -r --std=08 ec2k_batch_tb -gROUNDS=8 -gNFEED=192

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

use work.gf131_pkg.all;
use work.gf131_tb_pkg.all;

entity ec2k_batch_tb is
  generic (
    VECTORS : string := "vectors_ecc2k130.txt";
    LOG_W   : natural := 4;
    LOG_NB  : natural := 3;
    ROUNDS  : natural := 1;
    NFEED   : natural := 0                   -- 0: every STEP vector in the file
  );
end entity;

architecture sim of ec2k_batch_tb is

  constant TAG_W  : natural := 8;
  constant MAXVEC : natural := 256;

  signal clk : std_logic := '0';
  signal rst : std_logic := '1';
  signal running : boolean := true;

  signal in_valid, in_ready : std_logic := '0';
  signal in_x, in_y : gf_t := (others => '0');
  signal in_tag     : std_logic_vector(TAG_W - 1 downto 0) := (others => '0');
  signal out_valid, out_dp : std_logic;
  signal out_x, out_y : gf_t;
  signal out_hw     : hw_t;
  signal out_tag    : std_logic_vector(TAG_W - 1 downto 0);

  type vec_array_t is array (0 to MAXVEC - 1) of gf_t;
  type hw_array_t  is array (0 to MAXVEC - 1) of natural;
  signal vx, vy, vx3, vy3 : vec_array_t := (others => (others => '0'));
  signal vhw, vhw3 : hw_array_t := (others => 0);
  signal nvec   : natural := 0;
  signal loaded : boolean := false;

  signal errors  : natural := 0;
  signal checked : natural := 0;
  signal cycles  : natural := 0;
  signal timing  : boolean := false;
  signal seen    : hw_array_t := (others => 0);
  signal nfeed_s : natural := 0;
  signal total   : natural := 0;

begin

  clk <= not clk after 5 ns when running else '0';

  dut : entity work.ec2k_batch_pipe
    generic map (TAG_W => TAG_W, LOG_W => LOG_W, LOG_NB => LOG_NB,
                 FLUSH_CLK => 32, DP_WEIGHT => DP_WEIGHT_DEFAULT)
    port map (
      clk => clk, rst => rst,
      in_valid => in_valid, in_ready => in_ready,
      in_x => in_x, in_y => in_y, in_tag => in_tag,
      out_valid => out_valid, out_x => out_x, out_y => out_y,
      out_hw => out_hw, out_dp => out_dp, out_tag => out_tag);

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
    while not endfile(f) and n < MAXVEC loop
      readline(f, ln);
      if ln'length >= 5 then
        read(ln, kind, ok);
        if ok and kind = "STEP " then
          read(ln, h, ok);  vx(n)  <= hex_to_gf(h);  read(ln, ch, ok);
          read(ln, h, ok);  vy(n)  <= hex_to_gf(h);
          read(ln, iv, ok); vhw(n) <= iv;             read(ln, ch, ok);
          read(ln, h, ok);  vx3(n) <= hex_to_gf(h);  read(ln, ch, ok);
          read(ln, h, ok);  vy3(n) <= hex_to_gf(h);
          read(ln, iv, ok); vhw3(n) <= iv;
          n := n + 1;
        end if;
      end if;
    end loop;
    file_close(f);
    nvec <= n;
    wait for 1 ns;
    loaded <= true;
    report "loaded " & integer'image(n) & " STEP vectors";
    wait;
  end process;

  driver : process
    variable i, n, k : natural;
  begin
    wait until loaded;
    if NFEED = 0 or NFEED > nvec then
      n := nvec;
    else
      n := NFEED;
    end if;
    nfeed_s <= n;
    total   <= n * ROUNDS;
    wait until rising_edge(clk);
    rst <= '0';
    wait until rising_edge(clk);
    timing <= true;
    k := 0;
    while k < ROUNDS loop
      i := 0;
      while i < n loop
        in_valid <= '1';
        in_x     <= vx(i);
        in_y     <= vy(i);
        in_tag   <= std_logic_vector(to_unsigned(i, TAG_W));
        wait until rising_edge(clk);
        if in_ready = '1' then
          i := i + 1;
        end if;
      end loop;
      k := k + 1;
    end loop;
    in_valid <= '0';
    while checked < total loop
      wait until rising_edge(clk);
      assert cycles < 200 * total + 4000
        report "timeout: " & integer'image(checked) & " of "
               & integer'image(total) & " results"
        severity failure;
    end loop;
    timing <= false;
    wait until rising_edge(clk);

    if errors = 0 then
      report "ec2k_batch_tb: " & integer'image(checked) & " steps in "
             & integer'image(cycles) & " clk = "
             & integer'image((100 * cycles) / checked) & "/100 clk per step, W = "
             & integer'image(2 ** LOG_W) & ", " & integer'image(2 ** LOG_NB)
             & " batches in flight (" & integer'image(5 * 2 ** LOG_W + 5)
             & " multiplies per batch)";
      report "ec2k_batch_tb: PASS";
    else
      report "ec2k_batch_tb: FAIL -- " & integer'image(errors) & " mismatches"
        severity failure;
    end if;
    running <= false;
    wait;
  end process;

  checker : process (clk)
    variable idx : natural;
  begin
    if rising_edge(clk) then
      if timing then
        cycles <= cycles + 1;
      end if;
      if out_valid = '1' then
        idx := to_integer(unsigned(out_tag));
        if idx >= nfeed_s or seen(idx) >= ROUNDS then
          report "unexpected or duplicate result for tag " & integer'image(idx)
            severity error;
          errors <= errors + 1;
        end if;
        seen(idx) <= seen(idx) + 1;
        if out_x /= vx3(idx) or out_y /= vy3(idx) then
          report "step " & integer'image(idx) & " mismatch: got ("
                 & gf_to_hex(out_x) & ", " & gf_to_hex(out_y) & ") want ("
                 & gf_to_hex(vx3(idx)) & ", " & gf_to_hex(vy3(idx)) & ")"
            severity error;
          errors <= errors + 1;
        end if;
        if to_integer(out_hw) /= vhw3(idx) then
          report "step " & integer'image(idx) & " weight: got "
                 & integer'image(to_integer(out_hw)) & " want "
                 & integer'image(vhw3(idx))
            severity error;
          errors <= errors + 1;
        end if;
        if (out_dp = '1') /= (vhw3(idx) <= DP_WEIGHT_DEFAULT) then
          report "step " & integer'image(idx) & " dp flag wrong" severity error;
          errors <= errors + 1;
        end if;
        checked <= checked + 1;
      end if;
    end if;
  end process;

end architecture;
