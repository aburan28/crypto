-- gf131_mul_tb.vhd
-- Self-checking testbench for the pipelined GF(2^131) multiplier.
--
-- Drives every MUL record of vectors_ecc2k130.txt (from ecc2k_ref.py, which
-- takes its products from the ECC2K-130 client's own field model) through
-- the pipeline back-to-back, one per clock.  Each result is checked twice:
-- against the file, and against gf_mul_ref, the direct normal-basis product
-- computed in simulation, which shares nothing with the multiplier's
-- polynomial-basis route.  Tags must come back in order.
--
--   ghdl -a --std=08 gf131_pkg.vhd gf131_mul.vhd gf131_mul_tb.vhd
--   ghdl -e --std=08 gf131_mul_tb
--   ghdl -r --std=08 gf131_mul_tb

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

use work.gf131_pkg.all;
use work.gf131_tb_pkg.all;

entity gf131_mul_tb is
  generic (
    VECTORS : string := "vectors_ecc2k130.txt"
  );
end entity;

architecture sim of gf131_mul_tb is

  constant TAG_W  : natural := 10;
  constant MAXVEC : natural := 4096;

  signal clk : std_logic := '0';
  signal rst : std_logic := '1';
  signal running : boolean := true;

  signal in_valid  : std_logic := '0';
  signal in_a, in_b : gf_t := (others => '0');
  signal in_tag    : std_logic_vector(TAG_W - 1 downto 0) := (others => '0');
  signal out_valid : std_logic;
  signal out_r     : gf_t;
  signal out_tag   : std_logic_vector(TAG_W - 1 downto 0);

  type vec_array_t is array (0 to MAXVEC - 1) of gf_t;
  signal va, vb, vr : vec_array_t := (others => (others => '0'));
  signal nvec   : natural := 0;
  signal loaded : boolean := false;

  signal errors  : natural := 0;
  signal checked : natural := 0;

begin

  clk <= not clk after 5 ns when running else '0';

  dut : entity work.gf131_mul
    generic map (TAG_W => TAG_W)
    port map (
      clk => clk, rst => rst,
      in_valid => in_valid, in_a => in_a, in_b => in_b, in_tag => in_tag,
      out_valid => out_valid, out_r => out_r, out_tag => out_tag);

  loader : process
    file f        : text;
    variable ln   : line;
    variable st   : file_open_status;
    variable kind : string(1 to 4);
    variable ha, hb, hr : string(1 to HEXW);
    variable ch   : character;
    variable ok   : boolean;
    variable n    : natural := 0;
  begin
    file_open(st, f, VECTORS, read_mode);
    assert st = open_ok report "cannot open " & VECTORS severity failure;
    while not endfile(f) and n < MAXVEC loop
      readline(f, ln);
      if ln'length >= 4 then
        read(ln, kind, ok);
        if ok and kind = "MUL " then
          read(ln, ha, ok);
          read(ln, ch, ok);
          read(ln, hb, ok);
          read(ln, ch, ok);
          read(ln, hr, ok);
          va(n) <= hex_to_gf(ha);
          vb(n) <= hex_to_gf(hb);
          vr(n) <= hex_to_gf(hr);
          n := n + 1;
        end if;
      end if;
    end loop;
    file_close(f);
    nvec <= n;
    wait for 1 ns;
    loaded <= true;
    report "loaded " & integer'image(n) & " MUL vectors";
    wait;
  end process;

  driver : process
  begin
    wait until loaded;
    wait until rising_edge(clk);
    rst <= '0';
    wait until rising_edge(clk);
    for i in 0 to nvec - 1 loop
      in_valid <= '1';
      in_a     <= va(i);
      in_b     <= vb(i);
      in_tag   <= std_logic_vector(to_unsigned(i mod 2 ** TAG_W, TAG_W));
      wait until rising_edge(clk);
    end loop;
    in_valid <= '0';
    for i in 0 to MUL_LATENCY + 4 loop
      wait until rising_edge(clk);
    end loop;

    assert checked = nvec
      report "expected " & integer'image(nvec) & " results, saw "
             & integer'image(checked)
      severity failure;
    if errors = 0 then
      report "gf131_mul_tb: PASS -- " & integer'image(checked)
             & " multiplications, II=1, latency "
             & integer'image(MUL_LATENCY) & " clk";
    else
      report "gf131_mul_tb: FAIL -- " & integer'image(errors) & " mismatches"
        severity failure;
    end if;
    running <= false;
    wait;
  end process;

  checker : process (clk)
    variable idx : natural := 0;
    variable ref : gf_t;
  begin
    if rising_edge(clk) then
      if out_valid = '1' then
        idx := checked;
        if out_r /= vr(idx) then
          report "vector " & integer'image(idx) & " mismatch: got "
                 & gf_to_hex(out_r) & " want " & gf_to_hex(vr(idx))
            severity error;
          errors <= errors + 1;
        end if;
        ref := gf_mul_ref(va(idx), vb(idx));
        if ref /= vr(idx) then
          report "vector " & integer'image(idx) & ": gf_mul_ref disagrees with file: "
                 & gf_to_hex(ref) & " vs " & gf_to_hex(vr(idx))
            severity error;
          errors <= errors + 1;
        end if;
        assert out_tag = std_logic_vector(to_unsigned(idx mod 2 ** TAG_W, TAG_W))
          report "tag out of order at vector " & integer'image(idx)
          severity error;
        checked <= checked + 1;
      end if;
    end if;
  end process;

end architecture;
