-- sha1_tb.vhd
-- Self-checking testbench: drives the NIST FIPS 180-4 "abc" test vector
-- into the pipelined core and asserts the emitted digest matches
-- a9993e36 4706816a ba3e2571 7850c26c 9cd0d89d.
--
-- Run with GHDL:
--   ghdl -a --std=08 sha1_pkg.vhd sha1_pipe.vhd sha1_tb.vhd
--   ghdl -e --std=08 sha1_tb
--   ghdl -r --std=08 sha1_tb

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.sha1_pkg.all;

entity sha1_tb is
end entity;

architecture sim of sha1_tb is

  signal clk       : std_logic := '0';
  signal rst       : std_logic := '1';

  signal in_valid  : std_logic := '0';
  signal in_m      : block_t   := (others => (others => '0'));
  signal in_h      : state_t   := H_INIT;
  signal in_tag    : std_logic_vector(63 downto 0) := (others => '0');

  signal out_valid : std_logic;
  signal out_h     : state_t;
  signal out_tag   : std_logic_vector(63 downto 0);

  constant CLK_PER : time := 2 ns;  -- 500 MHz

  -- "abc" padded into a single 512-bit block.
  constant M_ABC : block_t := (
    0  => x"61626380",
    15 => x"00000018",
    others => (others => '0')
  );

  constant EXPECTED : state_t := (
    x"A9993E36", x"4706816A", x"BA3E2571", x"7850C26C", x"9CD0D89D"
  );

begin

  uut : entity work.sha1_pipe
    generic map (TAG_WIDTH => 64)
    port map (
      clk       => clk,
      rst       => rst,
      in_valid  => in_valid,
      in_m      => in_m,
      in_h      => in_h,
      in_tag    => in_tag,
      out_valid => out_valid,
      out_h     => out_h,
      out_tag   => out_tag
    );

  clk <= not clk after CLK_PER / 2;

  stim : process
  begin
    rst <= '1';
    wait for 4 * CLK_PER;
    rst <= '0';
    wait for CLK_PER;

    in_m     <= M_ABC;
    in_h     <= H_INIT;
    in_tag   <= x"00000000DEADBEEF";
    in_valid <= '1';
    wait for CLK_PER;
    in_valid <= '0';
    in_m     <= (others => (others => '0'));

    -- Latency: 1 (input) + 80 (rounds) + 1 (final add) = 82 cycles.
    wait until out_valid = '1' for 200 * CLK_PER;
    assert out_valid = '1'
      report "TIMEOUT: digest never produced"
      severity failure;

    for i in 0 to 4 loop
      assert out_h(i) = EXPECTED(i)
        report "MISMATCH at word " & integer'image(i) &
               ": got "      & to_hstring(out_h(i)) &
               ", expected " & to_hstring(EXPECTED(i))
        severity failure;
    end loop;

    assert out_tag = x"00000000DEADBEEF"
      report "tag passthrough failed"
      severity failure;

    report "PASS: SHA-1(""abc"") = a9993e36 4706816a ba3e2571 7850c26c 9cd0d89d";
    std.env.finish;
  end process;

end architecture;
