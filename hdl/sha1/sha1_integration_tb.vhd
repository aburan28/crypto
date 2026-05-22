-- sha1_integration_tb.vhd
-- Smoke tests for the wrapper modules:
--   1. sha1_pipe_prefix produces the same "abc" digest as the bare core.
--   2. sha1_cpc_pair with MSG_DELTA = 0 yields zero difference and accept='1'.
--   3. sha1_cpc_pair with a single-bit delta yields nonzero difference.
-- These are correctness sanity checks, NOT a CPC verification.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.sha1_pkg.all;

entity sha1_integration_tb is
end entity;

architecture sim of sha1_integration_tb is

  signal clk : std_logic := '0';
  signal rst : std_logic := '1';
  constant CLK_PER : time := 2 ns;

  constant M_ABC : block_t := (
    0  => x"61626380",
    15 => x"00000018",
    others => (others => '0')
  );
  constant EXPECTED : state_t := (
    x"A9993E36", x"4706816A", x"BA3E2571", x"7850C26C", x"9CD0D89D"
  );

  -- sha1_pipe_prefix
  signal pp_valid_in  : std_logic := '0';
  signal pp_valid_out : std_logic;
  signal pp_h_out     : state_t;
  signal pp_tag       : std_logic_vector(63 downto 0) := x"0123456789ABCDEF";
  signal pp_tag_out   : std_logic_vector(63 downto 0);

  -- sha1_cpc_pair (zero-delta and one-bit-delta runs use separate instances)
  signal p0_valid_in  : std_logic := '0';
  signal p0_valid_out : std_logic;
  signal p0_h         : state_t;
  signal p0_h_prime   : state_t;
  signal p0_diff      : state_t;
  signal p0_accept    : std_logic;
  signal p0_tag       : std_logic_vector(63 downto 0) := x"00000000DEADBEEF";
  signal p0_tag_out   : std_logic_vector(63 downto 0);

  signal p1_valid_in  : std_logic := '0';
  signal p1_valid_out : std_logic;
  signal p1_h         : state_t;
  signal p1_h_prime   : state_t;
  signal p1_diff      : state_t;
  signal p1_accept    : std_logic;
  signal p1_tag_out   : std_logic_vector(63 downto 0);

begin

  clk <= not clk after CLK_PER / 2;

  -- Instance 1: prefix wrapper with all-variable mask (functionally identical to bare core).
  uut_prefix : entity work.sha1_pipe_prefix
    generic map (
      TAG_WIDTH => 64,
      FIXED_W   => (others => (others => '0')),
      VAR_MASK  => x"FFFF"
    )
    port map (
      clk       => clk,
      rst       => rst,
      in_valid  => pp_valid_in,
      var_w     => M_ABC,
      in_h      => H_INIT,
      in_tag    => pp_tag,
      out_valid => pp_valid_out,
      out_h     => pp_h_out,
      out_tag   => pp_tag_out
    );

  -- Instance 2: CPC pair with zero delta.
  uut_pair_zero : entity work.sha1_cpc_pair
    generic map (
      MSG_DELTA      => (others => (others => '0')),
      EXPECTED_DELTA => (others => (others => '0')),
      CONDITION_MASK => (others => x"FFFFFFFF")        -- check all 160 bits
    )
    port map (
      clk         => clk,
      rst         => rst,
      in_valid    => p0_valid_in,
      in_m        => M_ABC,
      in_h        => H_INIT,
      in_tag      => p0_tag,
      out_valid   => p0_valid_out,
      out_h       => p0_h,
      out_h_prime => p0_h_prime,
      out_diff    => p0_diff,
      out_accept  => p0_accept,
      out_tag     => p0_tag_out
    );

  -- Instance 3: CPC pair with a single-bit message delta -> nonzero output diff.
  uut_pair_one : entity work.sha1_cpc_pair
    generic map (
      MSG_DELTA      => (0 => x"00000001", others => (others => '0')),
      EXPECTED_DELTA => (others => (others => '0')),
      CONDITION_MASK => (others => x"FFFFFFFF")
    )
    port map (
      clk         => clk,
      rst         => rst,
      in_valid    => p1_valid_in,
      in_m        => M_ABC,
      in_h        => H_INIT,
      in_tag      => p0_tag,
      out_valid   => p1_valid_out,
      out_h       => p1_h,
      out_h_prime => p1_h_prime,
      out_diff    => p1_diff,
      out_accept  => p1_accept,
      out_tag     => p1_tag_out
    );

  stim : process
    variable any_nonzero : boolean;
  begin
    rst <= '1';
    wait for 4 * CLK_PER;
    rst <= '0';
    wait for CLK_PER;

    pp_valid_in <= '1';
    p0_valid_in <= '1';
    p1_valid_in <= '1';
    wait for CLK_PER;
    pp_valid_in <= '0';
    p0_valid_in <= '0';
    p1_valid_in <= '0';

    -- Wait for digests to come out.
    wait until pp_valid_out = '1' for 200 * CLK_PER;
    assert pp_valid_out = '1' report "prefix wrapper timeout" severity failure;

    for i in 0 to 4 loop
      assert pp_h_out(i) = EXPECTED(i)
        report "prefix wrapper mismatch at word " & integer'image(i)
        severity failure;
    end loop;
    assert pp_tag_out = pp_tag report "prefix tag passthrough failed" severity failure;
    report "PASS (1): sha1_pipe_prefix matches bare core";

    -- CPC zero-delta: H == H', diff == 0, accept = '1'.
    wait until p0_valid_out = '1' for 200 * CLK_PER;
    assert p0_valid_out = '1' report "pair_zero timeout" severity failure;
    for i in 0 to 4 loop
      assert p0_h(i) = p0_h_prime(i)
        report "pair_zero h != h' at word " & integer'image(i) severity failure;
      assert p0_diff(i) = x"00000000"
        report "pair_zero diff nonzero at word " & integer'image(i) severity failure;
    end loop;
    assert p0_accept = '1' report "pair_zero accept should be 1" severity failure;
    report "PASS (2): sha1_cpc_pair zero-delta -> H == H', accept = 1";

    -- CPC one-bit delta: digest difference should be nonzero (avalanche).
    assert p1_valid_out = '1' report "pair_one timeout" severity failure;
    any_nonzero := false;
    for i in 0 to 4 loop
      if p1_diff(i) /= x"00000000" then any_nonzero := true; end if;
    end loop;
    assert any_nonzero
      report "pair_one diff is all-zero (SHA-1 avalanche failure?)"
      severity failure;
    assert p1_accept = '0'
      report "pair_one accept should be 0 (delta != expected_delta=0)"
      severity failure;
    report "PASS (3): sha1_cpc_pair 1-bit delta -> nonzero diff, accept = 0";

    report "ALL INTEGRATION TESTS PASSED";
    std.env.finish;
  end process;

end architecture;
