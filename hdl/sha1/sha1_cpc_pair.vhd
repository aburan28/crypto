-- sha1_cpc_pair.vhd
-- Chosen-prefix collision SCAFFOLD: paired SHA-1 pipelines with
-- differential-path conditioning hooks.
--
-- This is a STRUCTURAL SCAFFOLD, not a working CPC engine. The actual
-- attack requires:
--   * A pre-computed differential path (Wang/De Canniere/Mendel/Stevens
--     style) -- a set of XOR differences delta_W[t] on message words and
--     delta_State[t] on the working state that hold with high probability
--     across rounds.
--   * A list of "sufficient conditions" on intermediate state bits that,
--     when satisfied, force the path to actually hold rather than just
--     being probabilistically expected.
--   * Message-modification rules that turn a partial near-collision
--     trial into another by perturbing free message-word bits.
-- These artefacts come from a separate path-search tool (Stevens'
-- hashclash, or the Leurent-Peyrin "SHA-1 is a Shambles" toolchain).
-- Here we expose them as generics / constants so the rest of the engine
-- can be designed and timed without committing to a specific path.
--
-- The block runs two independent sha1_pipe instances on (M, M_prime),
-- with M_prime = M xor MSG_DELTA. After the round of interest, it
-- checks that the state XOR difference matches EXPECTED_DELTA on the
-- bits specified by CONDITION_MASK. Mismatches set out_reject (a host
-- can then either drop the trial or have an early-reject network of
-- intermediate checks running in parallel).

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.sha1_pkg.all;

entity sha1_cpc_pair is
  generic (
    TAG_WIDTH       : positive := 64;

    -- Differential path: XOR mask applied to M to produce M'.
    MSG_DELTA       : block_t  := (others => (others => '0'));

    -- Target state-XOR after compression and which bits we care about.
    EXPECTED_DELTA  : state_t  := (others => (others => '0'));
    CONDITION_MASK  : state_t  := (others => (others => '0'))   -- '1' bits = must match
  );
  port (
    clk         : in  std_logic;
    rst         : in  std_logic;

    in_valid    : in  std_logic;
    in_m        : in  block_t;
    in_h        : in  state_t;
    in_tag      : in  std_logic_vector(TAG_WIDTH - 1 downto 0);

    out_valid   : out std_logic;
    out_h       : out state_t;            -- digest from M
    out_h_prime : out state_t;            -- digest from M'
    out_diff    : out state_t;            -- H xor H'
    out_accept  : out std_logic;          -- '1' = diff matches path
    out_tag     : out std_logic_vector(TAG_WIDTH - 1 downto 0)
  );
end entity;

architecture rtl of sha1_cpc_pair is

  signal m_prime    : block_t;

  signal v0, v1     : std_logic;
  signal h0, h1     : state_t;
  signal t0, t1     : std_logic_vector(TAG_WIDTH - 1 downto 0);

  signal diff_r     : state_t;
  signal h0_r, h1_r : state_t;
  signal tag_r      : std_logic_vector(TAG_WIDTH - 1 downto 0);
  signal accept_r   : std_logic;
  signal valid_r    : std_logic := '0';

begin

  -- Apply the message delta.
  gen_delta : for i in 0 to 15 generate
    m_prime(i) <= in_m(i) xor MSG_DELTA(i);
  end generate;

  pipe_a : entity work.sha1_pipe
    generic map (TAG_WIDTH => TAG_WIDTH)
    port map (
      clk       => clk,
      rst       => rst,
      in_valid  => in_valid,
      in_m      => in_m,
      in_h      => in_h,
      in_tag    => in_tag,
      out_valid => v0,
      out_h     => h0,
      out_tag   => t0
    );

  pipe_b : entity work.sha1_pipe
    generic map (TAG_WIDTH => TAG_WIDTH)
    port map (
      clk       => clk,
      rst       => rst,
      in_valid  => in_valid,
      in_m      => m_prime,
      in_h      => in_h,
      in_tag    => in_tag,
      out_valid => v1,
      out_h     => h1,
      out_tag   => t1
    );

  -- Final-stage compare: digest XOR vs the path's expected difference,
  -- masked by CONDITION_MASK (so the host can ignore "don't care" bits).
  -- diff/accept/h/tag are all registered together so they stay in sync
  -- with valid_r (which lags v0/v1 by one cycle).
  process(clk)
    variable acc : std_logic;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        accept_r <= '0';
        valid_r  <= '0';
      else
        acc := '1';
        for i in 0 to 4 loop
          diff_r(i) <= h0(i) xor h1(i);
          if ((h0(i) xor h1(i) xor EXPECTED_DELTA(i)) and CONDITION_MASK(i)) /= x"00000000" then
            acc := '0';
          end if;
        end loop;
        h0_r     <= h0;
        h1_r     <= h1;
        tag_r    <= t0;
        accept_r <= acc;
        valid_r  <= v0 and v1;
      end if;
    end if;
  end process;

  out_valid   <= valid_r;
  out_h       <= h0_r;
  out_h_prime <= h1_r;
  out_diff    <= diff_r;
  out_accept  <= accept_r;
  out_tag     <= tag_r;

end architecture;
