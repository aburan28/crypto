-- gf131_mul.vhd
-- Pipelined GF(2^131) multiplier, one product per clock, fixed latency
-- MUL_LATENCY, operands and result in the permuted type-II normal basis.
--
--   stage 0                 latch a, b
--   stage 1                 a' = prep(a), b' = prep(b)     gamma -> c-powers
--   stages 2 .. 2+2L        a' b' over GF(2)[c]             Karatsuba, L levels
--   stage 3+2L              r = to_onb(h)                   c-powers -> gamma
--
-- There is no reduction: the back-conversion maps every c^k, k = 2..262,
-- straight to normal-basis coordinates.  Every stage is a few LUT levels
-- deep -- the widest single XOR is the 66-input column of prep -- and no
-- DSP is involved anywhere, so the clock is set by LUT-to-LUT routing and
-- nothing else.
--
-- The tag rides alongside in a shift register and is never inspected, so
-- one multiplier serves any number of independent contexts at II = 1.
-- ahead_tag is the tag of the product that will retire AHEAD clocks from
-- now, so that a consumer can start a synchronous memory read the result
-- will need and have the data the clock it arrives.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

entity gf131_mul is
  generic (
    TAG_W : natural := 8;
    AHEAD : natural := 0                     -- ahead_tag leads out_tag by this
  );
  port (
    clk       : in  std_logic;
    rst       : in  std_logic;
    in_valid  : in  std_logic;
    in_a      : in  gf_t;
    in_b      : in  gf_t;
    in_tag    : in  std_logic_vector(TAG_W - 1 downto 0);
    out_valid : out std_logic;
    out_r     : out gf_t;
    out_tag   : out std_logic_vector(TAG_W - 1 downto 0);
    ahead_valid : out std_logic;
    ahead_tag   : out std_logic_vector(TAG_W - 1 downto 0)
  );
end entity;

architecture rtl of gf131_mul is

  constant KM_LAT : natural := 2 * MUL_KARATSUBA + 1;

  subtype tag_t is std_logic_vector(TAG_W - 1 downto 0);
  type tag_arr is array (0 to MUL_LATENCY - 1) of tag_t;

  signal s0_a, s0_b : gf_t := (others => '0');
  signal pa, pb     : poly_t := (others => '0');
  signal h          : dpoly_t;

  -- valid and tag travel beside the datapath; index k is the value that
  -- entered k clocks ago
  signal vl : std_logic_vector(MUL_LATENCY - 1 downto 0) := (others => '0');
  signal tg : tag_arr := (others => (others => '0'));

begin

  km : entity work.gf2_kmul
    generic map (N => M, LEVELS => MUL_KARATSUBA)
    port map (clk => clk, a => pa, b => pb, r => h);

  datapath : process (clk)
  begin
    if rising_edge(clk) then
      s0_a  <= in_a;
      s0_b  <= in_b;
      pa    <= gf_prep(s0_a);
      pb    <= gf_prep(s0_b);
      out_r <= gf_to_onb(h);
    end if;
  end process;

  control : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        vl <= (others => '0');
      else
        vl <= vl(MUL_LATENCY - 2 downto 0) & in_valid;
        tg <= in_tag & tg(0 to MUL_LATENCY - 2);
      end if;
    end if;
  end process;

  assert AHEAD < MUL_LATENCY report "gf131_mul: AHEAD must be below the latency" severity failure;

  out_valid   <= vl(MUL_LATENCY - 1);
  out_tag     <= tg(MUL_LATENCY - 1);
  ahead_valid <= vl(MUL_LATENCY - 1 - AHEAD);
  ahead_tag   <= tg(MUL_LATENCY - 1 - AHEAD);

end architecture;
