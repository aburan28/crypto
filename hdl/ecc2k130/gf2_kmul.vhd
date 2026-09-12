-- gf2_kmul.vhd
-- Pipelined polynomial multiplier over GF(2), Karatsuba with a schoolbook
-- leaf, as a recursive entity.  N-bit operands, 2N-1 bit product, fixed
-- latency 2*LEVELS + 1 clocks, one product per clock.
--
--   a = a1 X^H + a0,  b = b1 X^H + b0,  H = ceil(N/2)
--   a b = a0 b0  +  (a0 b0 + a1 b1 + (a0+a1)(b0+b1)) X^H  +  a1 b1 X^2H
--
-- Each level is one registered pre-add stage (the two XORs, one LUT level),
-- three instances of itself on H bits, and one registered post-combine
-- stage (an XOR of at most four terms per output bit, one LUT level).  The
-- leaf at LEVELS = 0 is a schoolbook product in one stage: with a 33-bit
-- leaf that is 33 AND terms per output bit, which LUT6s absorb three at a
-- time with their XORs -- three or four LUT levels.
--
-- Over GF(2) Karatsuba is exact with no subtractions, so unequal halves
-- cost nothing: the high half is zero-extended to H bits and the product
-- bits above 2N-2 are provably zero and simply dropped.
--
-- Two levels on 131 bits: 9 products of 33 x 33 = 9801 AND gates instead
-- of 17161 for schoolbook, for about 850 extra XORs.  Deeper recursion
-- trades more XORs for fewer ANDs; where the LUT optimum sits is a
-- synthesis question (see README), which is why LEVELS is a generic.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

entity gf2_kmul is
  generic (
    N      : natural;
    LEVELS : natural
  );
  port (
    clk : in  std_logic;
    a   : in  std_logic_vector(N - 1 downto 0);
    b   : in  std_logic_vector(N - 1 downto 0);
    r   : out std_logic_vector(2 * N - 2 downto 0)
  );
end entity;

architecture rtl of gf2_kmul is
  constant H : natural := (N + 1) / 2;
begin

  leaf : if LEVELS = 0 generate
    process (clk)
    begin
      if rising_edge(clk) then
        r <= gf2_polymul(a, b);
      end if;
    end process;
  end generate;

  node : if LEVELS > 0 generate
    subtype half_t is std_logic_vector(H - 1 downto 0);
    subtype prod_t is std_logic_vector(2 * H - 2 downto 0);
    signal a0, a1, am, b0, b1, bm : half_t := (others => '0');
    signal p0, p1, p2 : prod_t;
  begin

    pre : process (clk)
      variable ah, bh : half_t;
    begin
      if rising_edge(clk) then
        ah := (others => '0');
        bh := (others => '0');
        ah(N - H - 1 downto 0) := a(N - 1 downto H);
        bh(N - H - 1 downto 0) := b(N - 1 downto H);
        a0 <= a(H - 1 downto 0);
        b0 <= b(H - 1 downto 0);
        a1 <= ah;
        b1 <= bh;
        am <= a(H - 1 downto 0) xor ah;
        bm <= b(H - 1 downto 0) xor bh;
      end if;
    end process;

    m0 : entity work.gf2_kmul
      generic map (N => H, LEVELS => LEVELS - 1)
      port map (clk => clk, a => a0, b => b0, r => p0);
    m1 : entity work.gf2_kmul
      generic map (N => H, LEVELS => LEVELS - 1)
      port map (clk => clk, a => am, b => bm, r => p1);
    m2 : entity work.gf2_kmul
      generic map (N => H, LEVELS => LEVELS - 1)
      port map (clk => clk, a => a1, b => b1, r => p2);

    post : process (clk)
      variable w   : std_logic_vector(4 * H - 2 downto 0);
      variable mid : prod_t;
    begin
      if rising_edge(clk) then
        w   := (others => '0');
        mid := p0 xor p1 xor p2;
        w(2 * H - 2 downto 0) := p0;
        w(3 * H - 2 downto H) := w(3 * H - 2 downto H) xor mid;
        w(4 * H - 2 downto 2 * H) := w(4 * H - 2 downto 2 * H) xor p2;
        r <= w(2 * N - 2 downto 0);
      end if;
    end process;

  end generate;

end architecture;
