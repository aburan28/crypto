-- gf2_dsp_leaf.vhd
-- Schoolbook product of two N-bit polynomials over GF(2) with the AND
-- terms counted by DSP48E2 integer multipliers.  Latency DSP_LEAF_LAT (4)
-- clocks, one product per clock.
--
-- An integer product of two packed polynomials counts, in every output
-- field, the AND terms that a GF(2) product would XOR:
--
--   A = sum a_u 2^(3u),  B = sum b_v 2^(3v)   =>   AB = sum_k c_k 2^(3k),
--   c_k = #{(u, v) : u + v = k, a_u = b_v = 1}
--
-- With the coefficients three bits apart and at most five on the narrow
-- side, c_k <= 5 never carries into the next field, and bit 3k of the
-- integer product is the GF(2) coefficient for free.  The DSP48E2 takes
-- 27 x 18 signed, so a product covers eight coefficients of a against five
-- of b (24 x 15 unsigned): the leaf is a grid of GA x GB such products,
-- registered A/B, M and P inside the block (three clocks), and one XOR of
-- the grid's parities plus what the grid does not cover (the top N mod 8
-- coefficients of a against all of b, the low part of a against the top
-- N mod 5 of b -- one or two AND terms per output bit, formed in LUTs a
-- clock earlier).  A 17-bit leaf is six DSPs and ~70 LUTs against ~140
-- LUTs for the all-LUT schoolbook (README, "The multiplier").
--
-- Nothing here is Xilinx-specific in the source: the multiply is an
-- unsigned integer product that GHDL simulates and Vivado puts in a DSP
-- (use_dsp = "yes").
-- With the coefficients three bits apart and at most six on the narrow
-- side, c_k <= 6 never carries into the next field, and bit 3k of the
-- integer product is the GF(2) coefficient for free.  The DSP48E2 takes
-- 27 x 18 signed, i.e. 26 x 17 unsigned: nine coefficients of a (25 bits)
-- against six of b (16 bits) per product, 54 AND terms.  The leaf splits a
-- into GA pieces of at most nine and b into GB of at most six, as even as
-- they go, and takes every pair: a 17-bit leaf is 9 + 8 against 6 + 6 + 5,
-- six DSPs covering all 289 terms, and the LUTs only XOR the grid's
-- parities into place (~35 against ~140 for the all-LUT schoolbook,
-- README "The multiplier").  The pieces were 8 x 5 (24 x 15 bits) with
-- the top coefficient of a and the top two of b -- 49 terms -- as a LUT
-- corner XORed in beside the parities, ~70 LUTs; the wider pieces cover
-- them in the same six DSPs.
--
-- Registered A/B, M and P inside the block (three clocks), then the XOR
-- stage.  Nothing here is Xilinx-specific in the source: the multiply is
-- an unsigned integer product that GHDL simulates and Vivado puts in a
-- DSP (use_dsp = "yes").

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

entity gf2_dsp_leaf is
  generic (
    N : natural
  );
  port (
    clk : in  std_logic;
    a   : in  std_logic_vector(N - 1 downto 0);
    b   : in  std_logic_vector(N - 1 downto 0);
    r   : out std_logic_vector(2 * N - 2 downto 0)
  );
end entity;

architecture rtl of gf2_dsp_leaf is
  constant SP : natural := 3;                -- bits between coefficients
  constant NA : natural := 8;                -- coefficients of a per product (24 bits)
  constant NB : natural := 5;                -- coefficients of b per product (15 bits)
  constant GA : natural := N / NA;           -- full pieces of a in DSPs
  constant GB : natural := N / NB;           -- full pieces of b in DSPs
  constant WA : natural := SP * NA;
  constant WB : natural := SP * NB;
  constant WP : natural := WA + WB;
  constant SP     : natural := 3;                -- bits between coefficients
  constant NA_MAX : natural := 9;                -- coefficients of a per product (25 of 26 bits)
  constant NB_MAX : natural := 6;                -- coefficients of b per product (16 of 17 bits)
  constant GA     : natural := (N + NA_MAX - 1) / NA_MAX;   -- pieces of a
  constant GB     : natural := (N + NB_MAX - 1) / NB_MAX;   -- pieces of b
  constant WA     : natural := SP * (NA_MAX - 1) + 1;
  constant WB     : natural := SP * (NB_MAX - 1) + 1;
  constant WP     : natural := WA + WB;

  -- piece i of g even pieces of n coefficients: the first n mod g pieces
  -- are one longer
  function piece_len (g, i : natural) return natural is
  begin
    if i < N mod g then
      return N / g + 1;
    else
      return N / g;
    end if;
  end function;

  function piece_off (g, i : natural) return natural is
  begin
    return i * (N / g) + minimum(i, N mod g);
  end function;

  type apack_t is array (0 to GA - 1) of unsigned(WA - 1 downto 0);
  type bpack_t is array (0 to GB - 1) of unsigned(WB - 1 downto 0);
  type grid_t  is array (0 to GA - 1, 0 to GB - 1) of unsigned(WP - 1 downto 0);

  signal pa   : apack_t := (others => (others => '0'));
  signal pb   : bpack_t := (others => (others => '0'));
  signal m, p : grid_t  := (others => (others => (others => '0')));
  attribute use_dsp : string;
  attribute use_dsp of m : signal is "yes";

  -- the AND terms outside the grid, one clock behind the pack stage so
  -- they are ready with the DSP products
  signal a1, b1 : std_logic_vector(N - 1 downto 0) := (others => '0');
  signal rest   : std_logic_vector(2 * N - 2 downto 0) := (others => '0');
  signal rest2  : std_logic_vector(2 * N - 2 downto 0) := (others => '0');

  function masked (v : std_logic_vector; lo, hi : natural) return std_logic_vector is
    -- v with only bits lo .. hi-1 kept
    variable o : std_logic_vector(v'range) := (others => '0');
  begin
    for i in v'range loop
      if i >= lo and i < hi then
        o(i) := v(i);
      end if;
    end loop;
    return o;
  end function;
begin

  assert GA >= 1 and GB >= 1 report "gf2_dsp_leaf: N too small for a DSP piece" severity failure;
  assert NB <= 2 ** SP - 1 report "gf2_dsp_leaf: field would overflow" severity failure;
begin

  assert GA >= 1 and GB >= 1 report "gf2_dsp_leaf: N too small for a DSP piece" severity failure;
  assert NB_MAX <= 2 ** SP - 1 report "gf2_dsp_leaf: field would overflow" severity failure;
  assert piece_len(GA, 0) <= NA_MAX and piece_len(GB, 0) <= NB_MAX
    report "gf2_dsp_leaf: piece too wide" severity failure;

  process (clk)
    variable acc : std_logic_vector(2 * N - 2 downto 0);
    variable idx : natural;
  begin
    if rising_edge(clk) then
      -- stage 0: pack the coefficients three bits apart
      for i in 0 to GA - 1 loop
        for k in 0 to NA - 1 loop
          pa(i)(SP * k) <= a(i * NA + k);
        end loop;
      end loop;
      for j in 0 to GB - 1 loop
        for k in 0 to NB - 1 loop
          pb(j)(SP * k) <= b(j * NB + k);
        end loop;
      end loop;
      a1 <= a;
      b1 <= b;
        for k in 0 to NA_MAX - 1 loop
          if k < piece_len(GA, i) then
            pa(i)(SP * k) <= a(piece_off(GA, i) + k);
          end if;
        end loop;
      end loop;
      for j in 0 to GB - 1 loop
        for k in 0 to NB_MAX - 1 loop
          if k < piece_len(GB, j) then
            pb(j)(SP * k) <= b(piece_off(GB, j) + k);
          end if;
        end loop;
      end loop;

      -- stage 1: the products (DSP M register)
      for i in 0 to GA - 1 loop
        for j in 0 to GB - 1 loop
          m(i, j) <= pa(i) * pb(j);
        end loop;
      end loop;
      -- what the grid does not cover: a's top coefficients against all of
      -- b, and a's covered coefficients against b's top ones
      rest <= gf2_polymul(masked(a1, GA * NA, N), b1)
          xor gf2_polymul(masked(a1, 0, GA * NA), masked(b1, GB * NB, N));

      -- stage 2: DSP P register; the rest waits
      p     <= m;
      rest2 <= rest;

      -- stage 3: the parities, XORed into place
      acc := rest2;
      for i in 0 to GA - 1 loop
        for j in 0 to GB - 1 loop
          for k in 0 to NA + NB - 2 loop
            idx := i * NA + j * NB + k;

      -- stage 2: DSP P register
      p <= m;

      -- stage 3: the parities, XORed into place
      acc := (others => '0');
      for i in 0 to GA - 1 loop
        for j in 0 to GB - 1 loop
          for k in 0 to piece_len(GA, i) + piece_len(GB, j) - 2 loop
            idx := piece_off(GA, i) + piece_off(GB, j) + k;
            acc(idx) := acc(idx) xor p(i, j)(SP * k);
          end loop;
        end loop;
      end loop;
      r <= acc;
    end if;
  end process;

end architecture;
