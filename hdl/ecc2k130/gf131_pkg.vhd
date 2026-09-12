-- gf131_pkg.vhd
-- GF(2^131) in the permuted type-II optimal normal basis, as used by the
-- ECC2K-130 client (ecc2k130/), and the ECC2K-130 iteration function.
-- Plain IEEE VHDL-2008; no vendor primitives.  Simulates under GHDL and
-- reads into any synthesis flow.  Reference part is the AMD Virtex
-- UltraScale+ VU47P (AWS f2), the same one hdl/ecc and hdl/sha1 target.
--
-- Representation.  An element is 131 bits; bit (i-1) is the coefficient of
-- gamma_i = zeta^i + zeta^-i, i = 1..131, zeta a primitive 263rd root of
-- unity.  Two facts drive everything here, exactly as they drive the GPU
-- client:
--
--   * squaring is the coordinate permutation i -> fold(2 i), so sigma^k is
--     wiring and costs nothing;
--   * the Hamming weight is invariant under squaring, so the iteration
--     function is well defined on Frobenius orbits.
--
-- Multiplication is not done in the normal basis.  The multiplier converts
-- both operands to the optimal polynomial basis {c, c^2, ..., c^131},
-- c = gamma_1 (Bernstein-Lange), multiplies them as plain GF(2)[c]
-- polynomials, and converts the 261-coefficient product back.  Both
-- conversions are constant GF(2)-linear maps, built here at elaboration:
--
--   gamma_i = T_i(c)          T_0 = 0, T_1 = c, T_i = c T_{i-1} + T_{i-2}
--   c^k     = sum_j C(k,j) zeta^(k-2j),   C(k,j) odd iff (j and (k-j)) = 0
--   zeta^e + zeta^-e = gamma_fold(e)
--
-- No reduction step exists: the back-conversion absorbs it.  gf_mul_ref is
-- the direct normal-basis product, gamma_i gamma_j = gamma_fold(i+j) +
-- gamma_|i-j|, kept as an independent in-simulation oracle for the
-- testbenches.  hdl/ecc2k130/ecc2k_ref.py proves both against the client's
-- own field model.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package gf131_pkg is

  constant M : natural := 131;                 -- field degree
  constant N : natural := 2 * M + 1;           -- 263, prime, ord_263(2) = 131

  subtype gf_t    is std_logic_vector(M - 1 downto 0);
  -- a'(c) with A(c) = c a'(c): bit t is the coefficient of c^(t+1)
  subtype poly_t  is std_logic_vector(M - 1 downto 0);
  -- a'(c) b'(c): bit t is the coefficient of c^(t+2) of A(c) B(c)
  subtype dpoly_t is std_logic_vector(2 * M - 2 downto 0);
  subtype hw_t    is unsigned(7 downto 0);     -- weights 0..131

  constant GF_ZERO : gf_t := (others => '0');
  constant GF_ONE  : gf_t := (others => '1');  -- 1 = sum of all gamma_i

  -- ECC2K-130 distinguished-point cutoff (ecc2k130/generated/eccF131.h)
  constant DP_WEIGHT_DEFAULT : natural := 34;

  -- Multiplier pipeline: one input latch, one basis-conversion stage, a
  -- Karatsuba tree of MUL_KARATSUBA levels (each one pre-add stage and one
  -- post-combine stage around its three halves) with a one-stage schoolbook
  -- leaf, and one back-conversion stage.  Three levels give a 17-bit leaf;
  -- synthesised on the VU47P that is 4855 LUTs against 5547 at two levels
  -- and 5019 at four, with the best slack of the three (README, "The
  -- multiplier").
  constant MUL_KARATSUBA : natural := 3;
  constant MUL_LATENCY   : natural := 2 * MUL_KARATSUBA + 4;          -- 10

  -- A point whose x has weight 1 sits on no walk but has d = x + sigma^3(x)
  -- = gamma_1 + gamma_8 /= 0, so it can pad a batch without zeroing the
  -- product tree.  Padded leaves never produce output.
  constant DUMMY_X : gf_t := (0 => '1', others => '0');

  -- Weight in two clocks: 22 groups of 6 bits (one LUT6 per output bit),
  -- then the sum; or in three, with the 22 groups first summed four at a
  -- time (six partial sums of up to 24), which halves the depth of the
  -- adder tree where the full weight is needed.
  constant HW_GROUPS : natural := (M + 5) / 6;
  type hw_parts_t is array (0 to HW_GROUPS - 1) of unsigned(2 downto 0);
  constant HW_NQUAD : natural := (HW_GROUPS + 3) / 4;
  type hw_quads_t is array (0 to HW_NQUAD - 1) of unsigned(4 downto 0);

  function fold (e : integer) return natural;

  function gf_add     (a, b : gf_t) return gf_t;
  -- a^(2^k): coordinate permutation, pure wiring for constant k
  function gf_frob    (a : gf_t; k : natural) return gf_t;
  -- sigma^(3+s), s in 0..7, as sigma^3 then conditional sigma^1,2,4
  function gf_sigma_j (a : gf_t; s : std_logic_vector(2 downto 0)) return gf_t;
  function gf_weight  (a : gf_t) return hw_t;
  function gf_weight_parts (a : gf_t) return hw_parts_t;
  function hw_sum (p : hw_parts_t) return hw_t;
  function hw_quads (p : hw_parts_t) return hw_quads_t;
  function hw_sum (q : hw_quads_t) return hw_t;

  -- the multiplier's two constant linear maps
  function gf_prep   (a : gf_t)    return poly_t;
  function gf_to_onb (h : dpoly_t) return gf_t;

  -- schoolbook product over GF(2)[c], any width: 2n-1 bits from two n-bit
  function gf2_polymul (a, b : std_logic_vector) return std_logic_vector;

  -- independent oracle: direct product in the normal basis, O(m^2) loop
  function gf_mul_ref (a, b : gf_t) return gf_t;

end package;

package body gf131_pkg is

  function fold (e : integer) return natural is
    variable r : integer;
  begin
    r := e mod N;
    if r <= M then
      return r;
    else
      return N - r;
    end if;
  end function;

  function gf_add (a, b : gf_t) return gf_t is
  begin
    return a xor b;
  end function;

  function gf_frob (a : gf_t; k : natural) return gf_t is
    variable e : natural := 1;
    variable r : gf_t := (others => '0');
  begin
    for i in 1 to k loop
      e := (2 * e) mod N;
    end loop;
    for i in 1 to M loop
      r(fold(i * e) - 1) := a(i - 1);
    end loop;
    return r;
  end function;

  function gf_sigma_j (a : gf_t; s : std_logic_vector(2 downto 0)) return gf_t is
    variable r : gf_t;
  begin
    r := gf_frob(a, 3);
    if s(0) = '1' then r := gf_frob(r, 1); end if;
    if s(1) = '1' then r := gf_frob(r, 2); end if;
    if s(2) = '1' then r := gf_frob(r, 4); end if;
    return r;
  end function;

  function gf_weight (a : gf_t) return hw_t is
    variable cnt : natural := 0;
  begin
    for i in a'range loop
      if a(i) = '1' then
        cnt := cnt + 1;
      end if;
    end loop;
    return to_unsigned(cnt, hw_t'length);
  end function;

  function gf_weight_parts (a : gf_t) return hw_parts_t is
    variable p   : hw_parts_t;
    variable cnt : natural;
  begin
    for g in 0 to HW_GROUPS - 1 loop
      cnt := 0;
      for k in 0 to 5 loop
        if 6 * g + k < M then
          if a(6 * g + k) = '1' then
            cnt := cnt + 1;
          end if;
        end if;
      end loop;
      p(g) := to_unsigned(cnt, 3);
    end loop;
    return p;
  end function;

  function hw_sum (p : hw_parts_t) return hw_t is
    variable s : natural := 0;
  begin
    for g in 0 to HW_GROUPS - 1 loop
      s := s + to_integer(p(g));
    end loop;
    return to_unsigned(s, hw_t'length);
  end function;

  function hw_quads (p : hw_parts_t) return hw_quads_t is
    variable q : hw_quads_t;
    variable s : natural;
  begin
    for k in 0 to HW_NQUAD - 1 loop
      s := 0;
      for g in 4 * k to 4 * k + 3 loop
        if g < HW_GROUPS then
          s := s + to_integer(p(g));
        end if;
      end loop;
      q(k) := to_unsigned(s, 5);
    end loop;
    return q;
  end function;

  function hw_sum (q : hw_quads_t) return hw_t is
    variable s : natural := 0;
  begin
    for k in 0 to HW_NQUAD - 1 loop
      s := s + to_integer(q(k));
    end loop;
    return to_unsigned(s, hw_t'length);
  end function;

  function gf2_polymul (a, b : std_logic_vector) return std_logic_vector is
    constant NA : natural := a'length;
    constant NB : natural := b'length;
    alias aa : std_logic_vector(NA - 1 downto 0) is a;
    alias bb : std_logic_vector(NB - 1 downto 0) is b;
    variable r : std_logic_vector(NA + NB - 2 downto 0) := (others => '0');
  begin
    for i in 0 to NA - 1 loop
      for j in 0 to NB - 1 loop
        r(i + j) := r(i + j) xor (aa(i) and bb(j));
      end loop;
    end loop;
    return r;
  end function;

  -- ------------------------------------------------------------------ --
  -- gamma -> c-powers: the Dickson polynomials T_i(c), i = 1..M, over GF(2).
  -- Bit k of DICKSON(i) is the coefficient of c^k; deg T_i = i, no constant
  -- term ever appears (T_0 = 0, T_1 = c and the recurrence preserves it).
  -- ------------------------------------------------------------------ --
  type dickson_t is array (1 to M) of std_logic_vector(M downto 0);

  function build_dickson return dickson_t is
    variable t        : dickson_t;
    variable tm2, tm1 : std_logic_vector(M downto 0);
  begin
    tm2 := (others => '0');
    tm1 := (others => '0');
    tm1(1) := '1';
    t(1) := tm1;
    for i in 2 to M loop
      t(i) := (tm1(M - 1 downto 0) & '0') xor tm2;
      tm2  := tm1;
      tm1  := t(i);
    end loop;
    return t;
  end function;

  constant DICKSON : dickson_t := build_dickson;

  function gf_prep (a : gf_t) return poly_t is
    variable p : std_logic_vector(M downto 0) := (others => '0');
  begin
    for i in 1 to M loop
      if a(i - 1) = '1' then
        p := p xor DICKSON(i);
      end if;
    end loop;
    return p(M downto 1);
  end function;

  -- ------------------------------------------------------------------ --
  -- c-powers -> gamma.  TOONB(t) holds the coordinates of c^(t+2):
  --   c^k = (zeta + zeta^-1)^k = sum_j C(k,j) zeta^(k-2j)
  -- The j and k-j terms pair into one gamma_fold(k-2j); the middle term
  -- C(k, k/2) is even for every k >= 2, so no constant appears.
  -- ------------------------------------------------------------------ --
  type toonb_t is array (0 to 2 * M - 2) of gf_t;

  function build_toonb return toonb_t is
    variable t : toonb_t;
    variable r : gf_t;
    variable k : natural;
    variable lo, hi : unsigned(8 downto 0);
  begin
    for idx in 0 to 2 * M - 2 loop
      k := idx + 2;
      r := (others => '0');
      for j in 0 to (k - 1) / 2 loop
        lo := to_unsigned(j, 9);
        hi := to_unsigned(k - j, 9);
        if (lo and hi) = 0 then
          r(fold(k - 2 * j) - 1) := not r(fold(k - 2 * j) - 1);
        end if;
      end loop;
      t(idx) := r;
    end loop;
    return t;
  end function;

  constant TOONB : toonb_t := build_toonb;

  function gf_to_onb (h : dpoly_t) return gf_t is
    variable r : gf_t := (others => '0');
  begin
    for t in 0 to 2 * M - 2 loop
      if h(t) = '1' then
        r := r xor TOONB(t);
      end if;
    end loop;
    return r;
  end function;

  function gf_mul_ref (a, b : gf_t) return gf_t is
    variable r : gf_t := (others => '0');
    variable d : natural;
  begin
    for i in 1 to M loop
      if a(i - 1) = '1' then
        for j in 1 to M loop
          if b(j - 1) = '1' then
            r(fold(i + j) - 1) := not r(fold(i + j) - 1);
            if i /= j then
              d := abs(i - j);
              r(d - 1) := not r(d - 1);
            end if;
          end if;
        end loop;
      end if;
    end loop;
    return r;
  end function;

end package body;
