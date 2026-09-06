-- fp_pkg.vhd
-- Shared constants and types for the secp256k1 field / Pollard-rho engine.
-- Plain IEEE.numeric_std VHDL-2008; targets AMD/Xilinx UltraScale+ (the
-- VU47P on AWS f2 is the reference part) but contains no vendor primitives,
-- so it simulates under GHDL and synthesises anywhere.
--
-- The arithmetic here mirrors gpu/ecc/fp256.cuh exactly: same prime, same
-- three-fold special reduction, same conditional-subtraction bound.  The
-- testbenches drive both with vectors emitted by gpu/ecc/ecref.py, so the
-- FPGA and GPU paths are checked against one oracle.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package fp_pkg is

  constant W        : natural := 256;                -- field width
  constant LIMB     : natural := 32;                 -- limb width
  constant NLIMB    : natural := W / LIMB;           -- 8 limbs

  subtype  fp_t     is unsigned(W - 1 downto 0);
  subtype  dbl_t    is unsigned(2 * W - 1 downto 0);

  -- secp256k1: p = 2^256 - 2^32 - 977
  constant P_MOD : fp_t :=
    x"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";

  -- The reduction constant C with 2^256 = C (mod p), C = 2^32 + 977.
  constant RED_C_LO : natural := 977;
  constant RED_C_SH : natural := 32;

  -- Multiplier pipeline depth: one input latch, NLIMB row-accumulation
  -- stages, three reduction stages.
  constant MUL_IN_STAGES  : natural := 1;
  constant MUL_ROW_STAGES : natural := NLIMB;   -- 8
  constant MUL_RED_STAGES : natural := 3;
  constant MUL_LATENCY    : natural :=
    MUL_IN_STAGES + MUL_ROW_STAGES + MUL_RED_STAGES;   -- 12

  -- Modular add / subtract, combinational helpers (one clock when
  -- registered by the caller).  Inputs must already be < p.
  function fp_add (a, b : fp_t) return fp_t;
  function fp_sub (a, b : fp_t) return fp_t;
  function fp_dbl (a    : fp_t) return fp_t;
  function fp_neg (a    : fp_t) return fp_t;

  -- True when a > (p-1)/2, the negation-map test.
  function fp_gt_half (a : fp_t) return boolean;

end package;

package body fp_pkg is

  function fp_add (a, b : fp_t) return fp_t is
    variable s : unsigned(W downto 0);
  begin
    s := ('0' & a) + ('0' & b);
    if s >= ('0' & P_MOD) then
      s := s - ('0' & P_MOD);
    end if;
    return s(W - 1 downto 0);
  end function;

  function fp_sub (a, b : fp_t) return fp_t is
    variable d : unsigned(W downto 0);
  begin
    d := ('0' & a) - ('0' & b);
    if d(W) = '1' then                    -- borrow
      d := d + ('0' & P_MOD);
    end if;
    return d(W - 1 downto 0);
  end function;

  function fp_dbl (a : fp_t) return fp_t is
  begin
    return fp_add(a, a);
  end function;

  function fp_neg (a : fp_t) return fp_t is
  begin
    if a = 0 then
      return a;
    end if;
    return P_MOD - a;
  end function;

  function fp_gt_half (a : fp_t) return boolean is
    -- (p-1)/2; p is odd so this is exact
    constant HALF : fp_t := shift_right(P_MOD - 1, 1);
  begin
    return a > HALF;
  end function;

end package body;
