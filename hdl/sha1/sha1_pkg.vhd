-- sha1_pkg.vhd
-- Shared types, round constants, and helper functions for the pipelined
-- SHA-1 compression core. Targets AMD/Xilinx Virtex UltraScale+ (e.g. the
-- VU47P on AWS f2.6xlarge), but is plain IEEE.numeric_std VHDL-2008 and
-- portable to other synthesis flows.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package sha1_pkg is

  subtype word_t       is unsigned(31 downto 0);
  type    word_array_t is array (natural range <>) of word_t;
  subtype state_t      is word_array_t(0 to 4);   -- a,b,c,d,e or H0..H4
  subtype block_t      is word_array_t(0 to 15);  -- 16 big-endian 32-bit words

  -- FIPS 180-4 round constants
  constant K0 : word_t := x"5A827999";  -- t in  0..19
  constant K1 : word_t := x"6ED9EBA1";  -- t in 20..39
  constant K2 : word_t := x"8F1BBCDC";  -- t in 40..59
  constant K3 : word_t := x"CA62C1D6";  -- t in 60..79

  constant H_INIT : state_t := (
    x"67452301", x"EFCDAB89", x"98BADCFE", x"10325476", x"C3D2E1F0"
  );

  function rotl  (x : word_t; n : natural) return word_t;
  function f_ch (b, c, d : word_t) return word_t;  -- t in  0..19
  function f_par(b, c, d : word_t) return word_t;  -- t in 20..39 and 60..79
  function f_maj(b, c, d : word_t) return word_t;  -- t in 40..59

  function k_of (t : natural) return word_t;
  function f_of (t : natural; b, c, d : word_t) return word_t;

end package;

package body sha1_pkg is

  function rotl(x : word_t; n : natural) return word_t is
  begin
    if n mod 32 = 0 then
      return x;
    end if;
    return x(31 - (n mod 32) downto 0) & x(31 downto 32 - (n mod 32));
  end function;

  function f_ch (b, c, d : word_t) return word_t is
  begin
    return (b and c) or ((not b) and d);
  end function;

  function f_par(b, c, d : word_t) return word_t is
  begin
    return b xor c xor d;
  end function;

  function f_maj(b, c, d : word_t) return word_t is
  begin
    return (b and c) or (b and d) or (c and d);
  end function;

  function k_of(t : natural) return word_t is
  begin
    if    t < 20 then return K0;
    elsif t < 40 then return K1;
    elsif t < 60 then return K2;
    else              return K3;
    end if;
  end function;

  function f_of(t : natural; b, c, d : word_t) return word_t is
  begin
    if    t < 20 then return f_ch (b, c, d);
    elsif t < 40 then return f_par(b, c, d);
    elsif t < 60 then return f_maj(b, c, d);
    else              return f_par(b, c, d);
    end if;
  end function;

end package body;
