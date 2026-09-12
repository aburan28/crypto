-- gf131_tb_pkg.vhd
-- Simulation-only helpers shared by the hdl/ecc2k130 testbenches: the hex
-- format of vectors_ecc2k130.txt (33 digits, 132 bits, top bit zero).

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

package gf131_tb_pkg is

  constant HEXW : natural := (M + 3) / 4;      -- 33 hex digits

  function hex_to_gf (s : string) return gf_t;
  function gf_to_hex (v : gf_t) return string;

end package;

package body gf131_tb_pkg is

  function hex_to_gf (s : string) return gf_t is
    variable r : unsigned(4 * HEXW - 1 downto 0) := (others => '0');
    variable d : natural;
  begin
    for i in s'range loop
      case s(i) is
        when '0' to '9' => d := character'pos(s(i)) - character'pos('0');
        when 'a' to 'f' => d := character'pos(s(i)) - character'pos('a') + 10;
        when 'A' to 'F' => d := character'pos(s(i)) - character'pos('A') + 10;
        when others     => d := 0;
      end case;
      r := shift_left(r, 4) or to_unsigned(d, r'length);
    end loop;
    assert r(4 * HEXW - 1 downto M) = 0
      report "hex field wider than the field: " & s severity failure;
    return std_logic_vector(r(M - 1 downto 0));
  end function;

  function gf_to_hex (v : gf_t) return string is
    constant DIG : string(1 to 16) := "0123456789abcdef";
    variable u : unsigned(4 * HEXW - 1 downto 0) := (others => '0');
    variable s : string(1 to HEXW);
    variable n : natural;
  begin
    u(M - 1 downto 0) := unsigned(v);
    for i in 0 to HEXW - 1 loop
      n := to_integer(u((HEXW - 1 - i) * 4 + 3 downto (HEXW - 1 - i) * 4));
      s(i + 1) := DIG(n + 1);
    end loop;
    return s;
  end function;

end package body;
