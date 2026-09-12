-- gf131_mul.vhd
-- Pipelined GF(2^131) multiplier, one product per clock, fixed latency
-- MUL_LATENCY, operands and result in the permuted type-II normal basis.
--
--   stage 0          latch a, b
--   stage 1          a' = prep(a), b' = prep(b)      (gamma -> c-powers)
--   stages 2..9      h ^= a' * b'[digit r] << r*D    (8 rows of 17 bits)
--   stage 10         r = to_onb(h)                   (c-powers -> gamma)
--
-- The rows are operand scanning over GF(2)[c]: row r ANDs a' with each of
-- the MUL_DIGIT bits of its digit of b', XORs the shifted copies together
-- and folds them into the 261-bit running product.  No carries, so a row is
-- 131 x 17 AND gates and an XOR tree five deep -- two or three LUT6 levels,
-- which is what keeps the clock high.  There is no reduction: the
-- back-conversion maps every c^k, k = 2..262, straight to normal-basis
-- coordinates.
--
-- Cost per multiplier (derived from the constant matrices, see README):
-- about 17k AND + 17k XOR for the product and 4.4k XOR for the two
-- conversions; zero DSPs.  The tag rides alongside and is never inspected,
-- so one multiplier serves any number of independent contexts at II = 1.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.gf131_pkg.all;

entity gf131_mul is
  generic (
    TAG_W : natural := 8
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
    out_tag   : out std_logic_vector(TAG_W - 1 downto 0)
  );
end entity;

architecture rtl of gf131_mul is

  constant D    : natural := MUL_DIGIT;
  constant ROWS : natural := MUL_ROWS;
  constant PBW  : natural := ROWS * D;          -- b' zero-extended to whole digits

  subtype tag_t is std_logic_vector(TAG_W - 1 downto 0);
  subtype pbx_t is std_logic_vector(PBW - 1 downto 0);
  subtype pp_t  is std_logic_vector(M + D - 2 downto 0);

  type acc_arr is array (0 to ROWS) of dpoly_t;
  type pa_arr  is array (0 to ROWS) of poly_t;
  type pb_arr  is array (0 to ROWS) of pbx_t;
  type tag_arr is array (0 to ROWS) of tag_t;
  type v_arr   is array (0 to ROWS) of std_logic;

  -- stage 0
  signal s0_v   : std_logic := '0';
  signal s0_a   : gf_t := (others => '0');
  signal s0_b   : gf_t := (others => '0');
  signal s0_tag : tag_t := (others => '0');

  -- rows: index r is the input of row r; index ROWS is the finished product
  signal acc   : acc_arr := (others => (others => '0'));
  signal pa    : pa_arr  := (others => (others => '0'));
  signal pb    : pb_arr  := (others => (others => '0'));
  signal tg    : tag_arr := (others => (others => '0'));
  signal vl    : v_arr   := (others => '0');

  -- a' times one digit of b', unshifted
  function row_pp (a : poly_t; dig : std_logic_vector(D - 1 downto 0)) return pp_t is
    variable r  : pp_t := (others => '0');
    variable sh : pp_t;
  begin
    for q in 0 to D - 1 loop
      if dig(q) = '1' then
        sh := (others => '0');
        sh(q + M - 1 downto q) := a;
        r := r xor sh;
      end if;
    end loop;
    return r;
  end function;

begin

  -- ---------------------------------------------------------------- --
  -- stage 0: input latch; stage 1: basis conversion
  -- ---------------------------------------------------------------- --
  front : process (clk)
    variable pbx : pbx_t;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        s0_v  <= '0';
        vl(0) <= '0';
      else
        s0_v   <= in_valid;
        s0_a   <= in_a;
        s0_b   <= in_b;
        s0_tag <= in_tag;

        vl(0)  <= s0_v;
        tg(0)  <= s0_tag;
        pa(0)  <= gf_prep(s0_a);
        pbx    := (others => '0');
        pbx(M - 1 downto 0) := gf_prep(s0_b);
        pb(0)  <= pbx;
        acc(0) <= (others => '0');
      end if;
    end if;
  end process;

  -- ---------------------------------------------------------------- --
  -- rows
  -- ---------------------------------------------------------------- --
  rows_g : for r in 0 to ROWS - 1 generate
    row : process (clk)
      variable p   : pp_t;
      variable nxt : dpoly_t;
    begin
      if rising_edge(clk) then
        if rst = '1' then
          vl(r + 1) <= '0';
        else
          p   := row_pp(pa(r), pb(r)(r * D + D - 1 downto r * D));
          nxt := acc(r);
          for t in 0 to M + D - 2 loop
            if r * D + t <= 2 * M - 2 then
              nxt(r * D + t) := nxt(r * D + t) xor p(t);
            end if;
          end loop;
          acc(r + 1) <= nxt;
          pa(r + 1)  <= pa(r);
          pb(r + 1)  <= pb(r);
          tg(r + 1)  <= tg(r);
          vl(r + 1)  <= vl(r);
        end if;
      end if;
    end process;
  end generate;

  -- ---------------------------------------------------------------- --
  -- last stage: back to the normal basis
  -- ---------------------------------------------------------------- --
  back : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        out_valid <= '0';
      else
        out_valid <= vl(ROWS);
        out_tag   <= tg(ROWS);
        out_r     <= gf_to_onb(acc(ROWS));
      end if;
    end if;
  end process;

end architecture;
