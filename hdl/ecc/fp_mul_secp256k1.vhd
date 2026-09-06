-- fp_mul_secp256k1.vhd
-- Fully pipelined 256-bit modular multiplier for secp256k1.
-- Initiation interval 1, latency MUL_LATENCY clocks.
--
-- Structure
-- ---------
--   in         one register stage latching the operands.
--
--   rows 0..7  operand scanning.  Stage i forms  t = acc + a(31:0)*b, emits
--              t(31:0) as the next product limb and carries t >> 32 into
--              stage i+1, with the multiplier operand shifted down 32 bits
--              each stage.  Each stage holds one 32 x 256 multiply -- eight
--              32x32 products, about 32 DSP48E2 slices -- and one 296-bit
--              add.  Eight stages, one full 256x256 product retired per
--              clock.
--
--              The emitted limbs are collected in a shift register: limb i
--              enters at the top and walks down, so after eight stages the
--              low half of the product is in natural order.  That keeps
--              every slice index locally static, which matters for both
--              simulation and synthesis.
--
--   red0..2    reduction by 2^256 = 2^32 + 977 (mod p), the same three folds
--              as Fp::secp_reduce in gpu/ecc/fp256.cuh:
--                U = Lo + Hi*977 + (Hi << 32)          < 2^289
--                V = U_lo + U_hi*977 + (U_hi << 32)    < 2^256 + 2^66
--                R = V_lo + V_hi*(2^32 + 977)          < 2^256
--              then one conditional subtraction of p.  The last fold cannot
--              carry out again because V_hi = 1 forces V_lo < 2^66.
--
-- The `tag` port travels with the operands so the caller can tell which of
-- its interleaved walks a result belongs to; the datapath never looks at
-- it.  That is what lets one multiplier serve many independent walks at
-- II = 1 with no stall logic at all.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.fp_pkg.all;

entity fp_mul_secp256k1 is
  generic (
    TAG_W : natural := 8
  );
  port (
    clk       : in  std_logic;
    rst       : in  std_logic;
    in_valid  : in  std_logic;
    in_a      : in  fp_t;
    in_b      : in  fp_t;
    in_tag    : in  std_logic_vector(TAG_W - 1 downto 0);
    out_valid : out std_logic;
    out_r     : out fp_t;
    out_tag   : out std_logic_vector(TAG_W - 1 downto 0)
  );
end entity;

architecture rtl of fp_mul_secp256k1 is

  -- 296 bits: the row accumulator never exceeds 2^289.
  constant ACC_W : natural := W + LIMB + 8;

  type acc_array_t is array (0 to MUL_ROW_STAGES) of unsigned(ACC_W - 1 downto 0);
  type fp_array_t  is array (0 to MUL_ROW_STAGES) of fp_t;

  signal acc  : acc_array_t := (others => (others => '0'));
  signal plo  : fp_array_t  := (others => (others => '0'));
  signal areg : fp_array_t  := (others => (others => '0'));
  signal breg : fp_array_t  := (others => (others => '0'));

  type tag_array_t is array (0 to MUL_LATENCY - 1) of std_logic_vector(TAG_W - 1 downto 0);
  signal vsr : std_logic_vector(MUL_LATENCY - 1 downto 0) := (others => '0');
  signal tsr : tag_array_t := (others => (others => '0'));

  signal r1_u   : unsigned(ACC_W - 1 downto 0) := (others => '0');
  signal r2_v   : unsigned(W + 8 downto 0)     := (others => '0');
  signal r3_res : fp_t                          := (others => '0');

begin

  -- ---------------------------------------------------------------- --
  -- input latch + operand-scanning rows
  -- ---------------------------------------------------------------- --
  rows : process (clk)
    variable t : unsigned(ACC_W - 1 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        acc  <= (others => (others => '0'));
        plo  <= (others => (others => '0'));
        areg <= (others => (others => '0'));
        breg <= (others => (others => '0'));
      else
        acc(0)  <= (others => '0');
        plo(0)  <= (others => '0');
        areg(0) <= in_a;
        breg(0) <= in_b;

        for i in 0 to MUL_ROW_STAGES - 1 loop
          t := acc(i) + resize(areg(i)(LIMB - 1 downto 0) * breg(i), ACC_W);
          plo(i + 1)  <= t(LIMB - 1 downto 0) & plo(i)(W - 1 downto LIMB);
          acc(i + 1)  <= shift_right(t, LIMB);
          areg(i + 1) <= shift_right(areg(i), LIMB);
          breg(i + 1) <= breg(i);
        end loop;
      end if;
    end if;
  end process;

  -- ---------------------------------------------------------------- --
  -- three reduction folds + conditional subtraction
  -- ---------------------------------------------------------------- --
  reduce : process (clk)
    variable uhi : unsigned(ACC_W - W - 1 downto 0);
    variable v   : unsigned(W + 8 downto 0);
    variable vlo : fp_t;
    variable wv  : unsigned(W downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        r1_u   <= (others => '0');
        r2_v   <= (others => '0');
        r3_res <= (others => '0');
      else
        -- red0: U = Lo + Hi*977 + (Hi << 32)
        r1_u <= resize(plo(MUL_ROW_STAGES), ACC_W)
              + resize(acc(MUL_ROW_STAGES)(W - 1 downto 0) * to_unsigned(RED_C_LO, 11), ACC_W)
              + shift_left(resize(acc(MUL_ROW_STAGES)(W - 1 downto 0), ACC_W), RED_C_SH);

        -- red1: V = U_lo + U_hi*977 + (U_hi << 32)
        uhi := r1_u(ACC_W - 1 downto W);
        v   := resize(r1_u(W - 1 downto 0), W + 9)
             + resize(uhi * to_unsigned(RED_C_LO, 11), W + 9)
             + shift_left(resize(uhi, W + 9), RED_C_SH);
        r2_v <= v;

        -- red2: R = V_lo + V_hi*(2^32 + 977), then conditional subtract
        vlo := r2_v(W - 1 downto 0);
        if r2_v(W) = '1' then
          wv := ('0' & vlo)
              + to_unsigned(RED_C_LO, W + 1)
              + shift_left(to_unsigned(1, W + 1), RED_C_SH);
        else
          wv := '0' & vlo;
        end if;
        if wv >= ('0' & P_MOD) then
          wv := wv - ('0' & P_MOD);
        end if;
        r3_res <= wv(W - 1 downto 0);
      end if;
    end if;
  end process;

  -- ---------------------------------------------------------------- --
  -- valid / tag travel alongside the datapath
  -- ---------------------------------------------------------------- --
  ctrl : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        vsr <= (others => '0');
        tsr <= (others => (others => '0'));
      else
        vsr(0) <= in_valid;
        tsr(0) <= in_tag;
        for i in 1 to MUL_LATENCY - 1 loop
          vsr(i) <= vsr(i - 1);
          tsr(i) <= tsr(i - 1);
        end loop;
      end if;
    end if;
  end process;

  out_valid <= vsr(MUL_LATENCY - 1);
  out_tag   <= tsr(MUL_LATENCY - 1);
  out_r     <= r3_res;

end architecture;
