-- sha1_pipe_prefix.vhd
-- Fixed-prefix variant of the pipelined SHA-1 compression core.
--
-- For chosen-prefix collision search, the attacker fixes the high-order
-- bytes of the candidate near-collision block (the prefix / birthday
-- conditioning) and varies only a few low-order message words across
-- billions of trials. Anything in the SHA-1 message schedule W[16..79]
-- that depends only on FIXED words is constant across all trials, and
-- synthesis can constant-fold it away -- shaving LUTs and shortening
-- the critical path of the affected stages.
--
-- This wrapper exposes:
--   * a generic FIXED_W[0..15]    : constant message words
--   * a generic VAR_MASK(15..0)   : '1' bits select positions overridden by
--                                   the runtime port var_w[].
--
-- Internally it assembles the full 512-bit block and instantiates the
-- vanilla sha1_pipe. The mux at each W[i] resolves to a wire when the
-- mask bit is '0' (constant) and to a 32-bit register/buffer when '1'.
-- Vivado/Synplify will then propagate the constants forward through the
-- message-schedule additions in W[16..79] automatically.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.sha1_pkg.all;

entity sha1_pipe_prefix is
  generic (
    TAG_WIDTH : positive       := 64;
    FIXED_W   : block_t        := (others => (others => '0'));
    VAR_MASK  : std_logic_vector(15 downto 0) := x"FFFF"
  );
  port (
    clk       : in  std_logic;
    rst       : in  std_logic;

    in_valid  : in  std_logic;
    var_w     : in  block_t;                                    -- positions where VAR_MASK='1'
    in_h      : in  state_t;
    in_tag    : in  std_logic_vector(TAG_WIDTH - 1 downto 0);

    out_valid : out std_logic;
    out_h     : out state_t;
    out_tag   : out std_logic_vector(TAG_WIDTH - 1 downto 0)
  );
end entity;

architecture wrap of sha1_pipe_prefix is
  signal m_full : block_t;
begin

  gen_w : for i in 0 to 15 generate
    fixed_g : if VAR_MASK(i) = '0' generate
      m_full(i) <= FIXED_W(i);
    end generate;
    var_g : if VAR_MASK(i) = '1' generate
      m_full(i) <= var_w(i);
    end generate;
  end generate;

  core : entity work.sha1_pipe
    generic map (TAG_WIDTH => TAG_WIDTH)
    port map (
      clk       => clk,
      rst       => rst,
      in_valid  => in_valid,
      in_m      => m_full,
      in_h      => in_h,
      in_tag    => in_tag,
      out_valid => out_valid,
      out_h     => out_h,
      out_tag   => out_tag
    );

end architecture;
