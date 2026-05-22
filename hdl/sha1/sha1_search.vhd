-- sha1_search.vhd
-- Counter-driven SHA-1 search engine with distinguished-point output.
--
-- Drives an incrementing counter through the variable portion of a
-- pipelined SHA-1 core and emits a hit on the AXI-Stream-style output
-- whenever the resulting digest has at least DP_ZERO_BITS leading
-- (high-order) zero bits. This is the building block for the birthday /
-- distinguished-point phase of a near-collision search: a host process
-- collects hits, follows trails, and looks for collisions across trails.
--
-- The output bus is 256 bits per beat: { counter[63:0] | digest[159:0] | reserved[31:0] }.
-- The host should drain it; backpressure stalls the message feeder via
-- the small built-in FIFO (depth = HIT_FIFO_DEPTH).
--
-- Throughput accounting:
--   * One trial per clock (after fill). At 500 MHz that's 5e8 trials/s
--     per engine, ~75 GH/s if ~150 engines fit on a VU47P.
--   * Hit rate ~ 2^-DP_ZERO_BITS per trial. For DP_ZERO_BITS=32,
--     ~0.12 hits/sec/engine -- well within a modest AXI-Stream link.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.sha1_pkg.all;

entity sha1_search is
  generic (
    DP_ZERO_BITS   : natural range 1 to 64 := 32;
    HIT_FIFO_DEPTH : positive              := 16;
    FIXED_W        : block_t               := (others => (others => '0'));
    VAR_MASK       : std_logic_vector(15 downto 0) := x"0001";   -- only W[0] varies
    COUNTER_LANE   : integer range 0 to 15 := 0                  -- which W[i] gets the counter
  );
  port (
    clk        : in  std_logic;
    rst        : in  std_logic;

    -- Engine control
    enable     : in  std_logic;
    seed       : in  unsigned(63 downto 0);     -- counter start value
    in_h       : in  state_t;                   -- chaining input (constant for the run)

    -- AXI-Stream-like hit output
    m_axis_tdata  : out std_logic_vector(255 downto 0);
    m_axis_tvalid : out std_logic;
    m_axis_tready : in  std_logic;

    -- Telemetry
    trials_q   : out unsigned(63 downto 0);
    hits_q     : out unsigned(31 downto 0)
  );
end entity;

architecture rtl of sha1_search is

  -- Message feeder
  signal counter   : unsigned(63 downto 0);
  signal var_w     : block_t;
  signal feed_valid: std_logic;
  signal feed_tag  : std_logic_vector(63 downto 0);

  -- Core output
  signal core_valid : std_logic;
  signal core_h     : state_t;
  signal core_tag   : std_logic_vector(63 downto 0);

  -- DP detection
  signal is_dp     : std_logic;

  -- Tiny ring-buffer hit FIFO
  type fifo_word_t is record
    digest : state_t;
    tag    : std_logic_vector(63 downto 0);
  end record;
  type   fifo_mem_t is array (0 to HIT_FIFO_DEPTH - 1) of fifo_word_t;
  signal fifo_mem  : fifo_mem_t;
  signal wr_ptr    : integer range 0 to HIT_FIFO_DEPTH - 1 := 0;
  signal rd_ptr    : integer range 0 to HIT_FIFO_DEPTH - 1 := 0;
  signal count_q   : integer range 0 to HIT_FIFO_DEPTH     := 0;
  signal fifo_full : std_logic;
  signal fifo_emp  : std_logic;

  -- Telemetry
  signal trials_r  : unsigned(63 downto 0) := (others => '0');
  signal hits_r    : unsigned(31 downto 0) := (others => '0');

begin

  -- ===== Message feeder ============================================================
  -- The counter is mapped into a single W-lane; other variable lanes (per VAR_MASK)
  -- are tied to the upper bits of the counter so they sweep on overflow.
  process(clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        counter    <= (others => '0');
        feed_valid <= '0';
      elsif enable = '1' then
        if feed_valid = '0' then
          counter <= seed;          -- first cycle of a run
        else
          counter <= counter + 1;
        end if;
        feed_valid <= '1';
      else
        feed_valid <= '0';
      end if;
    end if;
  end process;

  feed_tag <= std_logic_vector(counter);

  gen_lanes : for i in 0 to 15 generate
    is_counter : if i = COUNTER_LANE generate
      var_w(i) <= counter(31 downto 0) when VAR_MASK(i) = '1'
                  else FIXED_W(i);
    end generate;
    other_lanes : if i /= COUNTER_LANE generate
      var_w(i) <= counter(63 downto 32) when VAR_MASK(i) = '1'
                  else FIXED_W(i);
    end generate;
  end generate;

  -- ===== Pipelined core ============================================================
  engine : entity work.sha1_pipe_prefix
    generic map (
      TAG_WIDTH => 64,
      FIXED_W   => FIXED_W,
      VAR_MASK  => VAR_MASK
    )
    port map (
      clk       => clk,
      rst       => rst,
      in_valid  => feed_valid,
      var_w     => var_w,
      in_h      => in_h,
      in_tag    => feed_tag,
      out_valid => core_valid,
      out_h     => core_h,
      out_tag   => core_tag
    );

  -- ===== Distinguished-point detector =============================================
  -- DP if the top DP_ZERO_BITS bits of (H0 || H1 || ...) are zero. We check the
  -- high-order word(s) of the digest. For DP_ZERO_BITS <= 32 we only need H0.
  process(core_valid, core_h)
    variable z : std_logic;
    variable bits_left : integer;
    variable wi : integer;
    variable bi : integer;
  begin
    z := '1';
    if core_valid = '0' then
      z := '0';
    else
      bits_left := DP_ZERO_BITS;
      wi := 0;
      while bits_left > 0 and wi < 5 loop
        bi := 31;
        while bi >= 0 and bits_left > 0 loop
          if core_h(wi)(bi) /= '0' then
            z := '0';
          end if;
          bi := bi - 1;
          bits_left := bits_left - 1;
        end loop;
        wi := wi + 1;
      end loop;
    end if;
    is_dp <= z;
  end process;

  -- ===== Hit FIFO ==================================================================
  fifo_full <= '1' when count_q = HIT_FIFO_DEPTH else '0';
  fifo_emp  <= '1' when count_q = 0              else '0';

  process(clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        wr_ptr  <= 0;
        rd_ptr  <= 0;
        count_q <= 0;
      else
        -- Push: DP hit and room available
        if is_dp = '1' and fifo_full = '0' then
          fifo_mem(wr_ptr).digest <= core_h;
          fifo_mem(wr_ptr).tag    <= core_tag;
          if wr_ptr = HIT_FIFO_DEPTH - 1 then wr_ptr <= 0;
          else                                wr_ptr <= wr_ptr + 1;
          end if;
        end if;
        -- Pop: downstream took a beat
        if m_axis_tready = '1' and fifo_emp = '0' then
          if rd_ptr = HIT_FIFO_DEPTH - 1 then rd_ptr <= 0;
          else                                rd_ptr <= rd_ptr + 1;
          end if;
        end if;
        -- Occupancy
        if    is_dp = '1' and fifo_full = '0' and (m_axis_tready = '0' or fifo_emp = '1') then
          count_q <= count_q + 1;
        elsif (is_dp = '0' or fifo_full = '1') and m_axis_tready = '1' and fifo_emp = '0' then
          count_q <= count_q - 1;
        end if;
      end if;
    end if;
  end process;

  m_axis_tvalid <= not fifo_emp;
  m_axis_tdata  <=
      fifo_mem(rd_ptr).tag                                             -- [255:192] counter
    & std_logic_vector(fifo_mem(rd_ptr).digest(0))                     -- [191:160] H0
    & std_logic_vector(fifo_mem(rd_ptr).digest(1))                     -- [159:128] H1
    & std_logic_vector(fifo_mem(rd_ptr).digest(2))                     -- [127: 96] H2
    & std_logic_vector(fifo_mem(rd_ptr).digest(3))                     -- [ 95: 64] H3
    & std_logic_vector(fifo_mem(rd_ptr).digest(4))                     -- [ 63: 32] H4
    & x"00000000";                                                     -- [ 31:  0] reserved

  -- ===== Telemetry =================================================================
  process(clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        trials_r <= (others => '0');
        hits_r   <= (others => '0');
      else
        if core_valid = '1' then trials_r <= trials_r + 1; end if;
        if is_dp     = '1' and fifo_full = '0' then
          hits_r <= hits_r + 1;
        end if;
      end if;
    end if;
  end process;

  trials_q <= trials_r;
  hits_q   <= hits_r;

end architecture;
