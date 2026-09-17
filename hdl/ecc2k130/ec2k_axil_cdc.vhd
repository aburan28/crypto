-- ec2k_axil_cdc.vhd
-- AXI4-Lite clock-domain bridge: a slave on s_clk, a master on m_clk, one
-- transaction of each kind in flight, no relationship assumed between the
-- clocks.  On AWS F2 it sits between the shell's OCL port (clk_main_a0,
-- 250 MHz) and ec2k_axil on the engine clock, so the engines can run at
-- whatever the fabric closes rather than what the shell hands out.
--
-- Each channel is a four-phase level handshake: the slow side latches the
-- request (address, data, strobe) into cdc_*_src registers, then raises a
-- request flag one clock later; the fast side sees the flag through a
-- synchroniser, copies the request into cdc_*_cap registers, performs the
-- transaction on the master port, parks the response (read data) in a
-- cdc_r_src register and raises an ack flag; the slow side sees the ack,
-- copies the response, completes the AXI transaction and drops the
-- request; the fast side drops the ack when it sees the request fall.
-- Data therefore changes only while its flag is low and is sampled only
-- after the flag has been seen high through SYNC_FF stages, so the only
-- timing requirement on a data path is that it settle within a source
-- clock period plus the synchroniser latency.  Every register that a
-- cross-domain path ends at is named cdc_*_meta (flags) or cdc_*_cap
-- (data), which is what the constraints in
-- aws/cl_ecc2k130/build/constraints/cl_timing_user.xdc key on: false paths
-- into the _meta flops (ASYNC_REG), set_max_delay -datapath_only from
-- cdc_*_src to cdc_*_cap.
--
-- A round trip costs roughly 2 * (SYNC_FF + 1) clocks of each domain plus
-- the transaction itself, ~50 ns at 250 / 333 MHz; the host's PCIe access
-- is twenty times that, and it makes a dozen of them per distinguished
-- point, one every 2**26 steps per walk.
--
-- Write responses are always OKAY (ec2k_axil never errors); the read
-- response is carried across.  Write address and data may arrive on
-- different clocks, as AXI allows; reads and writes are independent.

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity ec2k_axil_cdc is
  generic (
    SYNC_FF : natural := 3               -- synchroniser stages per flag
  );
  port (
    -- slave side
    s_clk     : in  std_logic;
    s_rst     : in  std_logic;           -- synchronous, active high
    s_awaddr  : in  std_logic_vector(31 downto 0);
    s_awvalid : in  std_logic;
    s_awready : out std_logic;
    s_wdata   : in  std_logic_vector(31 downto 0);
    s_wstrb   : in  std_logic_vector(3 downto 0);
    s_wvalid  : in  std_logic;
    s_wready  : out std_logic;
    s_bresp   : out std_logic_vector(1 downto 0);
    s_bvalid  : out std_logic;
    s_bready  : in  std_logic;
    s_araddr  : in  std_logic_vector(31 downto 0);
    s_arvalid : in  std_logic;
    s_arready : out std_logic;
    s_rdata   : out std_logic_vector(31 downto 0);
    s_rresp   : out std_logic_vector(1 downto 0);
    s_rvalid  : out std_logic;
    s_rready  : in  std_logic;
    -- master side
    m_clk     : in  std_logic;
    m_rst     : in  std_logic;           -- synchronous, active high
    m_awaddr  : out std_logic_vector(31 downto 0);
    m_awvalid : out std_logic;
    m_awready : in  std_logic;
    m_wdata   : out std_logic_vector(31 downto 0);
    m_wstrb   : out std_logic_vector(3 downto 0);
    m_wvalid  : out std_logic;
    m_wready  : in  std_logic;
    m_bresp   : in  std_logic_vector(1 downto 0);
    m_bvalid  : in  std_logic;
    m_bready  : out std_logic;
    m_araddr  : out std_logic_vector(31 downto 0);
    m_arvalid : out std_logic;
    m_arready : in  std_logic;
    m_rdata   : in  std_logic_vector(31 downto 0);
    m_rresp   : in  std_logic_vector(1 downto 0);
    m_rvalid  : in  std_logic;
    m_rready  : out std_logic
  );
end entity;

architecture rtl of ec2k_axil_cdc is

  subtype word_t is std_logic_vector(31 downto 0);
  subtype sync_t is std_logic_vector(SYNC_FF - 2 downto 0);

  -- slow side: request sources and flags
  signal aw_got, w_got : std_logic := '0';
  signal cdc_w_src_addr : word_t := (others => '0');
  signal cdc_w_src_data : word_t := (others => '0');
  signal cdc_w_src_strb : std_logic_vector(3 downto 0) := (others => '0');
  signal cdc_wreq       : std_logic := '0';
  signal w_busy         : std_logic := '0';
  signal bvalid         : std_logic := '0';

  signal ar_got         : std_logic := '0';
  signal cdc_r_src_addr : word_t := (others => '0');
  signal cdc_rreq       : std_logic := '0';
  signal r_busy         : std_logic := '0';
  signal rvalid         : std_logic := '0';
  signal cdc_r_cap_data : word_t := (others => '0');
  signal cdc_r_cap_resp : std_logic_vector(1 downto 0) := "00";

  -- slow side: acks arriving from the fast side
  signal cdc_wack_meta, cdc_rack_meta : std_logic := '0';
  signal cdc_wack_sync, cdc_rack_sync : sync_t := (others => '0');
  signal wack_s, rack_s : std_logic;

  -- fast side: requests arriving from the slow side
  signal cdc_wreq_meta, cdc_rreq_meta : std_logic := '0';
  signal cdc_wreq_sync, cdc_rreq_sync : sync_t := (others => '0');
  signal wreq_m, rreq_m : std_logic;

  -- fast side: captured requests, transaction state, response source
  type wst_t is (W_IDLE, W_ISSUE, W_RESP, W_ACK);
  type rst_t is (R_IDLE, R_ISSUE, R_RESP, R_FLAG, R_ACK);
  signal wst : wst_t := W_IDLE;
  signal rs  : rst_t := R_IDLE;
  signal cdc_w_cap_addr : word_t := (others => '0');
  signal cdc_w_cap_data : word_t := (others => '0');
  signal cdc_w_cap_strb : std_logic_vector(3 downto 0) := (others => '0');
  signal cdc_r_cap_addr : word_t := (others => '0');
  signal awv, wv, arv   : std_logic := '0';
  signal cdc_r_src_data : word_t := (others => '0');
  signal cdc_r_src_resp : std_logic_vector(1 downto 0) := "00";
  signal cdc_wack, cdc_rack : std_logic := '0';

  attribute ASYNC_REG : string;
  attribute ASYNC_REG of cdc_wack_meta, cdc_rack_meta, cdc_wreq_meta, cdc_rreq_meta : signal is "TRUE";
  attribute ASYNC_REG of cdc_wack_sync, cdc_rack_sync, cdc_wreq_sync, cdc_rreq_sync : signal is "TRUE";

begin

  assert SYNC_FF >= 2 report "SYNC_FF must be at least 2" severity failure;

  wack_s <= cdc_wack_sync(SYNC_FF - 2);
  rack_s <= cdc_rack_sync(SYNC_FF - 2);
  wreq_m <= cdc_wreq_sync(SYNC_FF - 2);
  rreq_m <= cdc_rreq_sync(SYNC_FF - 2);

  ---------------------------------------------------------------------------
  -- slow side
  ---------------------------------------------------------------------------
  s_awready <= not aw_got;
  s_wready  <= not w_got;
  s_bvalid  <= bvalid;
  s_bresp   <= "00";
  s_arready <= not ar_got;
  s_rvalid  <= rvalid;
  s_rdata   <= cdc_r_cap_data;
  s_rresp   <= cdc_r_cap_resp;

  slow : process (s_clk)
  begin
    if rising_edge(s_clk) then
      cdc_wack_meta <= cdc_wack;
      cdc_wack_sync <= cdc_wack_sync(SYNC_FF - 3 downto 0) & cdc_wack_meta;
      cdc_rack_meta <= cdc_rack;
      cdc_rack_sync <= cdc_rack_sync(SYNC_FF - 3 downto 0) & cdc_rack_meta;

      -- write: latch address and data, raise the request once both are in,
      -- complete on the ack, and wait for the ack to fall before the next
      if s_awvalid = '1' and aw_got = '0' then
        cdc_w_src_addr <= s_awaddr;
        aw_got <= '1';
      end if;
      if s_wvalid = '1' and w_got = '0' then
        cdc_w_src_data <= s_wdata;
        cdc_w_src_strb <= s_wstrb;
        w_got <= '1';
      end if;
      if bvalid = '1' and s_bready = '1' then
        bvalid <= '0';
      end if;
      if cdc_wreq = '0' and w_busy = '0' and aw_got = '1' and w_got = '1'
         and (bvalid = '0' or s_bready = '1') then
        cdc_wreq <= '1';
        w_busy   <= '1';
      elsif cdc_wreq = '1' and wack_s = '1' then
        cdc_wreq <= '0';
        aw_got   <= '0';
        w_got    <= '0';
        bvalid   <= '1';
      elsif cdc_wreq = '0' and w_busy = '1' and wack_s = '0' then
        w_busy <= '0';
      end if;

      -- read: the same shape, with the response copied on the ack
      if s_arvalid = '1' and ar_got = '0' then
        cdc_r_src_addr <= s_araddr;
        ar_got <= '1';
      end if;
      if rvalid = '1' and s_rready = '1' then
        rvalid <= '0';
      end if;
      if cdc_rreq = '0' and r_busy = '0' and ar_got = '1'
         and (rvalid = '0' or s_rready = '1') then
        cdc_rreq <= '1';
        r_busy   <= '1';
      elsif cdc_rreq = '1' and rack_s = '1' then
        cdc_r_cap_data <= cdc_r_src_data;
        cdc_r_cap_resp <= cdc_r_src_resp;
        cdc_rreq <= '0';
        ar_got   <= '0';
        rvalid   <= '1';
      elsif cdc_rreq = '0' and r_busy = '1' and rack_s = '0' then
        r_busy <= '0';
      end if;

      if s_rst = '1' then
        aw_got   <= '0';
        w_got    <= '0';
        ar_got   <= '0';
        cdc_wreq <= '0';
        cdc_rreq <= '0';
        w_busy   <= '0';
        r_busy   <= '0';
        bvalid   <= '0';
        rvalid   <= '0';
      end if;
    end if;
  end process;

  ---------------------------------------------------------------------------
  -- fast side
  ---------------------------------------------------------------------------
  m_awaddr  <= cdc_w_cap_addr;
  m_wdata   <= cdc_w_cap_data;
  m_wstrb   <= cdc_w_cap_strb;
  m_awvalid <= awv;
  m_wvalid  <= wv;
  m_bready  <= '1' when wst = W_RESP else '0';
  m_araddr  <= cdc_r_cap_addr;
  m_arvalid <= arv;
  m_rready  <= '1' when rs = R_RESP else '0';

  fast : process (m_clk)
  begin
    if rising_edge(m_clk) then
      cdc_wreq_meta <= cdc_wreq;
      cdc_wreq_sync <= cdc_wreq_sync(SYNC_FF - 3 downto 0) & cdc_wreq_meta;
      cdc_rreq_meta <= cdc_rreq;
      cdc_rreq_sync <= cdc_rreq_sync(SYNC_FF - 3 downto 0) & cdc_rreq_meta;

      case wst is
        when W_IDLE =>
          if wreq_m = '1' then
            cdc_w_cap_addr <= cdc_w_src_addr;
            cdc_w_cap_data <= cdc_w_src_data;
            cdc_w_cap_strb <= cdc_w_src_strb;
            awv <= '1';
            wv  <= '1';
            wst <= W_ISSUE;
          end if;
        when W_ISSUE =>
          if awv = '1' and m_awready = '1' then
            awv <= '0';
          end if;
          if wv = '1' and m_wready = '1' then
            wv <= '0';
          end if;
          if (awv = '0' or m_awready = '1') and (wv = '0' or m_wready = '1') then
            wst <= W_RESP;
          end if;
        when W_RESP =>
          if m_bvalid = '1' then
            cdc_wack <= '1';
            wst <= W_ACK;
          end if;
        when W_ACK =>
          if wreq_m = '0' then
            cdc_wack <= '0';
            wst <= W_IDLE;
          end if;
      end case;

      case rs is
        when R_IDLE =>
          if rreq_m = '1' then
            cdc_r_cap_addr <= cdc_r_src_addr;
            arv <= '1';
            rs  <= R_ISSUE;
          end if;
        when R_ISSUE =>
          if m_arready = '1' then
            arv <= '0';
            rs  <= R_RESP;
          end if;
        when R_RESP =>
          if m_rvalid = '1' then
            cdc_r_src_data <= m_rdata;
            cdc_r_src_resp <= m_rresp;
            rs <= R_FLAG;
          end if;
        when R_FLAG =>                     -- data a clock ahead of its flag
          cdc_rack <= '1';
          rs <= R_ACK;
        when R_ACK =>
          if rreq_m = '0' then
            cdc_rack <= '0';
            rs <= R_IDLE;
          end if;
      end case;

      if m_rst = '1' then
        wst      <= W_IDLE;
        rs       <= R_IDLE;
        awv      <= '0';
        wv       <= '0';
        arv      <= '0';
        cdc_wack <= '0';
        cdc_rack <= '0';
      end if;
    end if;
  end process;

end architecture;
