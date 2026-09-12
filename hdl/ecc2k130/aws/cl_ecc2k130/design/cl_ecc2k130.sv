// cl_ecc2k130.sv -- AWS F2 custom logic: the ECC2K-130 rho engine behind
// the shell's OCL AXI4-Lite port.
//
// Everything the host touches is ec2k_axil (VHDL, ../../../ec2k_axil.vhd):
// NENG walker engines, a load port and a distinguished-point queue, on a
// 32-bit register map at BAR0 of the application PF.  This file is the
// CL_TEMPLATE tie-off with that one block wired to OCL.  No DDR, no HBM,
// no DMA, no interrupts: the engine's traffic is a dozen register accesses
// per distinguished point, one every 2**26 steps per walk.
//
// Clocking.  The shell fixes clk_main_a0 at 250 MHz and the engine closes
// well above that, so the register block and the engines run on their own
// clock, clk_eng, made from clk_main_a0 by one MMCM (MMCM_ENG):
//
//   f(clk_eng) = 250 MHz * ECC_MMCM_MULT / ECC_MMCM_DIV
//
// with the VCO at 250 * ECC_MMCM_MULT MHz, which must stay within 800 to
// 1600.  The defaults (4, 3) give 333.3 MHz; build_afi.sh's CLK_MHZ picks
// the pair.  The OCL port stays on clk_main_a0 and crosses to clk_eng in
// BRIDGE (ec2k_axil_cdc, VHDL), whose synchronisers and data captures are
// the only paths between the two domains; cl_timing_user.xdc constrains
// them by name.  The engine reset is the shell reset or the MMCM losing
// lock, asserted asynchronously and released synchronously to clk_eng,
// then registered once more before it fans out.
//
// Geometry is set at build time through cl_ecc2k130_defines.vh
// (ECC_NENG, ECC_ID_W, ECC_DP_WEIGHT, ECC_MMCM_*); build_afi.sh passes
// them as -verilog_define so one source builds every image.

module cl_ecc2k130
    #(
      parameter EN_DDR = 0,
      parameter EN_HBM = 0
    )
    (
      `include "cl_ports.vh"
    );

`include "cl_id_defines.vh"
`include "cl_ecc2k130_defines.vh"

//=============================================================================
// GLOBALS
//=============================================================================

  always_comb begin
     cl_sh_flr_done    = 'b1;
     cl_sh_status0     = 'b0;
     cl_sh_status1     = 'b0;
     cl_sh_status2     = 'b0;
     cl_sh_id0         = `CL_SH_ID0;
     cl_sh_id1         = `CL_SH_ID1;
     cl_sh_status_vled = 'b0;
     cl_sh_dma_wr_full = 'b0;
     cl_sh_dma_rd_full = 'b0;
  end

//=============================================================================
// SHELL-SIDE RESET
//=============================================================================

  // rst_main_n is synchronous to clk_main_a0; registered twice for fanout,
  // active high for the VHDL.
  logic [1:0] rst_sh_pipe = 2'b11;
  logic       rst_sh;
  always_ff @(posedge clk_main_a0) rst_sh_pipe <= {rst_sh_pipe[0], ~rst_main_n};
  assign rst_sh = rst_sh_pipe[1];

//=============================================================================
// ENGINE CLOCK
//=============================================================================

  localparam real    ECC_MMCM_MULT_R = `ECC_MMCM_MULT;
  localparam real    ECC_MMCM_DIV_R  = `ECC_MMCM_DIV;
  localparam integer ECC_CLK_KHZ     = int'(250000.0 * ECC_MMCM_MULT_R / ECC_MMCM_DIV_R);

  logic clk_fb, clk_fb_buf, clk_eng_raw, clk_eng, mmcm_locked;

  MMCME4_ADV
    #(
      .BANDWIDTH          ("OPTIMIZED"),
      .CLKIN1_PERIOD      (4.000),
      .DIVCLK_DIVIDE      (1),
      .CLKFBOUT_MULT_F    (ECC_MMCM_MULT_R),
      .CLKFBOUT_PHASE     (0.0),
      .CLKOUT0_DIVIDE_F   (ECC_MMCM_DIV_R),
      .CLKOUT0_DUTY_CYCLE (0.5),
      .CLKOUT0_PHASE      (0.0),
      .COMPENSATION       ("AUTO"),
      .STARTUP_WAIT       ("FALSE")
      )
  MMCM_ENG
    (
      .CLKIN1       (clk_main_a0),
      .CLKIN2       (1'b0),
      .CLKINSEL     (1'b1),
      .CLKFBIN      (clk_fb_buf),
      .CLKFBOUT     (clk_fb),
      .CLKFBOUTB    (),
      .CLKOUT0      (clk_eng_raw),
      .CLKOUT0B     (),
      .CLKOUT1      (),
      .CLKOUT1B     (),
      .CLKOUT2      (),
      .CLKOUT2B     (),
      .CLKOUT3      (),
      .CLKOUT3B     (),
      .CLKOUT4      (),
      .CLKOUT5      (),
      .CLKOUT6      (),
      .DADDR        (7'd0),
      .DCLK         (1'b0),
      .DEN          (1'b0),
      .DI           (16'd0),
      .DWE          (1'b0),
      .DO           (),
      .DRDY         (),
      .PSCLK        (1'b0),
      .PSEN         (1'b0),
      .PSINCDEC     (1'b0),
      .PSDONE       (),
      .LOCKED       (mmcm_locked),
      .CLKINSTOPPED (),
      .CLKFBSTOPPED (),
      .PWRDWN       (1'b0),
      .RST          (rst_sh),
      .CDDCDONE     (),
      .CDDCREQ      (1'b0)
      );

  BUFG BUFG_FB  (.I(clk_fb),      .O(clk_fb_buf));
  BUFG BUFG_ENG (.I(clk_eng_raw), .O(clk_eng));

//=============================================================================
// ENGINE-SIDE RESET
//=============================================================================

  // Asserted at once by the shell reset or loss of lock, released after
  // four clk_eng edges of neither; then registered twice more for fanout.
  logic rst_eng_async, rst_eng_a;
  logic [1:0] rst_eng_pipe = 2'b11;
  logic       rst_eng;

  assign rst_eng_async = rst_sh | ~mmcm_locked;

  xpm_cdc_async_rst
    #(
      .DEST_SYNC_FF    (4),
      .INIT_SYNC_FF    (1),
      .RST_ACTIVE_HIGH (1)
      )
  RST_ENG_SYNC
    (
      .src_arst  (rst_eng_async),
      .dest_clk  (clk_eng),
      .dest_arst (rst_eng_a)
      );

  always_ff @(posedge clk_eng) rst_eng_pipe <= {rst_eng_pipe[0], rst_eng_a};
  assign rst_eng = rst_eng_pipe[1];

//=============================================================================
// OCL -> BRIDGE (clk_main_a0 -> clk_eng) -> ec2k_axil
//=============================================================================

  logic [31:0] e_awaddr, e_wdata, e_araddr, e_rdata;
  logic [3:0]  e_wstrb;
  logic [1:0]  e_bresp, e_rresp;
  logic        e_awvalid, e_awready, e_wvalid, e_wready, e_bvalid, e_bready;
  logic        e_arvalid, e_arready, e_rvalid, e_rready;

  ec2k_axil_cdc
    #(
      .SYNC_FF (3)
      )
  BRIDGE
    (
      .s_clk     (clk_main_a0),
      .s_rst     (rst_sh),
      .s_awaddr  (ocl_cl_awaddr),
      .s_awvalid (ocl_cl_awvalid),
      .s_awready (cl_ocl_awready),
      .s_wdata   (ocl_cl_wdata),
      .s_wstrb   (ocl_cl_wstrb),
      .s_wvalid  (ocl_cl_wvalid),
      .s_wready  (cl_ocl_wready),
      .s_bresp   (cl_ocl_bresp),
      .s_bvalid  (cl_ocl_bvalid),
      .s_bready  (ocl_cl_bready),
      .s_araddr  (ocl_cl_araddr),
      .s_arvalid (ocl_cl_arvalid),
      .s_arready (cl_ocl_arready),
      .s_rdata   (cl_ocl_rdata),
      .s_rresp   (cl_ocl_rresp),
      .s_rvalid  (cl_ocl_rvalid),
      .s_rready  (ocl_cl_rready),
      .m_clk     (clk_eng),
      .m_rst     (rst_eng),
      .m_awaddr  (e_awaddr),
      .m_awvalid (e_awvalid),
      .m_awready (e_awready),
      .m_wdata   (e_wdata),
      .m_wstrb   (e_wstrb),
      .m_wvalid  (e_wvalid),
      .m_wready  (e_wready),
      .m_bresp   (e_bresp),
      .m_bvalid  (e_bvalid),
      .m_bready  (e_bready),
      .m_araddr  (e_araddr),
      .m_arvalid (e_arvalid),
      .m_arready (e_arready),
      .m_rdata   (e_rdata),
      .m_rresp   (e_rresp),
      .m_rvalid  (e_rvalid),
      .m_rready  (e_rready)
      );

  ec2k_axil
    #(
      .NENG      (`ECC_NENG),
      .ID_W      (`ECC_ID_W),
      .LOG_W     (`ECC_LOG_W),
      .LOG_NB    (`ECC_LOG_NB),
      .FLUSH_CLK (`ECC_FLUSH_CLK),
      .CNT_W     (32),
      .DP_WEIGHT (`ECC_DP_WEIGHT),
      .DP_FIFO_W (`ECC_DP_FIFO_W),
      .CLK_KHZ   (ECC_CLK_KHZ)
      )
  ENGINE
    (
      .clk       (clk_eng),
      .rst       (rst_eng),
      .s_awaddr  (e_awaddr),
      .s_awvalid (e_awvalid),
      .s_awready (e_awready),
      .s_wdata   (e_wdata),
      .s_wstrb   (e_wstrb),
      .s_wvalid  (e_wvalid),
      .s_wready  (e_wready),
      .s_bresp   (e_bresp),
      .s_bvalid  (e_bvalid),
      .s_bready  (e_bready),
      .s_araddr  (e_araddr),
      .s_arvalid (e_arvalid),
      .s_arready (e_arready),
      .s_rdata   (e_rdata),
      .s_rresp   (e_rresp),
      .s_rvalid  (e_rvalid),
      .s_rready  (e_rready)
      );

//=============================================================================
// PCIM (unused)
//=============================================================================

  always_comb begin
    cl_sh_pcim_awaddr  = 'b0;
    cl_sh_pcim_awsize  = 'b0;
    cl_sh_pcim_awburst = 'b0;
    cl_sh_pcim_awvalid = 'b0;

    cl_sh_pcim_wdata   = 'b0;
    cl_sh_pcim_wstrb   = 'b0;
    cl_sh_pcim_wlast   = 'b0;
    cl_sh_pcim_wvalid  = 'b0;

    cl_sh_pcim_araddr  = 'b0;
    cl_sh_pcim_arsize  = 'b0;
    cl_sh_pcim_arburst = 'b0;
    cl_sh_pcim_arvalid = 'b0;
  end

  always_comb begin
    cl_sh_pcim_awid    = 'b0;
    cl_sh_pcim_awlen   = 'b0;
    cl_sh_pcim_awcache = 'b0;
    cl_sh_pcim_awlock  = 'b0;
    cl_sh_pcim_awprot  = 'b0;
    cl_sh_pcim_awqos   = 'b0;
    cl_sh_pcim_awuser  = 'b0;

    cl_sh_pcim_wid     = 'b0;
    cl_sh_pcim_wuser   = 'b0;

    cl_sh_pcim_arid    = 'b0;
    cl_sh_pcim_arlen   = 'b0;
    cl_sh_pcim_arcache = 'b0;
    cl_sh_pcim_arlock  = 'b0;
    cl_sh_pcim_arprot  = 'b0;
    cl_sh_pcim_arqos   = 'b0;
    cl_sh_pcim_aruser  = 'b0;

    cl_sh_pcim_rready  = 'b0;
  end

//=============================================================================
// PCIS (unused)
//=============================================================================

  always_comb begin
    cl_sh_dma_pcis_bresp   = 'b0;
    cl_sh_dma_pcis_rresp   = 'b0;
    cl_sh_dma_pcis_rvalid  = 'b0;
  end

  always_comb begin
    cl_sh_dma_pcis_awready = 'b0;

    cl_sh_dma_pcis_wready  = 'b0;

    cl_sh_dma_pcis_bid     = 'b0;
    cl_sh_dma_pcis_bvalid  = 'b0;

    cl_sh_dma_pcis_arready  = 'b0;

    cl_sh_dma_pcis_rid     = 'b0;
    cl_sh_dma_pcis_rdata   = 'b0;
    cl_sh_dma_pcis_rlast   = 'b0;
    cl_sh_dma_pcis_ruser   = 'b0;
  end

//=============================================================================
// SDA (unused)
//=============================================================================

  always_comb begin
    cl_sda_bresp   = 'b0;
    cl_sda_rresp   = 'b0;
    cl_sda_rvalid  = 'b0;
  end

  always_comb begin
    cl_sda_awready = 'b0;
    cl_sda_wready  = 'b0;

    cl_sda_bvalid = 'b0;

    cl_sda_arready = 'b0;

    cl_sda_rdata   = 'b0;
  end

//=============================================================================
// SH_DDR (present but disabled; the shell requires the instance)
//=============================================================================

   sh_ddr
     #(
       .DDR_PRESENT (EN_DDR)
       )
   SH_DDR
     (
      .clk                       (clk_main_a0 ),
      .rst_n                     (            ),
      .stat_clk                  (clk_main_a0 ),
      .stat_rst_n                (            ),
      .CLK_DIMM_DP               (CLK_DIMM_DP ),
      .CLK_DIMM_DN               (CLK_DIMM_DN ),
      .M_ACT_N                   (M_ACT_N     ),
      .M_MA                      (M_MA        ),
      .M_BA                      (M_BA        ),
      .M_BG                      (M_BG        ),
      .M_CKE                     (M_CKE       ),
      .M_ODT                     (M_ODT       ),
      .M_CS_N                    (M_CS_N      ),
      .M_CLK_DN                  (M_CLK_DN    ),
      .M_CLK_DP                  (M_CLK_DP    ),
      .M_PAR                     (M_PAR       ),
      .M_DQ                      (M_DQ        ),
      .M_ECC                     (M_ECC       ),
      .M_DQS_DP                  (M_DQS_DP    ),
      .M_DQS_DN                  (M_DQS_DN    ),
      .cl_RST_DIMM_N             (RST_DIMM_N  ),
      .cl_sh_ddr_axi_awid        (            ),
      .cl_sh_ddr_axi_awaddr      (            ),
      .cl_sh_ddr_axi_awlen       (            ),
      .cl_sh_ddr_axi_awsize      (            ),
      .cl_sh_ddr_axi_awvalid     (            ),
      .cl_sh_ddr_axi_awburst     (            ),
      .cl_sh_ddr_axi_awuser      (            ),
      .cl_sh_ddr_axi_awready     (            ),
      .cl_sh_ddr_axi_wdata       (            ),
      .cl_sh_ddr_axi_wstrb       (            ),
      .cl_sh_ddr_axi_wlast       (            ),
      .cl_sh_ddr_axi_wvalid      (            ),
      .cl_sh_ddr_axi_wready      (            ),
      .cl_sh_ddr_axi_bid         (            ),
      .cl_sh_ddr_axi_bresp       (            ),
      .cl_sh_ddr_axi_bvalid      (            ),
      .cl_sh_ddr_axi_bready      (            ),
      .cl_sh_ddr_axi_arid        (            ),
      .cl_sh_ddr_axi_araddr      (            ),
      .cl_sh_ddr_axi_arlen       (            ),
      .cl_sh_ddr_axi_arsize      (            ),
      .cl_sh_ddr_axi_arvalid     (            ),
      .cl_sh_ddr_axi_arburst     (            ),
      .cl_sh_ddr_axi_aruser      (            ),
      .cl_sh_ddr_axi_arready     (            ),
      .cl_sh_ddr_axi_rid         (            ),
      .cl_sh_ddr_axi_rdata       (            ),
      .cl_sh_ddr_axi_rresp       (            ),
      .cl_sh_ddr_axi_rlast       (            ),
      .cl_sh_ddr_axi_rvalid      (            ),
      .cl_sh_ddr_axi_rready      (            ),
      .sh_ddr_stat_bus_addr      (            ),
      .sh_ddr_stat_bus_wdata     (            ),
      .sh_ddr_stat_bus_wr        (            ),
      .sh_ddr_stat_bus_rd        (            ),
      .sh_ddr_stat_bus_ack       (            ),
      .sh_ddr_stat_bus_rdata     (            ),
      .ddr_sh_stat_int           (            ),
      .sh_cl_ddr_is_ready        (            )
      );

  always_comb begin
    cl_sh_ddr_stat_ack   = 'b0;
    cl_sh_ddr_stat_rdata = 'b0;
    cl_sh_ddr_stat_int   = 'b0;
  end

//=============================================================================
// USER-DEFINED INTERRUPTS (unused)
//=============================================================================

  always_comb begin
    cl_sh_apppf_irq_req = 'b0;
  end

//=============================================================================
// VIRTUAL JTAG (unused)
//=============================================================================

  always_comb begin
    tdo = 'b0;
  end

//=============================================================================
// HBM MONITOR IO (unused)
//=============================================================================

  always_comb begin
    hbm_apb_paddr_1   = 'b0;
    hbm_apb_pprot_1   = 'b0;
    hbm_apb_psel_1    = 'b0;
    hbm_apb_penable_1 = 'b0;
    hbm_apb_pwrite_1  = 'b0;
    hbm_apb_pwdata_1  = 'b0;
    hbm_apb_pstrb_1   = 'b0;
    hbm_apb_pready_1  = 'b0;
    hbm_apb_prdata_1  = 'b0;
    hbm_apb_pslverr_1 = 'b0;

    hbm_apb_paddr_0   = 'b0;
    hbm_apb_pprot_0   = 'b0;
    hbm_apb_psel_0    = 'b0;
    hbm_apb_penable_0 = 'b0;
    hbm_apb_pwrite_0  = 'b0;
    hbm_apb_pwdata_0  = 'b0;
    hbm_apb_pstrb_0   = 'b0;
    hbm_apb_pready_0  = 'b0;
    hbm_apb_prdata_0  = 'b0;
    hbm_apb_pslverr_0 = 'b0;
  end

//=============================================================================
// PCIE (unused)
//=============================================================================

  always_comb begin
    PCIE_EP_TXP    = 'b0;
    PCIE_EP_TXN    = 'b0;

    PCIE_RP_PERSTN = 'b0;
    PCIE_RP_TXP    = 'b0;
    PCIE_RP_TXN    = 'b0;
  end

endmodule // cl_ecc2k130
