// cl_ecc2k130_defines.vh -- build-time geometry of the ECC2K-130 image.
//
// Each define yields to a -verilog_define on the synth_design command line,
// which is how build_afi.sh sets them (ECC_NENG=... etc.), so these are the
// defaults for a build that passes nothing.
//
// The host reads the geometry back from the GEOM register; nothing here
// needs to be repeated on the host side except the dp weight, which the
// campaign's --dp-weight must equal.

`ifndef CL_ECC2K130_DEFINES
`define CL_ECC2K130_DEFINES

// walker engines; ~9k LUTs and 22 RAMB36 each as synthesised, the VU47P
// has 1.30M LUTs and 2016 RAMB36
`ifndef ECC_NENG
`define ECC_NENG 48
`endif

// walks per engine = 2**ECC_ID_W; must comfortably cover 2**(LOG_W+LOG_NB),
// and 512 fills the FIFO's block RAM column exactly
`ifndef ECC_ID_W
`define ECC_ID_W 9
`endif

// step unit: walks per Montgomery batch and batches in flight (16 x 16
// fills the 512-deep block RAMs; 5.29 clocks per step against 5.54 at 8)
`ifndef ECC_LOG_W
`define ECC_LOG_W 4
`endif
`ifndef ECC_LOG_NB
`define ECC_LOG_NB 4
`endif

// idle clocks before a partial batch is padded and issued
`ifndef ECC_FLUSH_CLK
`define ECC_FLUSH_CLK 32
`endif

// distinguished point: normal-basis weight of x <= this (the challenge's 34)
`ifndef ECC_DP_WEIGHT
`define ECC_DP_WEIGHT 34
`endif

// report queue depth = 2**ECC_DP_FIFO_W
`ifndef ECC_DP_FIFO_W
`define ECC_DP_FIFO_W 6
`endif

// engine clock = 250 MHz * ECC_MMCM_MULT / ECC_MMCM_DIV, VCO = 250 * MULT
// within 800..1600 MHz; both may be fractional in eighths (4, 3 = 333.3 MHz;
// 4, 4 = 250; 6, 5 = 300; 6, 4 = 375; 4, 2.5 = 400).  build_afi.sh's
// CLK_MHZ sets them.
`ifndef ECC_MMCM_MULT
`define ECC_MMCM_MULT 4
`endif
`ifndef ECC_MMCM_DIV
`define ECC_MMCM_DIV 3
`endif

`endif
