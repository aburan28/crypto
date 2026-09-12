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

// walker engines; ~10k LUTs each, the VU47P has 2.85M
`ifndef ECC_NENG
`define ECC_NENG 32
`endif

// walks per engine = 2**ECC_ID_W; must comfortably cover 2**(LOG_W+LOG_NB)
`ifndef ECC_ID_W
`define ECC_ID_W 8
`endif

// step unit: walks per Montgomery batch and batches in flight
`ifndef ECC_LOG_W
`define ECC_LOG_W 4
`endif
`ifndef ECC_LOG_NB
`define ECC_LOG_NB 3
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

`endif
