// One ECC2K-130 Pollard rho core: walks circulate, one step retires per clock.
//
// The datapath blocks are generated from the same straight-line IR the CUDA
// client is built from (see ecc2k130/codegen/genverilog.py); this file is the
// sequencer around them, and it is hand written because it is control, not
// arithmetic.
//
//     R' = sigma^j(R) + R,    j = 3 + ((HW(x_R)/2) mod 8)
//
// Structure.  Every block has initiation interval 1 and no back pressure, so
// the "ring" the published FPGA designs describe needs no state machine at all:
//
//     load ->|                                                        |
//            +--> ecc_pre --> ecc_inv --> [state RAM] --> ecc_post ----+
//                   3 clk       96 clk        1 clk         24 clk
//
// A walk entering pre re-enters it RING_LAT clocks later, so exactly RING_LAT
// walks are resident and each advances one step per RING_LAT clocks.  The core
// retires one step per clock regardless, which is the number that matters.
//
// Why a RAM and not delay lines.  ecc_post needs x, y, d and e at the moment
// the inverse of d arrives, 96 clocks after ecc_pre produced them.  Holding
// them in registers would cost 4 x 131 x 96 = 50304 flip-flops per core, and
// flip-flops are the second-scarcest resource on the target part (see
// ecc2k130/FPGA-CEILING.md).  A RAM indexed by the walk's tag costs 2 block
// RAMs instead.
//
// The inverter is a port, not an instance.  One inverter is 8 multipliers, and
// Montgomery's trick means a whole device needs very few of them, so it is
// shared: see ecc2k130/FPGA-CEILING.md.  This core issues one inversion per
// step, which is the unbatched arrangement -- correct, and 10 multiplications
// per step instead of the 5.125 batching reaches.  ecc_batch.v is the
// optimisation, and this core is what it has to keep bit-exact.
//
// Distinguished points leave through dp*; the host supplies replacement start
// points through load*, exactly as in the published designs, because building a
// fresh start point costs 128 point additions and has no business in the ring.
`default_nettype none

module ecc_rho_core #(
    parameter M = 131,
    parameter TAGW = 8,
    parameter AW = 7,                      // state RAM address width
    parameter PRE_LAT = 3,
    parameter INV_LAT = 96,
    parameter POST_LAT = 24,
    parameter [7:0] DPCUT = 8'd34          // report when HW(x) <= DPCUT
) (
    input  wire            clk,

    // load / reload: inject a walk into the ring
    input  wire [M-1:0]    loadX,
    input  wire [M-1:0]    loadY,
    input  wire [TAGW-1:0] loadTag,
    input  wire            loadValid,

    // distinguished points out
    output wire [M-1:0]    dpX,
    output wire [TAGW-1:0] dpTag,
    output wire            dpValid,

    // every retired step, for verification and for checkpointing
    output wire [M-1:0]    stepX,
    output wire [M-1:0]    stepY,
    output wire [TAGW-1:0] stepTag,
    output wire            stepValid,

    // shared inverter, external
    output wire [M-1:0]    invA,
    output wire [TAGW-1:0] invTag,
    output wire            invValid,
    input  wire [M-1:0]    invR,
    input  wire [TAGW-1:0] invRTag,
    input  wire            invRValid
);

    localparam RING_LAT = PRE_LAT + INV_LAT + 1 + POST_LAT;

    // ---- ring input: a fresh walk, or the previous step's result ----------
    wire [M-1:0]    postX3;
    wire [M-1:0]    postY3;
    wire [TAGW-1:0] postTag;
    wire            postValid;

    wire [M-1:0]    ringX   = loadValid ? loadX   : postX3;
    wire [M-1:0]    ringY   = loadValid ? loadY   : postY3;
    wire [TAGW-1:0] ringTag = loadValid ? loadTag : postTag;
    wire            ringV   = loadValid ? 1'b1    : postValid;

    // ---- pass 1: weight, sigma^j, and the two addends (no multiplication) --
    wire [M-1:0]    preD;
    wire [M-1:0]    preE;
    wire [7:0]      preHw;
    wire [TAGW-1:0] preTag;
    wire            preV;

    ecc_pre131 uPre (
        .clk(clk), .x(ringX), .y(ringY),
        .tagIn(ringTag), .validIn(ringV),
        .d(preD), .e(preE), .hw(preHw),
        .tagOut(preTag), .validOut(preV));

    // x and y must reach the state RAM alongside d and e
    wire [M-1:0] preX;
    wire [M-1:0] preY;
    ecc_delay #(.WIDTH(M), .DEPTH(PRE_LAT)) uDelX (.clk(clk), .d(ringX), .q(preX));
    ecc_delay #(.WIDTH(M), .DEPTH(PRE_LAT)) uDelY (.clk(clk), .d(ringY), .q(preY));

    // ---- distinguished point test: one comparator ------------------------
    assign dpX     = preX;
    assign dpTag   = preTag;
    assign dpValid = preV && (preHw <= DPCUT);

    // ---- the inversion goes out to the shared inverter --------------------
    assign invA     = preD;
    assign invTag   = preTag;
    assign invValid = preV;

    // ---- in-flight state, indexed by the walk's tag -----------------------
    // Written when pass 1 finishes, read when the inverse comes back.
    localparam SW = 4 * M;
    reg [SW-1:0] state [0:(1<<AW)-1];
    reg [SW-1:0] stateQ;
    reg [TAGW-1:0] invRTagQ;
    reg            invRValidQ;
    reg [M-1:0]    invRQ;

    always @(posedge clk) begin
        if (preV)
            state[preTag[AW-1:0]] <= {preX, preY, preD, preE};
        stateQ     <= state[invRTag[AW-1:0]];
        invRQ      <= invR;
        invRTagQ   <= invRTag;
        invRValidQ <= invRValid;
    end

    wire [M-1:0] stX = stateQ[4*M-1:3*M];
    wire [M-1:0] stY = stateQ[3*M-1:2*M];
    wire [M-1:0] stD = stateQ[2*M-1:1*M];
    wire [M-1:0] stE = stateQ[1*M-1:0];

    // ---- pass 2: the affine addition (2 multiplications) ------------------
    ecc_post131 uPost (
        .clk(clk), .x(stX), .y(stY), .d(stD), .e(stE), .di(invRQ),
        .tagIn(invRTagQ), .validIn(invRValidQ),
        .x3(postX3), .y3(postY3),
        .tagOut(postTag), .validOut(postValid));

    assign stepX     = postX3;
    assign stepY     = postY3;
    assign stepTag   = postTag;
    assign stepValid = postValid;

endmodule
`default_nettype wire
