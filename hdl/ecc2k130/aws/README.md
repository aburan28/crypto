# Synthesising hdl/ecc2k130 on AWS

`synth.tcl` synthesises `gf131_mul`, `ec2k_batch_pipe`, `ec2k_step_pipe`
and `ec2k_walker` out-of-context for the VU47P and writes utilisation,
hierarchical utilisation and timing reports plus a one-line summary per
module. The launch script lives with the secp256k1 datapath and works here
unchanged:

```bash
# see the plan, spend nothing
HDLDIR=$(pwd)/.. ../../ecc/aws/run_aws_synthesis.sh \
    --key-name mykey --key-file ~/.ssh/mykey.pem --dry-run

# real run, sweeping the clock
for p in 3.333 2.5 2.0; do
  HDLDIR=$(pwd)/.. ../../ecc/aws/run_aws_synthesis.sh \
      --key-name mykey --key-file ~/.ssh/mykey.pem --period $p --out reports-$p
done
```

Read [`../../ecc/aws/README.md`](../../ecc/aws/README.md) first: it
explains the cost, the termination trap, and why this is synthesis rather
than an FPGA image. Everything there applies. **Nothing here has been run**;
there was no Vivado in the environment this was written in.

## What to compare against

The estimates in [`../README.md`](../README.md) that the results replace:

| Quantity | Estimate | Basis |
|---|---|---|
| LUT per multiplier | 5–6k | 9801 AND folded three per LUT6 with their XORs, plus the remaining ~9k XOR2 at five per LUT6 |
| FF per multiplier | ~3k | eight stages of 131–600 bits |
| DSP per multiplier | 0 | binary-field arithmetic has no carries |
| LUT per batched step unit | ~10k | multiplier, 640 × 131 bits of distributed RAM with two read ports on two of the arrays, `sigma^j` muxes on the operand path, two split popcounts |
| BRAM | 0 | every array is small enough for LUTRAM and has one write port |
| Clock | 300–400 MHz | leaf stage of the multiplier (33 AND terms per bit) and the operand-forming stage of the step unit are the deepest, three to four LUT levels |

If the clock lands low, sweep `MUL_KARATSUBA` in `gf131_pkg.vhd`: three
levels give a 17-bit leaf, halving the deepest stage, for two more clocks
of latency, which the scheduler absorbs by holding more batches in flight.
If the LUT count lands high, the same sweep in the other direction tells
whether Karatsuba's extra XORs are paying for themselves on this fabric.
The step unit's retire path is a RAM write and one XOR; it should not be
the limit.
