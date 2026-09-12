# Synthesising hdl/ecc2k130 on AWS

`synth.tcl` synthesises `gf131_mul`, `ec2k_step_pipe` and `ec2k_walker`
out-of-context for the VU47P and writes utilisation, hierarchical
utilisation and timing reports plus a one-line summary per module. The
launch script lives with the secp256k1 datapath and works here unchanged:

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
| LUT per multiplier | ~9k | 17161 AND folded three per LUT6 with their XORs, plus 22k XOR2 at five per LUT6 |
| FF per multiplier | ~5k | eight rows of 261 + 131 + 136 bits, plus the ends |
| DSP per multiplier | 0 | binary-field arithmetic has no carries |
| LUT per step unit | ~12k | multiplier plus slot memories, sigma multiplexers, two popcounts |
| Clock | 300–400 MHz | row stages are two to three LUT levels; the step unit's retire path is the likely limit |

If the clock lands low, look at the step unit before the multiplier: the
retire path reads a slot memory, permutes through `gf_sigma_j`, XORs and
writes back in one clock. Splitting it costs one more clock of step latency
and nothing in throughput, since the scheduler already tolerates latency by
holding more slots in flight.
