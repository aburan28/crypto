# Running hdl/ecc on AWS

Scripts to synthesise the secp256k1 datapath on an AWS build instance and
bring back real utilisation and timing numbers.

**Neither script has been run.** There was no AWS CLI and no Vivado in the
environment where they were written. The shell script is syntax-checked and
its argument handling and failure paths were exercised; the Tcl was reviewed
against the Vivado command set but not executed. Expect to fix something on
the first run, and read the cost warnings below before you start.

## Why synthesis and not an FPGA image

[`docs/ecc_fpga_cost_model.md`](../../../docs/ecc_fpga_cost_model.md)
estimates the multiplier at ~256 DSPs and ~300 MHz, and those two estimates
carry the entire FPGA-versus-GPU comparison. Vivado can settle both without
any FPGA being involved.

Actually *programming* an F2 instance is a different and much larger job.
`hdl/ecc` is a bare datapath: `clk`, `rst`, operands, results. It has no AWS
shell wrapper, no AXI4 interfaces, no DMA path, no host driver. The F2 shell
expects a custom-logic region that speaks its interfaces, so before an
Amazon FPGA Image could exist you would need to:

1. Wrap the datapath in the `aws-fpga` HDK custom-logic template, exposing
   an AXI4-Lite slave for control and an AXI4 master or stream for operands
   and distinguished-point output.
2. Build to a design checkpoint with the HDK build scripts, which run
   place-and-route against the shell's timing and pblock constraints.
3. Upload the tarball to S3 and call `aws ec2 create-fpga-image`, which
   takes hours and returns an AFI/AGFI pair.
4. Launch an F2 instance, `fpga-load-local-image` the AGFI, and write a
   host program against the FPGA management libraries.

That is a substantial project on its own, and none of it changes the
resource and timing numbers the cost model needs. Synthesis does.

## What the scripts do

| File | Purpose |
|---|---|
| `synth.tcl` | Out-of-context synthesis of `fp_mul_secp256k1` and `ec_add_pipe`, writing utilisation, hierarchical utilisation and timing reports plus a one-line summary |
| `run_aws_synthesis.sh` | Launches a build instance from the FPGA Developer AMI, copies the RTL up, runs `synth.tcl`, retrieves the reports, and terminates the instance |

Out-of-context means no I/O buffers and no board constraints, which is
right for a datapath destined to sit inside a shell. The numbers are
post-synthesis: utilisation is close to final for a DSP-bound design, but
treat the clock as an upper bound, since place-and-route on a part as large
as a VU47P typically costs another 10–20%.

## Cost, and how the script protects you

It launches an on-demand EC2 instance. A `z1d.2xlarge` is roughly $0.75 an
hour in `us-east-1`; this design should synthesise in well under an hour.

- Nothing launches until you type `y` at a prompt that shows the region,
  instance type, AMI and estimated rate.
- `--dry-run` prints that plan and exits without launching.
- The instance is terminated on **every** exit path — success, error,
  `Ctrl-C` — by a trap armed the moment the instance exists. A leaked
  instance is the expensive failure here, so the trap is the one piece of
  the script worth reading before you run it.
- `--keep` disables that, for debugging a failed synthesis. It prints the
  instance id and the exact command to terminate it. Use it deliberately.

Verify with `aws ec2 describe-instances` afterwards regardless. A script
promising to clean up is not the same as it having cleaned up.

## Use

```bash
# see the plan, spend nothing
./run_aws_synthesis.sh --key-name mykey --key-file ~/.ssh/mykey.pem --dry-run

# real run
./run_aws_synthesis.sh --key-name mykey --key-file ~/.ssh/mykey.pem \
                       --subnet subnet-0123 --sg sg-0123

# sweep the clock to find where timing actually closes
for p in 4.0 3.333 2.5 2.0; do
  ./run_aws_synthesis.sh --key-name mykey --key-file ~/.ssh/mykey.pem \
                         --period $p --out reports-$p
done
```

Prerequisites: awscli v2 configured, an EC2 key pair you hold the private
key for, a subnet with public IP and a security group allowing SSH from
your address, and a one-time free subscription to the AWS FPGA Developer
AMI in the Marketplace. The AMI is used because it ships Vivado already
licensed for AWS FPGA parts.

## What to do with the results

`reports/summary.txt` gives LUT, FF, DSP, BRAM, worst negative slack and an
implied Fmax per module. Compare against the cost model:

| Quantity | Cost model says | Source |
|---|---|---|
| DSP48E2 per multiplier | 256 | derived: 64 products of 32×32, 4 DSPs each |
| LUT per multiplier | ~10k | estimated |
| Clock | ~300 MHz | estimated from a 296-bit adder carry chain |

If the measured clock lands well below 300 MHz, the fix is already
identified in [`../README.md`](../README.md): split each row adder into two
pipeline stages, giving 16 row stages instead of 8. Throughput is
unaffected because the initiation interval stays at 1 and the engine
already tolerates latency by interleaving more walks.

Update the cost model with whatever comes back, and relabel those rows from
estimated to measured.
