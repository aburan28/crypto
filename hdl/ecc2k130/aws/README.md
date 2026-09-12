# hdl/ecc2k130 on AWS: synthesis, the F2 image, and an FPGA fleet

Two things live here. `synth.tcl` is the out-of-context synthesis probe
that brings back LUT, FF and Fmax for the engine's modules, spending nothing
on an FPGA. Everything else turns the engine into an F2 image and runs it
as part of the campaign in [`ecc2k130/aws/`](../../../ecc2k130/aws/):
the same bucket, the same slot registry, the same `dp/` corpus, the same
`merge.py`, with F2 instances and GPU instances side by side.

**Nothing here has been run.** There was no Vivado, no AWS account and no
F2 instance in the environment this was written in. The VHDL and the host
program are tested (GHDL, and the host against its own model of the
register block); the CL wrapper parses; the shell scripts pass `bash -n`.
Expect the first build to need a fix or two in the HDK glue, and read the
build log.

| File | What it is |
|---|---|
| `synth.tcl` | Out-of-context synthesis of `gf131_mul`, `ec2k_batch_pipe`, `ec2k_step_pipe`, `ec2k_walker`, `ec2k_axil` for the VU47P; run through `ecc/aws/run_aws_synthesis.sh` |
| `cl_ecc2k130/` | The F2 custom logic: `design/cl_ecc2k130.sv` ties off every shell interface except OCL and instantiates `ec2k_axil` there; `build/scripts/encrypt.tcl` and `synth_cl_ecc2k130.tcl` are the HDK's hooks, taught to read VHDL |
| `build_afi.sh` | From your laptop: push sources, launch a build instance, follow it, promote the resulting AFI |
| `build_afi_instance.sh` | What the build instance runs: HDK setup, `aws_build_dcp_from_cl.py`, upload, `create-fpga-image` |
| `bootstrap_f2.sh` | User data of an F2 worker: SDK, host program, load the AFI into every slot, one `worker.py` per slot |
| `f2.sh` | Launch template and EC2 Fleet for F2 instances, measured in FPGAs |

## The flow

Prerequisites: `ecc2k130/aws/infra.sh` has run (bucket, worker role,
security group, `campaign.json`), and the account is subscribed to the free
**FPGA Developer AMI** in the Marketplace (Ubuntu, Vivado preinstalled and
licensed for F2 builds).

```bash
cd hdl/ecc2k130/aws

./build_afi.sh push                     # engine + host sources -> s3://BUCKET/fpga/source.tar.gz
NENG=32 ./build_afi.sh launch           # r6i.4xlarge with the FPGA Developer AMI; hours
./build_afi.sh status                   # instance state, log tail, AFI state
./build_afi.sh wait                     # blocks until the AFI is available
./build_afi.sh promote                  # -> s3://BUCKET/fpga/afi.json, the image workers load

./f2.sh infra                           # launch template from bootstrap_f2.sh
./f2.sh up 8 --on-demand 8              # 8 FPGAs; f2.6xlarge/12xlarge/48xlarge, cheapest per AZ
./f2.sh status
../../../ecc2k130/aws/merge.py          # the same merge as for the GPUs; the points are indistinguishable
./f2.sh down
```

`launch` starts the instance with `--instance-initiated-shutdown-behavior
terminate` and the build script ends with `shutdown`, so a finished or
failed build costs nothing after it stops; `KEEP=1` keeps it for a
post-mortem over SSM. The build log is copied to
`fpga/builds/TAG/build.log` every five minutes, the utilisation and timing
reports to `fpga/builds/TAG/reports/`, and `afi.json` records the AFI and
AGFI ids with the geometry the image was built for.

## Geometry

The image is parameterised at build time through `-verilog_define`s that
`build_afi.sh` passes as `NENG`, `ID_W`, `DP_WEIGHT`
(`cl_ecc2k130_defines.vh` has the defaults). The host reads them back from
`GEOM` and refuses to run against an image with the wrong distinguished-point
weight, and `bootstrap_f2.sh` refuses to start if `campaign.json`'s
`dpWeight` differs from the image's, so a stale image cannot poison the
corpus.

`NENG` is the number to sweep. Each engine is one batched step unit plus
its walk memory, about 8.4k LUTs, 20 RAMB36 and 2 RAMB18 as synthesised
(`../README.md`, "Capacity"); the VU47P has 1.30M LUTs and 2 016 RAMB36.
The default 48 is a third of the device; 64 is 41% of the LUTs and 67% of
the block RAM, and the RAM is what runs out first, around 80. Read
`synth_utilization` and the post-route timing from the reports, then go to
what fits.

`ID_W` sets walks per engine, `2^ID_W`. Each walk is 304 bits of block
RAM, and the step unit holds `W · 2^LOG_NB` = 256 walks at once; 512 (the
default) keeps a full batch forming while reports and reloads drain and
fills the FIFO's block RAMs exactly, so fewer walks would save nothing.
More walks means more work lost on a restart and a longer time to the
first report, nothing else.

## Clocking

The F2 shell fixes `clk_main_a0` at 250 MHz, and the engine synthesises
with +1.98 ns of slack at 4.0 ns, so the CL makes its own clock: one
`MMCME4_ADV` (`MMCM_ENG` in `cl_ecc2k130.sv`) multiplies `clk_main_a0` by
`ECC_MMCM_MULT / ECC_MMCM_DIV` into `clk_eng`, and the register block and
every engine run on that. `build_afi.sh` takes `CLK_MHZ` (250, 300, 333 —
the default — 350, 375, 400) or the raw `MMCM_MULT` / `MMCM_DIV` pair;
the VCO at `250 · MULT` MHz must stay within 800–1600, and both values may
be fractional in eighths. The image reports the frequency it was built
for in the `CLOCK` register and the host prints it.

The shell's OCL port stays on `clk_main_a0` and crosses to `clk_eng` in
`BRIDGE`, an `ec2k_axil_cdc` (`../ec2k_axil_cdc.vhd`): one transaction of
each kind in flight, level handshakes through `ASYNC_REG` synchronisers,
data written a clock before its flag and copied a synchroniser after it.
Those synchroniser inputs and data copies are the only paths between the
two clocks, and `cl_timing_user.xdc` constrains exactly them by name
(false path into `cdc_*_meta`, `set_max_delay -datapath_only` from
`cdc_*_src` to `cdc_*_cap`). Because both clocks come from the same
source, Vivado times anything else that crosses against their edge
relationship, so a path that was forgotten fails timing rather than
sneaking through. The engine reset is the shell reset or the MMCM losing
lock, through an `xpm_cdc_async_rst` into `clk_eng` and one more register
per engine inside `ec2k_axil`, so no single flop drives every engine.

The MMCM and bridge cost 250 LUTs. Above what the routed design closes,
the build is flagged `timingViolated` in `afi.json`; such an image runs
(the host verifies every report by rewalking it, so a wrong step costs
throughput, not corpus integrity) but should be rebuilt one table entry
lower. The 250 MHz entry (`MMCM_MULT=4 MMCM_DIV=4`) is the baseline with
the same structure.

## What each F2 instance does

`bootstrap_f2.sh` clones `aws-fpga` (`f2` branch) and runs `sdk_setup.sh`
for the management tools and the `fpga_mgmt` library, fetches
`campaign.json`, `fpga/afi.json`, `fpga/source.tar.gz` and `aws/worker.py`
from the bucket, builds `ecc2k130-fpga` with `make pci`, and runs its
`--selftest`. Then per slot: `fpga-load-local-image -S slot -I agfi -R`,
wait for `loaded`, and probe the slot with `ecc2k130-fpga --launches 1`,
which checks `MAGIC` and `GEOM`. Only a slot that answers as this image
gets a worker.

The workers are `ecc2k130-worker@SLOT` systemd units running the
unmodified `worker.py` with `ECC_CLIENT=ecc2k130-fpga` and `ECC_GPU=SLOT`;
`worker.py` claims a campaign slot, hands the client `--run-id`,
`--dp-file`, `--checkpoint` and the rest, uploads `dp.bin` and the
checkpoint on its schedule, retires the slot on exit 6, and reports the
device as `ECC_DEVICE_NAME` (`f2.6xlarge fpga 32x256 TAG`) in the slot
registry. A spot interruption stops the workers with SIGTERM, the client
flushes its points and writes its checkpoint, and the replacement instance
resumes from S3; the walks that were in the FPGA at that moment are lost,
which costs at most `2^dpWeight` steps each.

## Cost, roughly

An F2 build instance runs for hours at r6i.4xlarge prices; the AFI
generation after it is free. `f2.6xlarge` on-demand is a few dollars an
hour for one VU47P; if the engine reaches the estimated 11–19 G steps/s
that is a small fraction of the GPU fleet's cost per step
(`ecc2k130/RTX-PRO6000.md` has the GPU side). None of that is measured
until the first build.

## Out-of-context synthesis without an image

`synth.tcl` synthesises the modules out-of-context and writes utilisation,
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
than an FPGA image.

### What came back

The estimates in [`../README.md`](../README.md) against what Vivado 2025.2
synthesised (out of context, `xcvu47p-fsvh2892-2-e`, 4.0 ns clock):

| Quantity | Estimate | Synthesised |
|---|---|---|
| LUT per multiplier | 5–6k | 5 547 at two Karatsuba levels, **4 855 at three** (now the default), 5 019 at four |
| FF per multiplier | ~3k | 2 892 / 4 647 / 7 262 at two / three / four levels |
| DSP per multiplier | 0 | 0 |
| LUT per engine (step unit + walker) | ~10k | 13 106 with the memories in LUTRAM (4 576 of them); **8 100 – 8 410** with them in block RAM, 180 LUTRAM left; 7 404 FF |
| Register block (`ec2k_axil`) | — | 647 LUTs, 2 032 FF for two engines, of which a spine stage of ~300 LUTs, ~870 FF per engine; bridge 228 LUTs, 201 FF |
| BRAM | 0 | **20 RAMB36 + 2 RAMB18 per engine** since the memories moved out of LUTRAM (below) |
| Clock | 300–400 MHz | +1.98 ns slack at 4.0 ns, +1.08 ns at 3.0 ns (synthesis, two engines with register block and bridge); the shell fixes `clk_main_a0` at 250 MHz, so the CL's own MMCM makes the engine clock, 333 MHz by default |

The first CL synthesis (32 engines) read 1.29M LUTs, four times this: a
two-writer counter array in the walker had become 8k flip-flops behind a
256:1 mux, and every distinct `array(expr)` in the step unit's operand
stage had become its own LUTRAM copy. Both are fixed
([`../README.md`](../README.md), "Capacity"); the per-module probe that
found them is a plain out-of-context `synth_design` of one entity with a
LUT histogram by driven signal, which is the tool to reach for if a future
build's numbers do not add up.

The first 48-engine *implementation* (250 MHz, the flat register block of
that revision) placed at **−2.5 ns** on the 4 ns clock while the engine
alone synthesised at +2 ns: the placer's log shows it replicating the
load and run registers across SLRs and re-placing the report holding
registers behind the 48:1 mux, i.e. the interconnect between the register
block and engines spread over three SLRs, not the arithmetic. That is what
the spine in `ec2k_axil` is for ([`../README.md`](../README.md), "The
host interface"): every inter-engine wire is register to neighbour, and
the block has no path that touches more than one engine.

The same log, read to the end, showed a second failure inside every
engine, and the larger one: total negative slack was −246 µs over 18 000
endpoints, "long congestion 32×32" in two SLRs, and the nets that could
not be fixed were the write address and write enable of each LUTRAM array
(`f_wr`, `leaf_addr`, `fb`, the multiplier tag's address decode — 1 400
to 1 800 loads each), which physical optimisation reports it cannot
replicate because the loads are memory primitives. 4 800 LUTRAM per
engine, 48 times, spread across every SLICEM of a region is a design
that cannot place tightly. Every wide memory is now block RAM
([`../README.md`](../README.md), "Capacity"), with the read latency of
two hidden by a look-ahead tag on the retire side and a prefetch buffer
in the walker; the engine lost 4 700 LUTs and gained 21 block RAMs,
and the free depth of a block RAM made 512 walks and 16 batches in flight
the default, 5.29 clocks per step against 5.54.
