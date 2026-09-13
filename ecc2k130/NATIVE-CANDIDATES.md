# Native arithmetic candidates: compiled, throughput unmeasured

The objective is to improve on the existing 14.637530 B scalar updates/s
benchmark and 14.106673 B/s DP34 collection audit. Those are measured reference
values from [SHARED-SIGMA.md](SHARED-SIGMA.md), not a fundamental 14 B/s ceiling
and not matched controls for these candidates. Neither new flag is promoted.

The direct AWS attempt authenticated successfully and passed EC2 launch permission
checks. GPU allocation remains blocked by **InsufficientInstanceCapacity**:
`g7e.2xlarge` failed in all four us-west-2 zones, and automatic placement failed
for both permitted single-GPU sizes. An earlier automatic-placement request
timed out; its exact client token was checked for an unacknowledged instance.
The [capacity attempt receipt](benchmarks/native-candidates/aws-capacity-attempt.json)
records the requests and final instance check. No GPU timings were obtained.

Earlier Actions runs
[34762622979](https://github.com/aburan28/crypto/actions/runs/34762622979) and
[34762827266](https://github.com/aburan28/crypto/actions/runs/34762827266)
failed authentication using repository secrets. Their
[receipt](benchmarks/native-candidates/aws-attempt.json) is historical; the
subsequent direct credentials worked. Neither candidate is promoted.

This is **engineering**. The generic-group work boundary and walk rules do not
change. At batch 16, every variant still uses `5 + 5/16 = 5.3125` field products
per scalar update, a field-product cost ratio of **1.0** to the control. There
is no claim of a discrete-logarithm exponent improvement. Wall-clock speed is
the practical objective; it remains unmeasured for the new variants.

## Two independent changes

`PACKED_CLMAD_SQUARE=1` implements [native squaring](NATIVE-SQUARE.md).

`PACKED_KARAT3=1` replaces the bit-mask 128-by-3-bit multiplication tail with
three-limb Karatsuba recombination. For 64-bit limbs `A0,A1` and three-bit `A2`,
write `A=A0+A1*t+A2*t^2`, with `t=X^64`, and likewise for `B`. Set
`Di=Ai*Bi` and `Mij=(Ai+Aj)*(Bi+Bj)` in characteristic two. The five limb
coefficients of the product are:

```
D0,
M01 + D0 + D1,
D1 + M02 + D0 + D2,
M12 + D1 + D2,
D2.
```

There are five full 64-bit carryless products and one tiny three-bit product.
This adds four CLMAD instructions per field product compared with the control,
but removes the longer mask/shift tail. The top cross term has degree at most
65, so its fourth 32-bit word cancels to zero. Inputs retain their existing
canonical three-bit top-limb contract. Reduction and paired-product semantics
are unchanged. The flag requires CLMAD and defaults to zero; Make, Modal image
and rebuild caches, metadata, runtime identity gates and GPU probes carry it.

## Compiler evidence

[compiler-review.json](benchmarks/native-candidates/compiler-review.json)
retains the CUDA 13.3.73 sm_120 build commands, device-source and binary hashes,
compiler output, resource reports and opcode counts. Four complete clients and
their arithmetic, storage and shared-sigma probes compiled successfully. All
use B16/T256/min2/CLMAD1/COMPACT1/WP2/TILE256/SHARED1.
The receipt names the exact source revision and the trailing-newline-only
adjustments made by local reconstruction; these reproduce its input hashes.
Client compiler output is retained for the packed walk, init and their helpers.

| Variant | Class | Static walk instructions | Ratio to control | Non-NOP instructions | Walk registers | Walk stack/local bytes | Correctness status |
|---|---|---:|---:|---:|---:|---:|---|
| Control | engineering | 4,513 | 1.0000 | 4,448 | 104 | 0/0 | upstream GPU baseline; this build unrun |
| Native square | engineering | 4,003 | 0.8870 | 3,885 | 104 | 0/0 | source proof and host checks; GPU pending |
| Three-limb Karatsuba | engineering | 4,403 | 0.9756 | 4,279 | 94 | 0/0 | field and raw-product host checks; GPU pending |
| Both | engineering | 3,889 | 0.8617 | 3,716 | 94 | 0/0 | CUDA probes compile; GPU pending |

These counts cover the static walk text section, including each out-of-line
helper once. They are **not dynamic instructions per iteration**, and the ratio
is not a throughput prediction. NOP scheduling slots are included in the first
count. All four walk builds report zero spill stores/loads and a 2,816-byte
compiled shared extent, including the 1,024-byte device reservation. Fewer
registers do not by themselves establish higher occupancy or faster execution.

`make test-karat3` passes the existing independent field reference checks and
27,170 unreduced polynomial products against bit convolution: all 131-by-131
basis pairs, 10,000 dense cases covering every pair of three-bit top limbs,
and nine edge pairs, with output canaries. The native-square Boolean proof
still covers all 2^32 inputs. These host checks do not execute CLMAD.

## Reproduce and measure

With CUDA 13.3.73 installed, from `ecc2k130/`:

```sh
make check-cli test-clmad-square test-karat3
python3 codegen/native_candidate_bench.py --compile-only --out /tmp/native-compile
python3 codegen/native_candidate_bench.py --out /tmp/native-gpu
```

Use a new output directory for each invocation. Full mode requires exactly one
RTX PRO 6000 Blackwell Server Edition. It binds all clients to that GPU's UUID,
disables PTX JIT, retains hashes and raw output, and executes all device probes
and client DP replay/restart/resume/guard checks before timing. It also compares
normal-basis checkpoints across all binaries at 8, 128 and 257 workers.

Warmups are excluded. Each candidate receives three alternating pairs with its
control for each of the complete scalar benchmark and DP34 collection. Every
sample must complete **201,863,462,912 scalar updates** at the unchanged
385,024-worker geometry. Missing or wrong mode markers, incomplete counts,
mismatches, dropped reports, nonfinite rates and differing sorted full DP
multisets invalidate the run. Duplicates remain in the corpus comparison.

The predeclared acceptance rule is at least 1% improvement in the ratio of
medians for **both** workloads, with every pair favoring that candidate.
Borderline results need a fresh allocation. Compilation alone cannot satisfy
the rule. Defaults remain off unless a candidate meets the GPU gates.

The dedicated `ECC2K-130 native candidate benchmark` Actions workflow attempts
one `g7e.2xlarge` using the existing AWS credentials and worker infrastructure.
It first performs an EC2 permission dry run, then uses a unique benchmark
source/result prefix and its own startup script. It never starts a campaign
worker or edits the existing template, IAM configuration or fleet. The parent
terminates its exact instance on completion/failure, and a separate 55-minute
instance shutdown deadline bounds an interrupted parent. The Docker workload
has a 45-minute timeout. No GPU execution is implied merely by adding this
workflow; its recorded result determines whether validation succeeded.

The first workflow run failed authentication before reaching EC2; it did not
exercise the EC2 startup, device gates, timings or termination path. Those paths
have offline checks and remain pending live validation.

The index-calculus scoreboard has no new throughput or operation-count result
to add: GPU timing is pending and the operation-count ratio is unchanged.

## Existing templates without instance roles

For a G7e template that already supplies its AMI, subnet security group and
100-GiB disposable root disk, use `--launch-template NAME --presigned-transfer`
with `aws/native_benchmark.py`. The controller uses its normal AWS credential
chain. It signs one-hour URLs for exactly three objects under the unique
benchmark prefix: source GET, results PUT and completion-status PUT. The VM
receives these limited transfer capabilities and no long-lived AWS credentials
or instance profile. The controller still checks launch permissions first,
retains results and terminates its exact instance. URLs are not placed in the
repository or result receipt. The existing instance-profile path remains the
default. Both startup variants pass offline shell and isolation checks.

## Capacity placement options

`--availability-zone us-west-2b` selects an existing eligible default subnet.
Alternatively, `--automatic-placement` lets EC2 select a zone, after verifying
that the existing VPC is default and its default subnets enable public IPs.
These options are mutually exclusive. `--instance-type g7e.4xlarge` is the only
alternative to the default `g7e.2xlarge`; AWS hardware metadata must confirm
exactly one GPU. Both sizes use the same GPU comparison geometry. Eight offline
benchmark and launch checks pass; startup, device execution and timing remain
pending a successful allocation.

## Cross-region measurement retry

A subsequent retry attempted automatic placement in three regions. The
existing us-west-2 single-GPU setup and us-east-1 `g7e.2xlarge` returned
`InsufficientInstanceCapacity`. us-east-2 `g7e.2xlarge` failed because the
regional instance-bucket vCPU quota was zero. us-east-1 `g7e.4xlarge` exceeded
the available quota within its 32-vCPU limit. These are allocation failures,
not failed device correctness tests; no GPU timing was obtained. The
[retry receipt](benchmarks/native-candidates/aws-cross-region-attempt.json)
records the final exact-token checks.

For a region without a worker template, `--launch-config PATH` accepts an
explicit JSON object containing `ImageId`, `SecurityGroupIds` and
`BlockDeviceMappings`, plus `IamInstanceProfile` when not using presigned
transfer. The existing one-GPU and disposable-disk checks still apply.
`--s3-region us-west-2` keeps the approved source/results bucket in its original
region while EC2 runs elsewhere. The retry selected Amazon-owned Deep Learning
Base OSS Nvidia Driver GPU AMIs (Ubuntu 24.04), an existing default security
group with outbound connectivity, and a disposable 100-GiB root disk. No IAM,
network, quota or production-template changes were made.
