# Dedicated two-socket ecbench qualification protocol

Status: **preregistered; no host launched and no hardware measurements taken**.
This is a follow-on to [the counted SAT/PDP panel](../ecbench_sat_pdp_20261009/README.md)
and PR #1611. Its purpose is to establish whether the harness can admit
paired wall measurements on real NUMA hardware. It does not promote the
SAT row into a complete-cost ECDLP result.

## Frozen question and inputs

Hypothesis: on a dedicated, two-socket host, the native runner can reserve a
whole core, bind its child and anonymous pages to that core's NUMA node, and
obtain at least L2 on 90% of measured executions in each of two separate
node-local sessions. A/A strong-rho pairs should measure the host's wall
noise floor. A/B pairs use the same eight public targets, round seeds,
factor base, attempt caps, and interleaving as the source panel.

- The exact spec is [`spec.json`](spec.json), SHA-256
  `c5f2bcf67257569daecea0dc96876a61caee2cefb5d2894a63d0d7928b78b891`.
  It is the budget-qualified five-arm spec: 8 public targets on
  `icv1-f2m17-tm101-00378d4e`, one warm-up plus five measured rounds,
  strong signed-Frobenius rho, its A/A control, MITM, Frobenius MITM, and
  native-XOR SAT. The SAT solver has a 20,000-conflict per-call budget and
  no wall budget; all exits, timeouts and failed answers remain in the record.
- Candidate source is PR #1611, initially
  `a148d750f0c3da0088ec4568df48f580394cabe7`. Its binary hash,
  compiler, dependency resolution and clean source commit must be frozen
  before a session. A later fix changes the source hash and requires a
  dated protocol amendment before measurement, never an in-place change
  to an existing session.
- Reference floor: `S_floor = sqrt(pi/(2A))` from each record's curve facts.
  The measured reference is the strong-rho arm on the same targets and
  rounds. Report `S`, `S/S_floor`, `S/S_rho`, correctness, and unpriced units
  in one table. A ratio involving an unpriced row is diagnostic only.

## Host and placements

Proposed host: one AWS `c6i.metal` in `us-west-2`, Linux, using the existing
EC2 key-pair name `meow34`. AWS lists 128 vCPUs and 256 GiB for this type.
The host must show **two physical sockets and at least two NUMA memory nodes**
under `lscpu -e` and `numactl --hardware`; the instance type alone is not
proof. Record kernel, BIOS, microcode, clock governor, turbo state, SMT,
online and isolated CPUs, NUMA distances, IRQ affinities, PMU permissions,
and `ecbench host` in a sealed pre-run manifest.

Before any timed run, choose and freeze one whole non-CPU-0 core from each
socket, including every SMT sibling in its `--cpus` list. Keep the same
binary, spec, governor and turbo policy across the two sessions. Exclude
other jobs from both reserved cores and their siblings. Configure L3
settings where the host permits them (`isolcpus` or `nohz_full`, performance
governor, turbo known off); record every setting and blocker. Run the two
sessions separately, one with the CPU and memory local to NUMA node 0 and
one local to node 1. `ecbench`'s child binds memory to its run CPU's node;
verify `get_mempolicy` and `/proc/self/numa_maps` from every record. A
remote-memory experiment needs a distinct, explicit memory-placement mode
and is **not** inferred from these two local sessions.

The native command on each node is, after recording `CPU_LIST` in the
pre-run manifest:

```sh
RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 ./target/release/ecbench run \
  --spec research/ecbench_hardware_20261009/spec.json \
  --out research/ecbench_hardware_20261009/sessions/NODE_AND_RUN_ID \
  --cpus CPU_LIST --wait
```

Use a new output directory for every attempt. The harness takes the shared
lock, interleaves the arms, pins children, checks preflight pressure, and
records why a run did or did not earn L2/L3. Do not run a build, audit, or
other benchmark on the host while timed work runs. Do not use `--allow-busy`
to override a failed preflight for the admitted panel.

If a registered isolab worker with matching hardware becomes available,
generate a strict `isolab.job/v1` from the **same frozen binary and spec**
with `ecbench isolab-job --spec ... --binary ... --cpus 2 --policy strict`.
Keep the generated job JSON and worker receipt. The Python isolab service
is legacy orchestration; all ECDLP execution, accounting, replay and result
generation in this protocol use native `ecbench`.

## Admission and analysis

1. Freeze the binary SHA-256 and `ecbench host` capsule before the first
   session. Refuse a dirty build. Abort before measurement if the host is
   not two-socket/two-node, PMU access cannot be characterized, the
   benchmark core's siblings cannot be reserved, or the preflight is busy.
2. Run node 0 and node 1 with the same spec. Retain all planned and actual
   executions, including non-completions. Require every completed scalar to
   satisfy `[d]G = Q`; any incorrect answer stops promotion. A comparison
   uses only matched `(workload, round)` pairs and discloses unpaired runs.
3. Audit each session with `ecbench verify --replay-all --exit-code` and
   store the receipt hash. An auditor on a **different host class** must
   independently replay every measured run from the same source commit.
   The receipt's SHA-256 is the independent replay certificate.
4. For each node, require at least 90% of the 200 measured arm executions
   to earn L2. Wall comparisons use only admitted matched pairs and at
   least five pairs. Report the A/A ratio-of-sums, paired 95% interval,
   median and minimum, then each A/B ratio-of-sums and cluster-bootstrap
   interval. A claimed wall difference must exceed the observed A/A spread
   and have a 95% paired interval excluding no improvement. Do not pool
   the two nodes or mix levels to improve a headline.
5. Report user-space instructions and cycles **only** for arms whose PMU
   values are present on every admitted run. Include the coverage and
   errors otherwise. Compare the full cold resource vector: phase GAE,
   native solver conflicts/table probes, CPU time, solve and process wall,
   peak RSS, and online window. Do not convert SAT conflicts or table work
   into GAE without a separately measured and pinned conversion.
6. Stop and retain the failure if the source build cannot compile, replay
   differs, correctness fails, a per-call solver budget exits, or the L2
   threshold fails. L3 is an additional achieved grade, not an assumed
   consequence of bare metal. A/A noise above 5% blocks a runtime claim;
   retain the measured counts and blockers.

## Cost and access gate (checked 9 October 2026)

The AWS Price List API reports `US West (Oregon)` on-demand Linux
`c6i.metal` at **USD 5.44 per instance-hour**. A hard two-hour instance
limit would cap compute at **USD 10.88**, plus EBS and any data transfer.
The type is advertised in `us-west-2`, and the `meow34` key-pair name exists
there. The account's Standard On-Demand vCPU quota (`L-1216C47A`) is **0**
in `us-west-2`, `us-east-1`, `us-east-2`, `us-west-1` and `eu-west-1`.
Thus no instance can be launched now. Requesting a 128-vCPU quota increase
and launching this paid host require explicit owner approval. If approved,
enforce an external two-hour termination guard, verify the instance is
terminated, and retain the billing/instance-id receipt without secrets.
Never copy or print the private key. If quota or capacity remains unavailable,
leave the hardware measurements pending and do not substitute VM wall time.

## Closeout

Commit the pre-run host manifest, exact commands, both session directories,
all receipts and failed attempts, resource reports, A/A and A/B comparisons,
and the decision in a follow-on PR. Update the canonical scoreboard in that
same PR if a measurement changes its numbers or verdict. Until then the
hardware columns remain null.
