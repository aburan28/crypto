# Stage 22: licensed Magma replay of the balanced Phase-B panel

Stage 22 prepares the missing matched Magma arm for the 160 public synthetic
Phase-B point-decomposition inputs. It does not execute Magma in CI, does not
contain a private target, and does not change the conclusion about Koblitz
index calculus. A prepared packet is execution readiness, not a result.

## Frozen source and information boundary

The packet builder verifies the complete successor-04 evidence bundle and the
original terminal run before creating output. It uses the archived
`inputs/phase-b-protocol.json` bytes from successor-04; the later top-level
Phase-B protocol is not a substitute. The live source run must have the exact
terminal seal SHA-256
`e528f9bbaa6b7554224379cb32f7b1ed5be0a3cb8bde15e56dee7d1a5f75b538`
and the exact 4,797-file inventory commitment
`bb4377054700edf3049a7590b9b9fdd8fe6a0c57b6a216a7928dd211dc74a12c`.

Every blind bundle entry and complete four-file source instance is checked by
the existing strict Phase-B validators. The builder then copies only the exact
`instance.magma` and `manifest.json` bytes for each of the 160 instances, in
blind-bundle order. It does not copy `score.json`, an oracle ledger, prior
solver output, target classes, planted witnesses, selection priorities, or
known discrete-log labels. The packet seal inventories all 320 instance files,
the frozen protocol and the self-hashed manifest.

Prepare and verify the packet on the source host:

```sh
python3 scripts/run_koblitz_phase_b_magma_replay.py prepare \
  --run-root /Volumes/SSD990/koblitz-balanced-pdp-phase-b-run-successor-01-20260910 \
  --output /Volumes/SSD990/koblitz-phase-b-magma-replay-packet-successor-02-20260910

python3 scripts/run_koblitz_phase_b_magma_replay.py verify-packet \
  --packet /Volumes/SSD990/koblitz-phase-b-magma-replay-packet-successor-02-20260910
```

Only that sealed packet and a clean checkout containing the runner, process
meter and exact point validator go to the licensed host. The successor-04
bundle itself must not be transferred because it contains the post-run score.
The runner imports `koblitz_phase_b_magma_core.py`, a small execution-only
module with no Stage-13 archive paths or import-time historical checks.

## Licensed-host preflight and execution

The licensed-host operator supplies a UTF-8 statement describing authorized
license access. Preflight records that statement together with Magma,
MiniSat, point-validator and process-meter hashes; version output; OS, CPU and
Python host identity; the runner sources; and the exact repository state.
Production execution requires a recognized Magma V2 banner and clean runner
source state.
The production point validator must match the charged Phase-B backend exactly:
1,018,784 bytes with SHA-256
`14ec9a3d976371b6b100fa32a4ab916f50ddc88aad18fb6c6df8e81de925f96d`.

```sh
python3 scripts/run_koblitz_phase_b_magma_replay.py preflight \
  --packet /path/to/koblitz-phase-b-magma-replay-packet-successor-02-20260910 \
  --magma /licensed/path/magma \
  --minisat /path/to/minisat \
  --backend /path/to/koblitz_pdp_backend \
  --license-access-statement-file /path/to/license-access.txt

python3 scripts/run_koblitz_phase_b_magma_replay.py run \
  --packet /path/to/koblitz-phase-b-magma-replay-packet-successor-02-20260910 \
  --magma /licensed/path/magma \
  --minisat /path/to/minisat \
  --backend /path/to/koblitz_pdp_backend \
  --license-access-statement-file /path/to/license-access.txt \
  --output /path/to/new-magma-return
```

The output directory must be new. There is no resume or retry option. Tasks run
sequentially with one worker. Each Magma command uses `-t 1`; each copied input
contains `SetNthreads(1)`, `SetGPU(false)` and direct sparse F4. Each task has
at most one separately metered primary F4 process, one F4-plus-`SAT(G)` process
when F4 returns a proper basis, and one exact backend point-validation process
when `SAT(G)` returns a model. Every process has a 120-second watchdog.

Identity-bound execution scripts emit the task ordinal, split blind ID and
split source SHA-256 beside the ordinary F4 terminal markers. This detects
swapped output without changing the sealed source input. A primary F4 `UNSAT`
is a candidate terminal outcome for the balanced panel. A validated point
model is `SAT`. The first nonlifting model, timeout, missing terminal, process
error or resource failure is inconclusive because the protocol allows one
attempt and no retry. `Unknown` and timeout never imply `UNSAT`.

Magma does not expose SAT-style conflict counts on this path, so conflicts are
recorded as unavailable rather than zero. The process meter reports total CPU
seconds, wall time and process high-water RSS. Its historical
`single_core_seconds` field aliases total CPU; it is not an independently
measured single-core elapsed time. Magma's internal F4 timers are contained in
the process receipt and must not be added to it.

The runner writes `return-seal.json` last. It retains exact commands, attempt
starts, raw output, raw metrics, model assignments and validation reports. A
failed task is sealed as inconclusive and the runner continues to the next
untouched task; it never executes an already-started task again.

## Post-seal scoring

Return the sealed result directory to the source side. Verify it before opening
the score:

```sh
python3 scripts/run_koblitz_phase_b_magma_replay.py verify-return \
  --packet /Volumes/SSD990/koblitz-phase-b-magma-replay-packet-successor-02-20260910 \
  --return-root /path/to/returned-magma-results

python3 scripts/run_koblitz_phase_b_magma_replay.py score \
  --packet /Volumes/SSD990/koblitz-phase-b-magma-replay-packet-successor-02-20260910 \
  --return-root /path/to/returned-magma-results \
  --output /Volumes/SSD990/koblitz-phase-b-magma-replay-score-20260910
```

The scorer first verifies the write-once return seal and every retained process
receipt. Only then does it open successor-04 `score.json`, authenticated by the
original score seal, and classify the blind Magma outcomes. A local synthetic
run or a self-attested licensed-host return does not by itself establish an
unaffiliated reproduction. Novelty review, the remaining full-cost stages and
the Koblitz SOTA claim all stay false.

## Prepared successor packet

The hardened successor packet prepared after rebasing to the relocatable
Phase-B evidence bundle contains exactly 160 Magma inputs, 160 matching
manifests, the protocol, and the packet manifest. Its terminal identities are:

- packet-seal SHA-256: `1548c85abaecd61a4079235bb75705eaf6d085d5fb11660a575e04b620bb4342`;
- inventory SHA-256: `3a160c4bbfbcbcc9be36e20d739923e81212beb3ad3e2ad9cad22b94d3d4c856`;
- manifest SHA-256: `dcd2a602e6e9e2073483d6a9a1efd84bf5c5d983256b41ff20ce1f1d10e2d886`;
- protocol SHA-256: `3afd6baa75be22639b451a9fb4d021976baaf56e6c6d2edd5b35113656f1c6cd`.

The predecessor packet and successor-01 remain retained. They are superseded:
the hardened verifier rejects the predecessor preflight, and successor-01 is
bound to the older non-relocatable terminal-evidence bundle. No Magma task has
run on any packet.

The class-blind successor-02 packet is distributed with this branch as
`stage-22-phase-b-magma-replay-packet-successor-02-20260910.tar.gz`. Its archive
SHA-256 is
`97c87bce7da810674ead478f5a04ef6739c5da7e269287c42b21c91f1887fe54`.
CI checks the archive digest, rejects absolute paths, parent traversal, links,
and special files, extracts it into a fresh directory, and reruns
`verify-packet`. The archive contains no oracle, target class, known witness,
or previous solver output.
