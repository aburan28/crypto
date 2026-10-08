# Inherited-F4 overhead round: result

**Outcome: success under the preregistered rule. Class: engineering.**
On 48 frozen inputs every output is identical. The candidate executes
**0.847×** the baseline's instructions on the confirmatory m = 19 inputs
(threshold 0.90), and no input regresses. The ecbench `S` is unchanged by
construction: no row, order, layout or charge moved, so neither did any
counted unit. What fell is the machine work the counted unit does not
see. At m = 19 the inherited-F4 arm went from about 120 to about 95
instructions per counted word XOR.

Protocol: [`PROTOCOL.md`](PROTOCOL.md), committed with the candidate code
before any confirmatory run (`4d9355cb`). Raw results: `results/run-1/`.

## Host and binaries

From `results/run-1/host.txt` and `binaries.sha256`:

- Intel Xeon @ 2.80 GHz (cloud VM), 4 logical CPUs, 16 GiB, Linux 6.18
  x86-64, with `popcnt avx2 avx512f bmi2 pclmulqdq`.
- `rustc` 1.97.0, Valgrind 3.22.0.
- **Baseline**: `main` at `fa80835a`, `ecbench` SHA-256 `48d34563…`.
  Its A/A copy is byte-identical.
- **Candidate**: `4d9355cb`, `ecbench` SHA-256 `6eddf7b4…`.

Both are release builds of the default profile. The binaries are not
committed: they rebuild from those revisions.

## The gate: identical outputs

`results/run-1/gate.tsv` has 48 inputs and 48 `IDENTICAL`.

| set | inputs | what they are |
| --- | --: | --- |
| L13 | 16 | ladder m = 13, round 1, both arms |
| L19 | 16 | ladder m = 19, round 1, both arms |
| L23 | 4 | ladder m = 23, two workloads, both arms |
| H13 | 8 | fresh holdout targets, m = 13 (seed 202610052113) |
| H19 | 4 | fresh holdout targets, m = 19 (seed 202610052119) |

"Identical" means `jq -S` byte equality after removing only timing and
host-placement fields (`strip.jq`). That covers the recovered scalar, the
verification, every phase's counted operations, `solver_ops`, every solver
statistic and the total GAE. The table records a SHA-256 of each stripped
output.

## The table: instructions (Callgrind `Ir`, whole process)

One run per (binary, input); `Ir` is deterministic. Ratio is candidate /
baseline. In the ecbench unit, `S / floor` and `S / rho` are unchanged
(identical counts) and are those of the
[F6-IC ladder](../f6_ic_ecbench_ladder_20261005/RESULT.md).

| m | arm | set | inputs | baseline `Ir` | candidate `Ir` | ratio (pooled) | per-input range | `S` | class | correct |
| --: | --- | --- | --: | --: | --: | --: | --: | --- | --- | --- |
| 13 | inherited F4 | ladder | 8 | 1.036 × 10¹⁰ | 8.932 × 10⁹ | **0.863** | 0.856–0.878 | unchanged | engineering | 8/8 identical |
| 13 | inherited F4 | holdout | 4 | 5.251 × 10⁹ | 4.533 × 10⁹ | **0.863** | 0.858–0.871 | unchanged | engineering | 4/4 identical |
| 13 | F6-IC | ladder | 8 | 5.868 × 10⁹ | 5.558 × 10⁹ | **0.947** | 0.945–0.954 | unchanged | engineering | 8/8 identical |
| 13 | F6-IC | holdout | 4 | 2.995 × 10⁹ | 2.832 × 10⁹ | **0.946** | 0.943–0.948 | unchanged | engineering | 4/4 identical |
| 19 | inherited F4 | ladder | 2 | 1.225 × 10¹¹ | 9.722 × 10¹⁰ | **0.793** | 0.792–0.795 | unchanged | engineering | 2/2 identical |
| 19 | inherited F4 | holdout | 2 | 1.388 × 10¹¹ | 1.101 × 10¹¹ | **0.793** | 0.793–0.794 | unchanged | engineering | 2/2 identical |
| 19 | F6-IC | ladder | 2 | 6.977 × 10¹⁰ | 6.395 × 10¹⁰ | **0.917** | 0.916–0.917 | unchanged | engineering | 2/2 identical |
| 19 | F6-IC | holdout | 2 | 7.924 × 10¹⁰ | 7.260 × 10¹⁰ | **0.916** | 0.916–0.916 | unchanged | engineering | 2/2 identical |

**Decision statistic.** The confirmatory m = 19 inputs are L19 seq 24, 27
and 28 and H19 0–3; L19 seq 25 was used in exploration. Their pooled
ratio is 291,811,013,044 / 344,582,852,768 = **0.8469**, which is ≤ 0.90.
The worst single input is L13-39 (F6-IC) at 0.954, so no input regresses.
Pooled over all 32 `Ir` inputs the ratio is 0.841.

**Instructions per counted word XOR.** At m = 19 the inherited-F4 arm drops
from 119.3–119.8 to 94.8–95.1. The F6-IC arm drops from 91.9–92.4 to
84.2–84.6. The gain grows with m (inherited F4: 0.863 at m = 13, 0.793 at
m = 19), because the completion-row bookkeeping it removes grows with the
layout. F6-IC gains less because it does fewer reductions, and the
geometric work in its gate is untouched.

## Time (practicality note)

`results/run-1/wall.tsv`: process CPU time from the child output, pinned
with `taskset -c 3`, run alone on the host. There were five rounds, and
each ran baseline, candidate, the byte-identical baseline copy (A/A), then
the candidate again. The A/B ratio is the round's mean candidate over its
baseline. Scheduling delay was at most 9.9 ms per run (mean 2.2 ms), so no
run is marked contended. The VM's neighbours and clock are outside the
pinning (AGENTS.md §10).

| input | arm | baseline CPU, median (min) | candidate CPU, median (min) | A/B, median (range) | A/A, median (range) |
| --- | --- | --: | --: | --: | --: |
| L19 seq 25 | inherited F4 | 8.58 s (8.29) | 7.12 s (6.97) | **0.840** (0.810–0.846) | 1.011 (0.989–1.061) |
| L19 seq 24 | F6-IC | 5.03 s (4.93) | 4.77 s (4.64) | **0.937** (0.897–0.960) | 0.991 (0.953–1.011) |
| L23 seq 25 | inherited F4 | 110.4 s (108.9) | 91.7 s (90.3) | **0.834** (0.822–0.845) | 1.003 (1.000–1.023) |

Every inherited-F4 A/B round falls outside its A/A spread. The F6-IC rounds
overlap the A/A spread at its low end (0.953), so that arm's time gain is
not resolved from noise on this host; its instruction gain (0.917) is. At
m = 23, where `Ir` was not measured, the time ratio (0.834) tracks the
m = 19 ratios. Wall time is a practicality note here and decides nothing.

## What changed, and where the instructions went

The exploratory profile of L19 seq 25 (`PROTOCOL.md`, exploratory record)
showed that word XORs are a small part of the process's instructions:

| baseline share of `Ir`, L19 seq 25 | what it is |
| --: | --- |
| 19% | `rewrite`: the per-bit column remap that materialises a row in a new layout (counted as specialisation, per word) |
| 18% | `insert`: re-reducing a displaced row, which includes the counted XORs |
| ~13% | sorting, almost all of it completion rows and layout extension |
| ~8% | `malloc` / `free` |
| 8% | `specialise_shared`'s own bookkeeping |
| 5% | substituting the system at every node |

The candidate removes most of the sorting and three of the four hash
lookups per completion monomial, rebuilds no column index after a layout
grows, inlines the pivot is-current test, and vectorises the pivot XOR.
Everything left is in `rewrite`, `insert`, allocation and substitution.

Rejected, with its number: a sink-word scatter in `rewrite`, which replaced
the deleted-column branch with a clamped store. It cost 2.5 × 10⁹ more
instructions on L19 seq 25, because the branch skips the deleted bits
cheaply.

## What this does not establish

- **No `S`, ratio or crossover moves.** The ecbench `S` is a lower bound
  in counted units and is identical here. The ladder's verdict stands: IC
  is slower than rho at m = 13, 19 and 23. A cheaper instruction stream
  lowers only the *measured* solver price (`ns_per_op`) and any
  calibration ratio re-measured from it. Re-pinning
  `docs/ic/calibration.json` would be a separate accounting round.
- **One host class.** All of this ran on x86-64 Linux in a cloud VM. The
  instruction counts are host-independent for this binary, but the binary
  is compiled for the baseline x86-64 target, so other targets (Arm64,
  other x86 builds) are not measured. GPUs are not relevant to this code.
- **No asymptotic claim.** The engine's exponent is set by how many
  candidate tuples must be ruled out, not by its per-node constant.
