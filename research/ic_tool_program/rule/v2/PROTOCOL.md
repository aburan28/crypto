# The rule's comparison at baseline v2

Declared 2026-10-02, before any run below. What ran first, in scratch
space, with no number quoted or kept in the record:

- `icprog rule claims` and `icprog rule analyse` (N3) reproduced ledger
  §23's claims, manifests and `analysis.json` from §23's run tree, byte
  for byte (`tests/icprog.rs`).
- A smoke run of the pin below: v2 priced `T01` at each of the six sizes
  untimed, and decided each as §23's binary did.
- The strong rho fixture walked `T01` at each size untimed, to see that
  it runs at all six and hashes §23's targets.

## Question

IC_TOOL_PROGRAM.md asks, at every new baseline, for the rule's
comparison: one unseen public point, the index calculus against Pollard
rho on that same point, online and cold. It runs at §23's six sizes with
64 targets, as §23 did, and it moves the scoreboard's §23 row.

v2 is the baseline R05 made (#1187): `edcb0bec`, binary
`76a2a2fd05acd47d6dc9a0851a013d8056cc9d7d9b7a34320d1d92089a794597`.
v1 has no comparison of its own; v2's supersedes it.

## What is §23's, unchanged

Everything §23's protocol
([`research/ic_single_target_20260930/PROTOCOL.md`](../../../ic_single_target_20260930/PROTOCOL.md))
fixes, except the rows under "What changes":

- **The sizes**, in the order of `r`:
  - `icv1-f2m47-t22705043-f4e44623`;
  - `icv1-f2m57-tm747311035-c1f545af`;
  - `icv1-f2m41-tm2308219-7f48b14a`;
  - `icv1-f2m53-tm56619371-dac20a85`;
  - `icv1-f2m59-tm943548413-98844ecc`;
  - `icv1-f2m61-t158598901-ab42b6c5`.
- **The rows.** At each size, targets `T01`–`T64` once (`R1`), then
  `T01`–`T04` again (`R2`, the A/A).
  - Each row prices §23's own parameter file, copied from §23's run tree
    and checked byte for byte on every use.
  - Its rho seed is `0x230000 + i`.
- **The arms.** One `ic price --single-target` process a row: three
  repetitions, the index calculus then rho, both rebuilt each time. The
  online interval, its five exclusive phases and the reusable set-up are
  §23's.
- **The retries.** A contended or failed process runs again, at most
  twice. The first clean, complete attempt is the row, and every attempt
  stays.
- **The stop rule.** A row that does not complete stops its size.
- **The claims.** For each row:
  - the IC1 candidate, the workload and the run ID;
  - both arms' certificates replayed outside `ic`;
  - rho's EC1 reference identity;
  - the `vs_rho` claim, checked by the repository's checker.
- **The reading.** The same figures, bootstrap intervals (10,000
  resamples, §23's seeds), fits and diagnostics.
- **The inadmissible list.**

## What changes

- **The binary** is v2.
- **The tools are native** (AGENTS.md, no Python):
  - `icprog rule` (N3) runs the rows, writes the claims and reads the
    figures, in place of `run.py`, `claims.py` and `analyse.py`;
  - `isolated_bench` (N1) isolates each process, in place of
    `tools/isolated_bench.py`.
- **The identities differ from §23's, as they should.**
  - The resource envelope names the native isolation tool, so every
    workload ID is new.
  - `koblitz_index_calculus.rs` changed in v1 and v2, so the candidate
    IDs are new too.
- **The replay's checker** is `icprog`'s oracle: field arithmetic of its
  own, in a process apart from `ic`.
- **The pin.** It runs before any size.
  - `T01` at each size is priced untimed under `taskset -c 2`.
  - Each must give the outputs a pin compares equal to §23's `T01-R1`:
    the status, the counts, both arms' logarithms, rho's counts and both
    checks.
  - Every name must be the curve's slug.
  - If the pin fails, the comparison stops.
- **Each row's logarithms** are checked against §23's row: the same
  target, and the same logarithm from both arms.
- **What each size is read against.** §23 set its figures beside a
  prediction. This comparison sets them beside §23's measurement of the
  same figure (`s23_then_v2`).
- **The host** is this container, kernel build `6.18.44-fc-v51`. Its
  manifest is `runs/host.json`.

## The reference check (`docs/ic/BOUNDARY_TARGETS.md`, 2026-10-01)

Since §23, the ledger requires every `vs_rho` comparison to measure
against a strong rho: `koblitz_rho_fixture <n> <a> signed_frobenius 1
strong`. A different rho is admissible only with a measured per-step
cost, on the same target, no worse than the strong one's. `ic price`'s
rho is a different walk, so this comparison measures that.

- **The strong walk.**
  - It is `koblitz_rho_fixture`'s `strong` backend with its defaults: 32
    lanes, 4 distinguished-point bits.
  - It walks targets `T01`–`T08` at every size: `hash:<23000 + i>`,
    which is `ic workflow`'s domain, so the point is the row's own.
  - Its walk seed is `0x230000 + i`.
  - Each walk is one process through `isolated_bench`, with the same
    retries.
- **The fixture's rungs** now include 47, 57 and 61, the three sizes of
  §23 they missed. Nothing else in the fixture changes.
- **Per step.**
  - The strong walk's cost is `walk_ms / walk_steps`.
  - `ic price`'s rho cost is the row's median online wall over its
    steps.
  - Each is the median over the eight targets.
- **Admissible** at a size when all three hold:
  - `ic price`'s rho costs no more per step than the strong walk;
  - every strong walk verified its logarithm;
  - every strong target is the row's.

  At a size where it is not admissible, that size's online and cold
  ratios against rho are reported, but marked inadmissible under the
  ledger's minimum. They do not move the scoreboard.

## Accounting, classes and scope

- **The headline figures** are §23's:
  - the online speedup;
  - the cold ratio, index calculus over rho;
  - the break-even count.
- **Class.** This comparison changes no algorithm. It re-measures the
  method at a new baseline, so it is not a round and has no class of its
  own: v2's class is R05's (engineering). It reports where v2 stands
  against rho.
- **What moves.** The scoreboard's §23 row moves to v2's figures,
  keeping §23's as the before mark. So do the ledger's online/rho and
  cold/rho columns for v2.
- **Non-claims** are §23's:
  - one public target per row;
  - set-up reported apart from the online interval;
  - no precomputed rho table;
  - nothing past `n = 61`, and `m = 83` is not run (AGENTS.md §8a);
  - one x86-64 container.

## Cost, as an estimate

§23 took about 1.5 hours on this host. v2's set-up is 20–30% cheaper at
the top sizes, so about 1.2–1.5 hours. The reference check's 48 walks
take minutes.

## Commands

    icprog rule all --comparison v2 --ic <v2 ic> --ic-commit edcb0bec948ffab566a3c8998dc3c92d52375da2 --isolate <isolated_bench>
    icprog rule reference --comparison v2 --fixture <koblitz_rho_fixture> --fixture-commit <commit> --isolate <isolated_bench>
    icprog rule claims --comparison v2
    icprog rule analyse --comparison v2 > research/ic_tool_program/rule/v2/analysis.json
