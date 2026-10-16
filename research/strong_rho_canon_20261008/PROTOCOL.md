# Strong-rho canonicalization and batched-inversion scratch: protocol

Declared before the evidence runs. Two engineering changes are measured
together against their parent:

| commit | change | code it touches |
|:--|:--|:--|
| `23ed8d3a` | `StrongRhoG::canonicalize` finds the least rotation from the longest zero runs instead of trying all `n − 1` rotations | the strong rho reference, every width (`StrongRho`, `WideStrongRho`) |
| `a8f4bb60` | `Gf2::batch_inv` and the AVX-512 `add_many_lazy` stop zero-filling their scratch on every call | every one-word batched inversion: the IC scan and build, and the narrow rho walks |

No output may change: same canonical forms, walks, relations and logarithms.

## Arms

- **Baseline:** the library at `7ce52fa2`, the parent of both changes. The
  measurement example `examples/strong_rho_step_rate.rs` is added to that tree
  unchanged. It uses only public API that both arms share.
- **Candidate:** `a8f4bb60`.

Both arms are release builds in one target dir, on the same host and
toolchain. The binary hashes are in `runs/host.json`.

## Measurements

Every timed process runs under `isolated_bench run --cpus 3 --max-other-cpu
0.3` with `RAYON_NUM_THREADS=1`. The 0.3 allowance is because this container's
session harness uses about 0.1 CPU on another core. Each measurement has an
A/A pass (baseline twice, 3 rounds), then 5 interleaved A/B rounds.

1. **M1: m = 83 reference** (`icv1-f2m83-tm6151469093347-debefd74`).
   - Command: `strong_rho_step_rate 262144`.
   - Recorded: canonicalization ns per state over a fixed sequence of 262,144
     states, the digest of every canonical form, and the wall time of a walk
     capped at 10^6 steps (default 32 lanes, 4 distinguished-point bits).
   - This isolates `23ed8d3a`: the wide field's batched inversion is not
     touched by `a8f4bb60`.
2. **M2: narrow strong rho.**
   - Command: `koblitz_rho_fixture n a signed_frobenius 4 strong` on the six
     §23 curves: n = 41, 53, 57, 61 with a = 0, and n = 47, 59 with a = 1.
     That is 4 published synthetic fixtures per process, each solved to the end.
   - Recorded: the process wall time from the isolation record, and every
     output line.
   - Both changes act here.
3. **M3: index-calculus pipeline.**
   - Command: `ic price --params P --repeats 1 --repeats-fast 1 --json` on the
     frozen §23 files for n53 and n61, portable (`KIC_SCAN_SIMD=0`) and AVX-512.
   - Only `a8f4bb60` acts here.
   - Also: Callgrind instruction counts for the whole process, portable n61,
     one run per arm.

## Identity gate

Any difference stops the comparison:

- **M1:** equal digests and both walks reaching the cap.
- **M2:** identical output lines, timing aside. That means charges, restarts,
  fruitless cycles and the recovered scalar.
- **M3:** identical `counts`, `recovered` and `all_verified`.

## Success and stop conditions

- **M1:** the candidate's capped walk is faster than the baseline's by more
  than the A/A spread.
- **M2, M3:** report the ratio against the A/A spread. A gain inside the spread
  is not established. A loss beyond it is a regression and is reported as one.
- **Stop:** an identity failure, or a contended run that repeats.

## What this changes, and does not

- **A faster reference.**
  - The strong rho is the `vs_rho` reference of `docs/ic/boundary_targets.json`
    and of ecbench's strong methods. Making it faster makes every future
    index-calculus-against-rho ratio measured with it harder for index
    calculus.
  - Committed sessions keep their numbers. Where a dashboard cites a strong-rho
    figure measured before this change, that figure prices rho with the slower
    canonicalization.
  - This round reports by how much (M1, M2). It does not re-measure any
    dashboard row.
- **No index calculus at m = 83.** There is no pair-table index calculus at
  m = 83 (`FrobeniusCanon` stops at n = 63), so M1 measures the reference only.
- **No other hosts.** No Arm64, AMD or GPU host is measured. The GPU ECC2K-130
  rho is separate code and is untouched.
- **Class:** engineering (§3). Constants only, no algorithmic change.
