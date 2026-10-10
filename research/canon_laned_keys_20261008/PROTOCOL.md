# Laned Frobenius-orbit keys: protocol

Declared before the evidence runs. The change is commit `a98ba95a`
(`src/cryptanalysis/koblitz_fast.rs`). It is an engineering change to how a
folded pair table keys a point. No key value changes.

## Question

The m = 3 pair-table scan and the folded pair-table build key every point by
the least rotation of its abscissa's normal-basis coordinates
(`FrobeniusCanon`). Without AVX-512F, `canon_in_place` computes those keys in
one of two ways:

- one scalar search per key;
- in `point_keys`, every one of the `n − 1` rotations of eight words.

The candidate computes the same keys with `least_rotations::<8>`: doubling,
binary lifting and selects, with the lanes in lockstep. Does that lower the
cold cost of the pipeline on that path, with identical results?

## Arms

- **Baseline:** `ic` built from `7424539b`, the parent of the change, unmodified.
- **Candidate:** `ic` built from `a98ba95a`.

Both are release builds from the same toolchain on the same host. The host
manifest is in `runs/host.json`.

## Inputs (frozen)

The six §23 parameter files, one public target each:

| file | sha256 (first 16) |
|:--|:--|
| `research/ic_single_target_20260930/runs/k0n41/T01.params.json` | `8ff9b683716a9359` |
| `research/ic_single_target_20260930/runs/k1n47/T01.params.json` | `02a4902b29bbcb90` |
| `research/ic_single_target_20260930/runs/k0n53/T01.params.json` | `a65d549a056bb626` |
| `research/ic_single_target_20260930/runs/k0n57/T01.params.json` | `a7d19389a74809fa` |
| `research/ic_single_target_20260930/runs/k1n59/T01.params.json` | `d32fecc809cf4db5` |
| `research/ic_single_target_20260930/runs/k0n61/T01.params.json` | `34a69de14f899439` |

## Paths

- **portable** (`KIC_SCAN_SIMD=0`): the path on hosts without AVX-512F, and
  the path Callgrind measures, since Valgrind hides AVX-512.
- **avx512** (default on this host): the AVX-512 kernel, which the change does
  not touch. This is a control, expected unchanged.

## Measurements

1. **Native phases.** Each run is
   `ic price --params P --repeats 1 --repeats-fast 1 --json`, one target and
   no rho arm. It runs under `isolated_bench run --cpus 3` with
   `RAYON_NUM_THREADS=1`.
   - Recorded per run: `collect`, `build`, `total_ns`, `total_units`, counts,
     recovered logarithms and verification.
   - A/A first: the baseline against itself, 3 rounds per size and path.
   - Then A/B interleaved, baseline then candidate, 5 rounds per size and path.
2. **Instruction counts.** Callgrind `Ir` for the whole process, baseline and
   candidate, at n = 53 and n = 61. This is deterministic.
3. **Rho reference.** `ic price --single-target --rho-seed 7` at n = 61, 3
   interleaved rounds per arm. It checks that rho's online time and step count
   are unchanged. The change leaves the scalar search that rho calls untouched.

## Identity gate

Per size and path, every run of both arms must give identical `counts`,
`recovered` and `all_verified`. Any difference stops the comparison: the
candidate is then a different algorithm, not a speedup.

## Success and stop conditions

- **Portable, success:** at every size, the median candidate `collect` and
  `total_ns` are below the baseline's by more than the A/A spread (max/min − 1
  of the baseline's own rounds).
- **avx512:** the candidate stays within the A/A spread. A difference beyond it
  is reported as a regression or as unexplained.
- **Rho:** step counts must be identical. Online time must stay within its A/A
  spread.
- **Stop:** an identity failure, a failed verification, or a contended run that
  repeats twice (the contended attempt is kept).

## Not claimed

- **No AVX-512 gain.** No gain on AVX-512 hosts, the class this host and §23's
  belong to.
- **No Arm64 measurement.** Apple silicon, Graviton and AMD/consumer x86 hosts
  run the portable path, but none is measured here. Their gain is unmeasured.
- **No ECC2K-130 transfer.** `FrobeniusCanon` exists only for n ≤ 63
  (single-word fields). The m = 83 gate (§8a) and m = 131 use wide fields and
  are untouched.
- **No change to the dashboards.** No IC/rho ratio in Table A, the scoreboard
  or the progress chart changes, because those are §23 measurements on the
  avx512 host class.
- **Not run: the frozen WDSat regression suite (§8).** It measures a SAT solver
  stage this change does not touch. The matched suite here is the six frozen
  §23 parameter files.
- **Class:** engineering (§3). Constants on one key path, no algorithmic change.
