# LLL degeneracy probes, re-run 2026-10-16

This report records observations only: it re-runs the committed probes behind
`papers/lll-degeneracy-note` to correct that note. A sheriff audit
(crypto-autoresearcher `research/sheriff/SHF-20261016-3caa6f`, blocker P1)
found that the note reports, as curve properties, LLL failures that
`research/notes/lattice-hnp/RESEARCH_LLL_GS_ANALYSIS.md` §10 had already
traced to floating-point Gram–Schmidt bugs and fixed in May 2026.

## What ran

| Item | Value |
| --- | --- |
| Code | commit `97d89f0fd52a`. This is `origin/main` 3e7bcb7 plus a `src/zk`-only change; `src/cryptanalysis/lattice.rs` and `hnp_ecdsa.rs` are untouched. |
| Probes | `cargo test --release --test lll_degeneracy_probe -- --include-ignored --nocapture --test-threads=1 --skip p521` (`logs/probes.log`) |
| P-521 probe | `cargo test --release --test lll_degeneracy_probe probe_p521_lll_hp -- --exact --ignored --nocapture` (`logs/p521_hp.log`) |
| Parameters | `m = 8` signatures, `k_bits = 3/4 · n_bits`, `δ = 0.75`, seeds as in the logs |
| Host | `logs/host.txt`: 4 vCPU cloud VM, not isolated |

Timings are wall time on that host, one run each. They are not a benchmark.

## Results

These come from `results.json`, which `parse.py` builds from the logs, and are
drawn in `recovery.svg` by `make_figure.py`.

| Curve | Gram–Schmidt path | Recovered / probes | Note's pre-fix report |
| --- | --- | --- | --- |
| secp192k1 | scaled f64 | 3 / 3 | recovered |
| secp224k1 | scaled f64 | 3 / 3 | recovered |
| secp256k1 | scaled f64 | **6 / 6** | **0 / 3**, "curve-specific" |
| P-256 | scaled f64 | 6 / 6 | 3 / 3 |
| brainpoolP256r1 | scaled f64 | 1 / 1 | recovered |
| P-384 | scaled f64 | **4 / 4** | **iteration cap** |
| brainpoolP384r1 | scaled f64 | **4 / 4** | **iteration cap** |
| P-521 | scaled f64 | 0 / 4 | iteration cap |
| P-521 | 2048-bit fixed point (`lll_reduce_hp`) | **3 / 3** | not run |

Each recovered probe returned the planted key `d`. The probe compares the
recovered key with `d` (`tests/lll_degeneracy_probe.rs`, `probe_once`). Every
probe at 256 or 384 bits finished in 13–27 ms. The note reported about 10 s and
about 20 s "iteration cap" failures.

![recovery](recovery.svg)

## Reading

- **Failure mode A (secp256k1 at 256 bits).** It does not reproduce. As
  `RESEARCH_LLL_GS_ANALYSIS.md` §10.1 records, the cause was an `f64` overflow
  in the Gram–Schmidt squared norms: `n⁴` exceeds `f64::MAX` for secp256k1's
  order but not for P-256's or brainpoolP256r1's. The sheriff case checks this
  with exact arithmetic (`checks/f64-overflow.txt` there). The scaling fix in
  `src/cryptanalysis/lattice.rs` removes it.
- **Failure mode B (384 bits and above).** It does not reproduce at 384 bits:
  the same scaling fix suffices, so the iteration cap was not the cause. P-521
  still fails on the scaled `f64` path, through the catastrophic cancellation
  of §10.2, and recovers 3 of 3 with the 2048-bit path of §10.3.
- **The `probe_koblitz_lll_degeneracy_hypothesis` closing line.** It prints
  "Hypothesis confirmed", but its near-power-of-two bucket is empty (0/0), so
  the line is vacuous. `RESEARCH_GLV_HNP.md` already says so. This re-run
  inherits that message, and it supports nothing.

## Limits

- **Small samples.** One to six probes per curve, at `m = 8` only, on one
  unisolated host.
- **Scope.** The re-run establishes that the note's failures do not reproduce
  with the fixed code. It establishes nothing about LLL cost across curves, and
  it is not a security statement about any curve.
