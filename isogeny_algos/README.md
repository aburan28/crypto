# isogeny_algos — baseline isogeny-computation algorithms + benchmark harness

Standalone crate (no dependencies, own `[workspace]`). Every algorithm is a deliberately
**naive baseline** (u128 `%` modular multiplication, extended-Euclid inversion, schoolbook
polynomials, dense Gaussian elimination, affine curve arithmetic) so there is something to
iterate against. Nothing here is constant-time or intended for production.

## What is implemented

Three different problems are solved by different algorithms; they are benchmarked per problem, not against each other.

| # | Algorithm | Reference | Problem | Module |
|---|---|---|---|---|
| 1 | Vélu | Vélu 1971 | kernel (generator) → codomain + map | `kernel::velu` |
| 2 | Kohel formulas | Kohel thesis 1996 | kernel polynomial → codomain + rational x-map | `kernel::kohel` |
| 3 | √élu | Bernstein–De Feo–Leroux–Smith 2020 (adapted to short Weierstrass) | generator → codomain + x-map | `kernel::sqrt_velu` |
| 4 | Division-polynomial factoring | Schoof / Elkies / Couveignes ℓ-torsion style | E, ℓ → all rational ℓ-kernels | `find::divpoly` |
| 5 | Φ_ℓ roots + Elkies codomain + BMSS | Elkies 1998, Bostan–Morain–Salvy–Schost 2008 | E, ℓ → all rational ℓ-isogenies | `find::{modpoly, elkies}` |
| 6 | Bidirectional BFS | Galbraith 1999 | E1, E2 → path in ℓ-graph | `path::galbraith` |
| 7 | Random-walk collision (+ Kohel volcano normalisation) | Galbraith–Hess–Smart 2002 | E1, E2 → path | `path::ghs` |
| 8 | Volcano navigation + crater walk | Kohel 1996 | E1, E2 in one ℓ-volcano → path | `path::volcano` |
| 9 | Class-group-action meet-in-the-middle | Couveignes 2006 / Rostovtsev–Stolbunov | E2 = [a]E1 → exponent vector + path | `path::couveignes` |
| 10 | Delfs–Galbraith | Delfs–Galbraith 2016 | supersingular E1, E2 over F_p² → j-path | `path::delfs_galbraith` |

Φ_ℓ is computed mod p from q-expansions (`find::modpoly`); checked against the known Φ₂.
`path::graph::explicit_chain` converts a j-path into explicit normalised isogenies (Elkies + BMSS per step).

## Correctness checks (all pass: `cargo test --release`)

- Vélu, Kohel and √élu agree on codomain and x-map for ℓ ∈ {3,5,7,11,13,29,61,101}; Vélu is checked to be a group homomorphism on random F_p-points.
- Division-polynomial kernels ⊇ the known kernel and each passes the identity (x³+ax+b)·f′² = f³+a′f+b′; Elkies/BMSS reproduces Kohel's numerator, codomain and kernel exactly for ℓ ≤ 13.
- Every path is checked edge-by-edge (Φ_ℓ(j_i, j_{i+1}) = 0); explicit chains are applied to random points and checked to land on the final curve.
- The benchmark re-verifies each result before timing and stores `verified` in every record.

## Run

```
cargo test --release
cargo run --release --bin bench -- --out results/run.jsonl     # all groups; ~30 min
cargo run --release --bin bench -- --quick kernel              # smoke test
python3 scripts/report.py results/run.jsonl > results/run.md
```

## Baseline

`results/baseline.jsonl` (raw), `results/baseline.md` (tables), `results/env.txt` (machine:
4-vCPU Xeon @ 2.1 GHz VM, single thread, no pinning, shared host, so expect noise).
`results/baseline-v1-all.*` is the first complete run, kept as evidence; its path-finding
instances were weaker (some start curves had one splitting prime, giving tiny components) and
were regenerated in `baseline-path-v2.jsonl`; `baseline.jsonl` = v1 kernel/find + v2 path.

Measured (medians, see tables for everything): at 32-bit p the codomain from a generator takes
Vélu 340 ns (ℓ=3) → 142.6 µs (ℓ=1009); Kohel with h given 491 ns → 3.87 ms; √élu 441 ns → 80.3 µs.
Finding all kernels at ℓ=23 (61-bit p): division-poly factoring 2.10 s, Φ-roots+Elkies+BMSS 23.0 ms
(plus a one-off 150 ms to build Φ₂₃).

## Known gaps and limitations (not yet implemented / not established)

- **√élu has no asymptotic saving yet**: the I/J/K structure and series-ring products are implemented,
  but the resultant is evaluated by naive Horner over the I points (O(ℓ), not Õ(√ℓ)); no multipoint evaluation.
- **BMSS Padé step uses dense ℓ×ℓ Gaussian elimination** (O(ℓ³)), not the Newton/Euclid versions.
- Kernel algorithms handle **odd ℓ only**; Elkies/BMSS rejects j ∈ {0, 1728} and degenerate Φ derivatives;
  `explicit_chain` fails for paths through such j.
- Φ_ℓ construction is practical to ℓ ≈ 23–31 only (dense linear algebra); no Φ_ℓ for larger ℓ.
- Couveignes is the **class-group-action (HHS) algorithm**, with exponents planted in a box; it is not
  the 1996 "ℓ-isogenies from the p-torsion" method for small characteristic, and is not directly
  comparable with generic path finding. Prime selection requires ℓ ∤ conductor and distinct eigenvalue classes.
- Kohel volcano path only connects curves in the **same ℓ-volcano**. GHS with volcano normalisation
  needs the Frobenius trace; it failed on 3 of 15 conductor>1 instances where the remaining split primes did not connect the two curves.
- Delfs–Galbraith returns a **j-invariant path only** (no explicit F_p² isogeny), requires p ≡ 7 mod 8 here
  (for p ≡ 3 mod 8 the F_p-rational 2-graph is a forest of stars and BFS cannot connect), and 2 of 5
  instances at 38 bits found no path within the F_p subgraph search (raw records kept).
- Path-finding instances are small: path lengths ≤ 35 at 20–32 bits (E2 = random walk from E1 over primes
  {3,5,7}); class-group sizes were not measured, so **do not read scaling exponents from these runs**.
- Benchmark workloads: kernel group uses curves with a rational ℓ-torsion point (BSGS point counting);
  61-bit runs stop at ℓ=401 because generating ℓ=1009 curves is too slow with this arithmetic.
- No Conductor task could be obtained for this work (see PR/report notes); all edits are confined to `isogeny_algos/`.
