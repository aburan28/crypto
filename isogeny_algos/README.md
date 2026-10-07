# isogeny_algos — isogeny-computation algorithms + benchmark harness

Standalone crate (no dependencies, own `[workspace]`). Every algorithm is a deliberately
**naive baseline** (u128 `%` modular multiplication, extended-Euclid inversion, schoolbook
polynomials and series, dense Gaussian elimination, affine curve arithmetic) so there is something
to iterate against. Nothing here is constant-time or intended for production.

**What is covered, what is not, and why: [`docs/SURVEY.md`](docs/SURVEY.md)** (status of every
algorithm in the literature I could identify, with the paper-level complexity where known).

## Implemented (V1 = first delivery, V2 = extension)

Different problems, different algorithms; they are benchmarked per problem, not against each other.

| Problem | Algorithms (module under `src/`) |
|---|---|
| kernel → isogeny | Vélu: odd (V1), **any subgroup incl. even degree / non-cyclic** (V2) `kernel/velu.rs` · Kohel formulas: odd (V1), **even / 2-torsion** (V2) `kernel/kohel.rs` · √élu (structure only) `kernel/sqrt_velu.rs` · **x-only Vélu from x(P)**, works for kernels over extension fields `kernel/xonly.rs` · **Montgomery x-only Vélu** `kernel/montgomery.rs` · **ℓⁿ / smooth kernels as chains, naive / balanced / cost-model strategies** `kernel/chain.rs` |
| (E, Ẽ) → isogeny | Padé on the ℘-series (V1) `find/elkies.rs` · **the BMSS family: linear algebra, Stark 1972, Atkin 1992, Atkin+modular composition, Elkies 1992, Elkies 1998, fastElkies, fastElkies′**, and σ from Φ's second derivatives (V2) `find/bmss.rs` |
| (E, ℓ) → ℓ-isogenies | Φ_ℓ roots + Elkies codomain (V1) · division-polynomial factoring (V1) `find/divpoly.rs` · Φ_ℓ mod p from q-expansions `find/modpoly.rs` · end-to-end with any BMSS method `find/bmss.rs::isogenies_via_phi` (V2) |
| (E₁, E₂) → isogeny | Galbraith BFS, GHS (+ volcano normalisation), Kohel volcano crater walk, Couveignes hard-homogeneous-space MITM (ordinary), Delfs–Galbraith (supersingular) (V1) · **Galbraith–Stolbunov weighted walk**, **CSIDH-style class-group action + MITM + ideal orders from cycles** (V2) `path/*` |
| auxiliary | **dual isogeny** `find/dual.rs`, **Kohel's End(E) conductor** `path/endo.rs` (V2) |

## Correctness checks (`cargo test --release`: 28 tests, all pass)

Each algorithm is checked against an independent computation, not only against itself:
Vélu = Kohel = √élu = x-only = Montgomery (via Weierstrass Vélu) on common kernels; E/E[2] has the j of E;
all eight BMSS methods reproduce Kohel's kernel polynomial and numerator for ℓ ≤ 101 (also Galois-stable
kernels); σ from Φ equals the true σ for ℓ ≤ 23; dual∘φ = [ℓ] and φ∘dual = [ℓ] on random points; SIDH-style
key exchange over F_{p²} agrees for every chain strategy; the CSIDH action commutes, is invertible, and every
step satisfies Φ_ℓ(j, j′) = 0; ideal-class orders divide an independently counted class number; the End(E)
conductor drops by exactly ℓ per ascent step; every path is verified edge by edge against Φ_ℓ. The benchmark
re-verifies each result before timing and stores `verified` in every record.

## Run

```
cargo test --release
cargo run --release --bin bench -- --out results/run.jsonl                 # V1 groups: kernel find path
cargo run --release --bin bench -- v2 --out results/run-v2.jsonl           # V2 groups: kernel2 chain bmss csidh v2path
cargo run --release --bin bench -- --quick kernel2                         # smoke test of one group
python3 scripts/report.py results/run-v2.jsonl > results/run-v2.md
```

## Baselines (single thread, 4-vCPU Xeon @ 2.1 GHz VM, shared host, no pinning: expect noise)

| file | content |
|---|---|
| `results/baseline.jsonl`, `baseline.md` | V1: kernel/find from the first run + regenerated path instances; `env.txt` |
| `results/baseline-v1-all.*` | the first complete V1 run, kept as evidence (its path instances were weaker) |
| `results/baseline-v2.jsonl`, `baseline-v2.md` | V2 groups; `env-v2.txt` |
| `results/baseline-v2-run1.jsonl` | the first, complete V2 run |
| `results/baseline-v2-csidh-rerun.jsonl` | `csidh` group re-run after a fix: the Montgomery ladder ran 128 iterations whatever the exponent size (per-step cost was flat ≈ 8 µs; now 2.0 µs at 20-bit p to 5.0 µs at 50-bit) |
| `results/baseline-v2-chain-rerun*.jsonl` | `chain` group re-run after two fixes: the strategy cost model charged n multiplications at every split instead of n−k, and the benchmark's calibration timed an evaluation at a kernel point (which returns at once). `…-rerun1-bad-calibration.jsonl` is the intermediate run with only the first fix |

`baseline-v2.jsonl` = `run1` for kernel2/bmss/walks/endo + the two re-runs. All 346 records verified.

Selected measurements (medians, this VM; no scaling or end-to-end claims are made):

* codomain from a generator, 32-bit p: general Vélu 287 ns (ℓ=3) to 144 µs (ℓ=1009); x-only Vélu 189 ns to 145 µs;
  Kohel with h given 1.4 µs to 4.06 ms.
* Montgomery x-only Vélu vs Weierstrass Vélu on the same CSIDH kernel: 1.3 µs vs 7.1 µs at ℓ=43.
* (E, Ẽ) → isogeny at ℓ=101, 40-bit p: Elkies 1998 74 µs, Elkies 1992 147 µs, Atkin 1992 2.18 ms, fastElkies 2.59 ms,
  Stark 73 ms, fastElkies′ 146 ms, linear algebra 1.65 s, V1 Padé 3.62 s (Kohel with the kernel given: 47 µs).
  The paper's O(M(ℓ)) bounds are not what schoolbook series give: fastElkies is slower than the O(ℓ²) Elkies 1998 here.
* ℓ² chain, e = 24, over F_{p²}: naive 263 µs (276 ℓ-multiplications), balanced 81 µs (60), calibrated cost model 96 µs
  (43 multiplications, 75 evaluations); for ℓ = 3 and e ≥ 10 the calibrated strategy is within 4 % of balanced (slightly faster); for ℓ = 2 it is 13–19 % slower.
* CSIDH-style action: 2.0 µs (20-bit p) to 5.0 µs (50-bit p) per isogeny step; MITM at n = 8 primes, m = 2:
  13 122 nodes, 123 ms; the order of each of 𝔩₃, 𝔩₅, 𝔩₇ at p = 1021019 is 1905 = h(−4p).
* identical 32-bit instances (5): GHS over all neighbours 17 ms (95 steps), Galbraith–Stolbunov uniform weights 2.4 ms (48),
  weights (16,8,4,2,1) 17 ms (732), Galbraith BFS 17.5 ms (93 nodes). Five instances: the spread is large.

## Known gaps and limitations (details and the full list in `docs/SURVEY.md`)

* **Not implemented** (implementation gaps): supersingular endomorphism-ring methods (quaternion orders, KLPT, ideal → isogeny,
  SQIsign), Mestre/Brandt graphs, higher-dimensional isogenies (Richelot, theta, Kani-lemma methods), characteristic-2/small-characteristic
  algorithms (Couveignes 1996 p-torsion, Lercier, Lercier–Sirvent), isogeny cycles for Atkin primes / full SEA, radical isogenies, other
  curve models (Edwards, Hessian, …), extension fields beyond F_{p²}. Not implementable here (resource limit): quantum algorithms.
* **√élu has no asymptotic saving yet**: the I/J/K structure is implemented but the product is evaluated by naive Horner (O(ℓ)).
  Likewise the BMSS methods use schoolbook arithmetic, so the paper's complexity bounds are targets, not achieved.
* Odd ℓ only for √élu and the BMSS family; Elkies/BMSS reject j ∈ {0, 1728}; `explicit_chain` fails through such j.
* Φ_ℓ from q-expansions is practical to ℓ ≈ 23–31 only (dense linear algebra).
* Couveignes (ordinary) is the class-group-action algorithm with planted exponents, not the 1996 p-torsion method. Kohel volcano path
  only connects curves in one ℓ-volcano. GHS with volcano normalisation needs the trace and failed on 3 of 15 conductor > 1 instances (kept in the data).
  Galbraith–Stolbunov is the weighted-prime walk only.
* Delfs–Galbraith returns a j-invariant path (no explicit F_{p²} isogeny), needs p ≡ 7 mod 8 here, and 2 of 5 instances at 38 bits found no path.
* End(E) conductor: ℓ = 2 dividing the conductor of Z[π] is unsupported; every ℓ | f_π needs Φ_ℓ.
* The CSIDH benchmark has no 14-prime set (no prime p = 4∏ℓ − 1 < 2⁶² among the candidates tried); sizes are 20–50 bits.
* Path-finding instances are small (paths ≤ 35 at 20–32 bits, primes {3,5,7}); class-group sizes were not measured, so do not read
  scaling exponents from these runs.
* No Conductor task could be obtained for this work (the control plane needs Docker/Postgres, unavailable here); all edits are confined to `isogeny_algos/`.
