# isogeny_algos — isogeny-computation algorithms + benchmark harness

Standalone crate (no dependencies, own `[workspace]`). V1 and V2 were deliberately naive baselines
(u128 `%` multiplication, Euclid inversion, schoolbook polynomials); V3 replaced the arithmetic
underneath (Montgomery fields up to 512 bits with an assembly multiplier, binary-GCD inversion,
Karatsuba, product/remainder trees, GF(2ⁿ) with carry-less multiplication) and added the missing
algorithm families. Nothing here is constant-time or intended for production.

**What is covered, what is not, and why: [`docs/SURVEY.md`](docs/SURVEY.md)** (status of every
algorithm in the literature I could identify, with the module and the test that checks it).

## Implemented

V1 = first delivery, V2 = kernel-only/BMSS extension, V3 = speed and coverage round. Different
problems, different algorithms; they are benchmarked per problem, not against each other.

| Problem | Algorithms (module under `src/`) |
|---|---|
| kernel → isogeny | Vélu, odd (V1) and any subgroup (V2) `kernel/velu.rs` · Kohel, odd (V1) and even (V2) `kernel/kohel.rs` · x-only Vélu `kernel/xonly.rs` · Montgomery x-only Vélu, affine (V2) and projective (A24:C24) (V3) `kernel/montgomery.rs` · ℓⁿ chains with naive / balanced / cost-model strategies `kernel/chain.rs` (V2) · **√élu with product and remainder trees, Weierstrass `kernel/sqrt_velu.rs` and Montgomery `kernel/sqrt_velu_mont.rs`** · **twisted Edwards and Huff isogenies (Moody–Shumow)** `kernel/models.rs` · **radical 3- and 5-isogenies (Castryck–Decru–Vercauteren)** `kernel/radical.rs` · **twisted Hessian isogenies** `kernel/hessian.rs` · **Montgomery 2-/4-isogeny 2ᵉ chains with optimal strategies (SIKE formulas)** `kernel/two_power.rs` · **char-2 Vélu and Kohel** `binary.rs` (V3) |
| (E, Ẽ) → isogeny | Padé on the ℘-series (V1) `find/elkies.rs` · the BMSS family: linear algebra, Stark, Atkin, Atkin + modular composition, Elkies 1992, Elkies 1998, fastElkies, fastElkies′, and σ from Φ's second derivatives (V2) `find/bmss.rs` |
| (E, ℓ) → ℓ-isogenies | Φ_ℓ roots + Elkies codomain (V1) · division-polynomial factoring (V1; char 2 V3) · **Φ_ℓ by Hecke operators and Newton's identities; integer Φ_ℓ by CRT, Φ_ℓ mod 2** `find/modpoly.rs` (V3) · **Schoof–Elkies–Atkin point counting with isogeny cycles (t mod ℓᵏ) and BSGS recombination** `find/sea.rs` (V3) |
| (E₁, E₂) → isogeny | Galbraith BFS, GHS, Kohel volcano walk, Couveignes MITM (ordinary), Delfs–Galbraith (V1) · Galbraith–Stolbunov weighted walk, CSIDH action + MITM + ideal orders (V2) `path/*` · **CSIDH-512 with CLMPR batching and a projective tree strategy; relation lattice / class-group structure; paths on binary curves through any neighbour oracle** (V3) `path/csidh.rs`, `path/relation.rs`, `path/graph.rs` |
| supersingular, endomorphism side (V3) | **B_{p,∞}, O₀, ideals, LLL, Fincke–Pohst** `quat/mod.rs` · **KLPT (ℓ = 2)** `quat/klpt.rs` · **class sets, Brandt matrices, Mestre's graph, Eichler mass formula** `quat/brandt.rs` · **Deuring correspondence ideal ↔ kernel over F_{p⁴}** `quat/deuring.rs` |
| genus 2 (V3) | **Richelot (2,2)-isogenies (codomain, points), splitting J(C) → E₁ × E₂, gluing E₁ × E₂ → J(C), Igusa–Clebsch invariants, superspecial Richelot graph** `genus2.rs` · **theta-model (2,2)-isogenies, theta gluing, split detection, theta doubling, Kani-lemma (2ᵃ, 2ᵃ)-chains** `theta.rs` |
| higher dimension (V3) | **general-dimension theta (2,…,2)-chains and dimension-4 Kani embeddings** (auxiliary α ∈ M₂(ℤ[i]) of any degree via four squares; split decision) `theta_g.rs` |
| auxiliary | dual isogeny `find/dual.rs`, Kohel's End(E) conductor `path/endo.rs` (V2) |
| arithmetic (V3) | `fpm.rs` Montgomery F_p (1–8 limbs, MULX/ADCX/ADOX assembly for 512 bits, Pornin inversion) · `fp2.rs` F_{p²} over any of them · `gf2n.rs` GF(2ⁿ) · `ext.rs` F_{p⁴} · `int.rs`, `bigint.rs` big integers · `poly.rs`, `series.rs` Karatsuba, Newton |

## Correctness checks (`cargo test --release`: 86 tests, all pass)

Each algorithm is checked against an independent computation, not only against itself. From V1/V2:
Vélu = Kohel = √élu = x-only = Montgomery on common kernels; all eight BMSS methods reproduce Kohel's
kernel polynomial for ℓ ≤ 101; dual∘φ = [ℓ]; SIDH-style key exchange agrees for every chain strategy;
the CSIDH action commutes and every step satisfies Φ_ℓ(j, j′) = 0; paths are verified edge by edge.
Added in V3:

* field arithmetic against plain big-integer arithmetic (30–127-bit, P-256, secp256k1, CSIDH-512,
  2⁵¹² − 569); the assembly multiplier against the portable one on 3000+ products per field;
  binary-GCD inversion against Fermat;
* √élu (both models), Vélu, Kohel and BMSS agree at 256 and 511 bits; Montgomery √élu = Vélu on all
  74 CSIDH-512 primes; CSIDH-512 batched = stepwise = tree strategy, and two keys commute;
* char 2: Φ_ℓ mod 2 neighbours = neighbours from factored division polynomials; paths verified edge by edge;
* KLPT: N(J) = 2ᵉ, J ⊂ O₀, J = Iξ exactly; Brandt: mass formula, tr(Bᵏ) = tr(Aᵏ), and B(2) equals the
  Φ₂ multiplicity matrix entrywise under the Deuring bijection (classes → supersingular j, p = 1259);
* genus 2: Richelot preserves the L-polynomial; L(C) = L(E₁)L(E₂) for splitting and gluing; the
  superspecial graph's vertex counts equal Ibukiyama–Katsura–Oort (11 primes) and h(h+1)/2;
* Hecke/Newton Φ_ℓ = linear-algebra Φ_ℓ (ℓ ≤ 23 in tests, ≤ 43 in the bench) and = the CRT integer
  coefficients; SEA = BSGS group order at 40 and 61 bits;
* radical chains on CSIDH-512 = the CSIDH action 𝔩₃ᵏ, 𝔩₅ᵏ; relation-lattice vectors act trivially;
* Edwards and Huff codomains have the j of the Weierstrass/Montgomery Vélu codomain; images lie on
  the codomain; the maps are homomorphisms; Hessian isogenies (ℓ = 5..17) send the kernel to the
  identity, land on the codomain given by the derived formula, are homomorphisms, and the codomain
  has the same point count;
* theta: gluing = Howe–Leprévost–Poonen gluing (Igusa–Clebsch invariants); every theta step is one
  of the 15 Richelot neighbours; Kani chains split exactly at E₀ × X with X computed by Vélu, and a
  twisted isotropic kernel does not split; with an endomorphism auxiliary isogeny (126-bit p) the chain
  splits off j = 1728; the strategy (theta doublings) and push-everything chains give identical codomains;
* theta_g: the general-dimension code with g = 2 reproduces the dimension-2 Kani split; a dimension-4
  embedding of a 3ᵇ-isogeny (auxiliary degree 2ⁿ − 3ᵇ as four squares) splits for the true kernel and
  not for a twisted one;
* isogeny cycles: t mod 3⁴, 5³, 7², 11² equal the BSGS trace; Montgomery 4-isogeny chain = 2-isogeny
  chain = Weierstrass Vélu chain; generic F_{p²} = the u64 F_{p²}.

The benchmark re-verifies each result before timing and stores `verified` in every record.

## Run

```
cargo test --release
cargo run --release --bin bench -- --out results/run.jsonl                 # V1 groups: kernel find path
cargo run --release --bin bench -- v2 --out results/run-v2.jsonl           # V2 groups: kernel2 chain bmss csidh v2path
cargo run --release --bin bench -- --out results/run-v3.jsonl p1kernel big char2 quat genus2 phi radical relation models
cargo run --release --bin bench -- --out results/run-v3b.jsonl theta twopow sea                # sea = SEA part of phi
cargo run --release --bin bench -- --quick kernel2                         # smoke test of one group
cargo run --release --bin micro -- --out results/micro.jsonl               # field / polynomial / GF(2^n) micro benchmarks
python3 scripts/report.py results/run-v2.jsonl > results/run-v2.md
```

`ISOGENY_NO_ADX=1` disables the assembly multiplier; `ISOGENY_ADX4=1` enables it for 256-bit fields.

## Results (single thread, 4-vCPU Xeon @ 2.1 GHz VM, shared host, no pinning: expect noise)

| file | content |
|---|---|
| `results/baseline*.jsonl`, `*.md`, `env*.txt` | V1 and V2 baselines (see the V2 notes below) |
| `results/micro-*.jsonl` | micro benchmarks: `micro-v2-baseline` (before V3), `micro-p1-step1..4` (V3 steps), `micro-p3` (with GF(2ⁿ) and Karatsuba thresholds) |
| `results/p1-kernel.jsonl`, `p1-csidh.jsonl` | V3 kernel → isogeny old vs new; CSIDH incl. CSIDH-512 before the assembly multiplier |
| `results/p2-big.jsonl` | 256- and 511-bit workloads |
| `results/p3-char2.jsonl` | characteristic 2 |
| `results/p4-quat.jsonl` | KLPT, class sets/Brandt, Deuring |
| `results/p5-genus2.jsonl` | Richelot, gluing, superspecial graph |
| `results/p6-*.jsonl` | Φ_ℓ and SEA, radical isogenies, relation lattice, Edwards/Huff |
| `results/p7-csidh-adx.jsonl` | `csidh` group re-run with the assembly multiplier |
| `results/p8-theta.jsonl`, `p8-models.jsonl` | theta (2,2)-isogenies, Kani chains, Richelot on the same field; `models` with Hessian |
| `results/p8-sea.jsonl`, `p8-sea-walk.jsonl` | SEA with cycle bounds 0/20/40/80, BSGS recombination; the earlier run with the linear walk |
| `results/p8-twopow.jsonl` | 2^216-isogeny chains over the SIKEp434 F_{p²} |

Selected V3 measurements (medians unless stated; all records `verified: true`):

* **CSIDH-512, exponents in [−5, 5]⁷⁴**: 81 ms (CLMPR batched, before the field work; commit message) →
  49.3 ms (CLMPR) → 34.2 ms (projective tree strategy, `p1-csidh.jsonl`) → 32.2 ms (tree strategy with the
  assembly multiplier, `p7-csidh-adx.jsonl`; CLMPR 46.2 ms in the same run).
* **Radical vs Vélu steps on CSIDH-512**: ℓ = 3: 59 µs per step vs 1.16 ms (stepwise CSIDH: fresh
  point, cofactor ladder, Vélu), 19.6×; ℓ = 5: 70 µs vs 987 µs, 14.1×. Radical steps apply only to
  N = 3, 5 here.
* **Montgomery √élu vs projective Vélu, codomain only**: 511 bits ℓ = 1009 1.02×, 4001 1.41×,
  10007 1.62×; 256 bits ℓ = 10007 1.51×; below ℓ ≈ 1000 Vélu is faster (ℓ = 101: 40 vs 62 µs).
  Weierstrass √élu (three power sums) is slower than x-only Vélu at 256 bits for every ℓ measured and
  at 511 bits up to ℓ = 4001; at 511 bits ℓ = 10007 it is faster (13.9 vs 18.6 ms). The commit
  message of `774cdfd6` said "never" and "1.67×"; the recorded run says the above.
* **Field**: CSIDH-512 multiplication 85 ns (portable) / 78 ns (assembly); inversion 75.5 → 8.0 µs;
  512-bit sqrt 370 → 57 µs; P-256 inversion 13.1 → 3.2 µs. 256-bit assembly multiplication was not faster.
* **Φ_ℓ mod p (61 bits), Hecke/Newton vs dense linear algebra**: ℓ = 11 0.56 vs 1.94 ms (3.5×), ℓ = 23
  10.4 vs 86.8 ms (8.4×), ℓ = 31 31 vs 414 ms (13.4×), ℓ = 43 75 ms vs 2.52 s (33.4×); ℓ = 127 6.1 s
  (linear algebra not run). The commit message of `856b162c` quotes 3.8×–36.6× from an earlier run.
* **SEA** with Φ_ℓ precomputed (`p8-sea.jsonl`, one run per setting): with baby-step giant-step on the
  candidate progression (above 1024 candidates) and isogeny cycles up to degree 40, 61-bit curves take
  1.2–2.2 ms (1.3–4.5 ms without cycles; 6–33 ms with the earlier linear walk; BSGS point counting
  175–204 ms); 40-bit 0.3–0.8 ms; 127-bit 0.26–0.90 s with cycles vs 0.28–5.0 s without (3 curves).
  Cycle bound 80 was slower at 61 bits (2.4–6.6 ms). Precomputing Φ_ℓ for ℓ ≤ 89 over the 127-bit
  field took 34 s and is not included.
* **Theta model** (`p8-theta.jsonl`, 50-bit F_{p²}): (2,2) codomain 0.64 µs and image 0.17 µs vs Richelot
  codomain 1.81 µs and point image 13.0 µs on the same field; Kani (2ᵃ, 2ᵃ)-chain including the split
  test 26 µs (a = 8) to 77 µs (a = 16). With an endomorphism γ = u + v·i of E₀ as the auxiliary isogeny
  (no smoothness needed) the split test runs at cryptographic sizes: a = 64 / 126-bit p 1.09 ms, a = 128 /
  261-bit 6.8 ms, a = 200 / 360-bit 15.6 ms (optimal strategy with theta doublings; pushing every multiple
  instead: 3.1 / 28.9 / 87.4 ms; the first version, which also recomputed the multiples by affine doubling,
  took 620 ms at a = 200).
* **2^216-isogeny over the SIKEp434 F_{p²}** (`p8-twopow.jsonl`, 3 points pushed): Montgomery 4-isogeny
  chain 2.67 ms with the optimal strategy (18.6 ms multiply-only, 13.1 ms push-only), 2-isogeny chain
  3.13 ms, affine Weierstrass Vélu chain 20.8 ms; F_{p²} multiplication 259 ns.
* **Hessian** (61-bit, `p8-models.jsonl`): kernel + codomain 0.44–1.70 µs vs Weierstrass Vélu 0.57–3.96 µs
  (ℓ = 5..31; projective, one inversion); evaluation 0.18–0.89 µs (projective output) vs 0.35–0.96 µs (affine).
* **KLPT (ℓ = 2)**, mean of 5: e/log₂p = 4.61 (31-bit p), 4.10 (60), 3.87 (100), 3.79 (128); 12.6–33.2 ms.
* **Class sets of O₀** (BFS on 2-neighbours): p = 3499, 292 classes, 0.93 s; Deuring ideal → curve
  16 ms per class (p = 1259), 30 ms (p = 3499).
* **Genus 2**: Richelot codomain 0.65 µs (61-bit), 3.8 µs (P-256); point image 3.3 µs; superspecial
  Richelot graph at p = 199: 3077 Jacobians + 153 products, 48 412 edges, 4.8 s.
* **Relation lattice**: h = 1905 (20-bit p) to h = 59 617 257 (50-bit p) in 0.04–1.5 s; a uniformly
  random class [𝔩₁]ᵃ evaluated from the Babai-reduced exponent vector: n = 8 primes, ℓ₁ norm 7.6 and
  21 µs vs ℓ₁ norm ≈ 24 106 and 95 ms as 𝔩₁ᵃ.
* **Char 2** (GF(2⁶¹)): Vélu codomain 0.28–1.35 µs for ℓ = 3..13; Φ_ℓ-mod-2 neighbours 21–101 µs vs
  division-polynomial factoring 28 µs–87 ms; BFS with the Φ oracle 16–84× faster than with the kernel oracle
  (instances with ≥ 3 nodes).
* **Edwards (61-bit)**: one-variable X evaluation 0.68–0.91× the time of Weierstrass Vélu evaluation; the
  kernel + codomain step is slower than Weierstrass Vélu from ℓ = 7 (ℓ = 31: 5.9 vs 4.0 µs).

V1/V2 measurements (u64 fields; kept for the record): codomain from a generator, 32-bit p: Vélu 287 ns (ℓ = 3)
to 144 µs (ℓ = 1009); (E, Ẽ) → isogeny at ℓ = 101, 40 bits: Elkies 1998 74 µs … V1 Padé 3.62 s; ℓ² chain
e = 24 over F_{p²}: naive 263 µs, balanced 81 µs, cost model 96 µs; GHS / Galbraith–Stolbunov / BFS on five
32-bit instances 2.4–17.5 ms. V2 re-runs: the CSIDH ladder ran 128 iterations whatever the exponent
(`baseline-v2-csidh-rerun.jsonl`); the chain cost model charged n instead of n − k multiplications per split
(`baseline-v2-chain-rerun*.jsonl`).

## Known gaps and limitations (details in `docs/SURVEY.md`)

* **Not implemented** (implementation gaps): reading the secret isogeny off the embedding (Robert's
  evaluation of F on torsion → full key recovery; needs extracting points from a generic split codomain),
  the Castryck–Decru digit-guessing recovery for 2ᵃ < 3ᵇ, dimension-8 embeddings exercised on a concrete
  instance, (ℓ,ℓ)-isogenies for odd ℓ, SQIsign2D. The Kani *split decision* is implemented in dimension 2
  (up to 360-bit p) and dimension 4 (arbitrary auxiliary degree via four squares); KLPT output → isogeny at cryptographic size (Deuring here needs u64 p and torsion
  over F_{p⁴}); isogeny cycles for Atkin primes; Couveignes 1996 p-torsion, Lercier, Lercier–Sirvent for
  (E, Ẽ) → isogeny in small characteristic; radical isogenies other than N = 3, 5; Jacobi-quartic models;
  characteristic 3; Sutherland-style Φ_ℓ. Not implementable here (resource limit): quantum algorithms.
* The BMSS methods and √élu use Karatsuba, not FFT multiplication, so the papers' M(ℓ) bounds are not reached.
* The assembly multiplier gains 9 % at 512 bits and nothing at 256 bits; the CSIDH-512 action is variable-time.
* KLPT is for left O₀-ideals with ℓ = 2 and p ≡ 3 mod 4; e/log₂p ≈ 3.8 at 128 bits, above the ≈ 3.5 heuristic.
* Hecke/Newton Φ_ℓ needs char > ℓ + 1; its cost grows quickly (ℓ = 127 in 6.1 s at 61 bits, ℓ ≤ 89 in 34 s
  at 127 bits); SEA is therefore benchmarked with Φ_ℓ precomputed. The SEA stopping rule and cycle bound
  are tuned on 13 curves; the per-curve spread is large.
* Theta: formulas derived here from the duplication formula and validated against the Mumford-side code;
  the Rosenhain formula as remembered had t₁ and t₃ exchanged in μ and ν (the test caught it).
* Relation lattices are computed for 20–50-bit p only (class-group computation is BSGS, not subexponential).
* From V1/V2, unchanged: Elkies/BMSS reject j ∈ {0, 1728}; Couveignes (ordinary) uses planted exponents;
  Kohel volcano paths stay in one ℓ-volcano; GHS with volcano normalisation failed on 3 of 15 conductor > 1
  instances; Delfs–Galbraith returns a j-path, needs p ≡ 7 mod 8, and found no path in 2 of 5 38-bit instances;
  ℓ = 2 dividing the conductor of Z[π] is unsupported in End(E); path instances are small, so do not read
  scaling exponents from them.
* No Conductor task could be obtained (the control plane needs Docker/Postgres, unavailable here); all
  edits are confined to `isogeny_algos/`.
