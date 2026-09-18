# Pair-enumeration index calculus on one G7e GPU

**Question.** Can the cheapest three-summand decomposition oracle this
repository has — exhaustive pair enumeration — be made, by running it on
an RTX PRO 6000 Blackwell (`g7e.2xlarge`), into an index-calculus attack
on ECC2K-130 that beats Pollard rho?

**Answer.** No. The GPU is a constant in front of the same product law
already priced in `RESEARCH_ECC2K130_DECOMPOSITION.md`. Wall-clock moves;
the ratio to the floor does not. The work is an engineering calibration
of that oracle in the packed type-II normal-basis arithmetic the rho
client already uses, plus verified planted decompositions and two toy
discrete logs recovered by the same algorithm on the CPU. CryptoMiniSat
on this host's CPU recovers the same toy logs and loses to pair
enumeration on planted triples; SAT internals stay uncalibrated and do
not enter `S`.

## 1. The boundary, stated before measuring

Per `AGENTS.md` §1 this thread has both kinds of boundary.

The instance is `K_0 : y² + xy = x³ + 1` over `F_2^131`, cofactor 4,
prime subgroup order

```
r = 680564733841876926932320129493409985129,   log2 r = 129.0
#E = 4r ≈ 2^131
```

**The reference** is Pollard rho with the `⟨−1⟩ × ⟨π⟩` speed-up on the
same subgroup: `2^60.809` walk iterations. On this SKU the shipping
packed client charges **5.3125 field products per iteration**
(`ecc2k130/THROUGHPUT-29B.md`). In the product unit that is

```
rho products = 5.3125 × 2^60.809 ≈ 2^63.21
S_rho        = (5.3125 × 2^60.809) / 2^64.5 ≈ 0.411
```

**The floor** is the product law for an `m`-summand factor base `F`.
A random target in a group that contains every sum of `m` base points
decomposes with probability `C(|F|, m) / #G`. Streaming pair
enumeration spends `C(|F|, 2)` affine additions on every trial. Charging
`|F|` relations (one logarithm, no Frobenius collapse) therefore costs,
at leading order in `|F|`,

```
|F| · (#G / C(|F|, 3)) · C(|F|, 2) · 20   field products
    ≈  60 · #G
```

independent of `|F|`. For `G = E(F_2^131)` that is `60 · 2^131 ≈ 2^136.91`
field products, `2^73.70` times the rho product count. Applying the
Frobenius collapse (`|F|/131` orbit unknowns) divides by `131 ≈ 2^7.03`
and leaves `2^66.66` times rho — still not a crossover, and still a
counting identity rather than a measured runtime.

Affine addition is priced at **10 field products**: eight for the
Itoh–Tsujii inverse `inv131` and two for the slope and the `y`
coordinate. Each pair probe does `P_i + P_j` and `R − S`, so 20
products. Squares and Frobenius permutations are not charged; that
favours the oracle.

**Falsification target, stated in advance.** This thread would be a
success if some GPU kernel answered "does `R` decompose over this `F`?"
in fewer than `2^{70}`-th of the `C(|F|, 2)` pair cost, on a base with
`m · log2 |F| ≥ 131`, with every hit verified on the curve, and with the
whole logarithm (collection, linear algebra, verification) below
`2^60.809` walk-equivalent products. It is abandoned if measured yields
match `C(|F|, 3)/r` and the kernel still visits a constant fraction of
all pairs.

Inadmissible: quoting pairs/second as an attack cost, dropping the
second addition (`R − S`), charging a pair table that does not fit in
the 96 GiB device, changing `m` or `n`, or converting units after the
run.

## 2. What the GPU is for

`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md` §7 already showed that
brute-force pair enumeration is cheaper than the Riemann–Roch SAT
encodings on the real curve, at dimension 6. Algebra does not beat
search here. The only remaining question that a `g7e` can answer is
whether search, run in the same field arithmetic as rho, becomes
*practically* interesting.

It does not, and the reason is the floor above, not occupancy. A 17 B/s
rho client on this SKU is a 17e9-iteration/s walk. Pair enumeration is
inverse-heavy affine addition. Even granting the walk's iteration rate
as an upper bound on pair-adds — which it is not, because a walk step
is 5.3125 products and an affine add is 10 plus the second add — a
complete logarithm still needs `~2^131` probes.

What the GPU *can* do, and what this thread measures:

1. Build the Hamming-weight-`≤ 2` even-trace factor base in the packed
   ONB (the same coordinates the walk already popcounts).  Even weight is
   the necessary odd-order condition `Tr(x) = 0`; this thread does not
   spend a `[r]P` certificate on every base point.
2. Recover planted 3-sums from that base, on device, checked against
   the independent `Ref` adder.
3. Scan a natural subgroup point (the client generator) to completion
   and confirm a hit count consistent with `C(|F|, 3)/r ≈ 2^{-90}`.
4. Report affine-adds/s and pairs/s on this SKU, as a practicality
   column, not as `S`.
5. Recover planted discrete logs at `n = 5` and `n = 9` by the same
   pair-enumeration algorithm on the CPU, so the oracle is not a GPU
   artefact.

Weight 3 is admitted by the binary (`--weight 3`) and is not required
for the verdict: it enlarges `|F|` and does not change the leading
`60 · #E` product count.

## 3. One table, one unit

Unit: **field products**. Conversion from rho iterations: 5.3125
products/iteration, measured on this SKU. `S = products / √r` with
`√r = 2^{64.5}`. Frozen receipt:
`ecc2k130/benchmarks/indexcalc-g7e/summary.json`.

| Variant | log2 products | S | vs rho | vs floor | Correctness | Class |
|---|---:|---:|---:|---:|---|---|
| Pollard rho `⟨−1⟩×⟨π⟩` *(reference)* | 63.22 | 0.411 | `2^+0.00` | — | shipping client | baseline |
| Product-law floor, `m = 3`, any `\|F\|` | 136.91 | `2^{72.4}` | `2^+73.69` | `2^+0.00` | derived | floor |
| + Frobenius collapse `\|F\|/131` | 129.87 | `4.78×10^{19}` | `2^+66.66` | `2^+0.00` | derived | accounting |
| Streaming pair enum, G7e, `\|F\|=8384` | 129.87 | `4.78×10^{19}` | `2^+66.66` | `2^+0.00` | 8/8 planted; 0 hits on the generator, matching `3.6×10^{-29}` | engineering |

The measured row is the Frobenius-accounting row, run. Its ratio to the
floor is 1. The GPU does not appear in the product column.

**Practicality, not the metric.** On this RTX PRO 6000 Blackwell
(`sm_120`, 97,252 MiB) the occupied pair scan did `3.5141536×10^7`
pairs in 12.73 ms (`2.76×10^9` pairs/s). The same scan on the eight
host cores took 40.23 s, a **3161×** wall-clock speedup and an
engineering constant. At that occupied rate a streaming logarithm is
still `2^{94.2}` seconds. A 2048-point add microbench, too small to
fill the device, measured only `2.52×10^8` affine adds/s and is not
the rate used above.

SAT does not get a row in that table. Its conflict counts have no
measured conversion to field products, so putting them in `S` would be
relabelling. The diagnostic below is wall-clock on this host, frozen in
`ecc2k130/benchmarks/indexcalc-g7e/sat.json` (pycryptosat 5.14.7,
python-sat 1.9.dev15, PEP 668 venv). Type-II ONB exists
only for degrees with `2m+1` prime and `ord_{2m+1}(2) ∈ {m, 2m}`; the
ladder is 5, 9, 11, 23.

| Cell | Result | Correctness | Class |
|---|---|---|---|
| Growth `m=5`, `k=3`, `w=3` | 5/5 solved, median 0.0018 s, 507 gates | planted, 0 unsat, 0 spurious | engineering |
| Growth `m=9` | 5/5 solved, median 0.096 s, 1587 gates | planted, 0 unsat, 0 spurious | engineering |
| Growth `m=11` | 5/5 solved, median 1.25 s, 2311 gates | planted, 0 unsat, 0 spurious | engineering |
| Growth `m=23` | 0/1, budget at 30.05 s, 8029 gates | timeout, 0 unsat | engineering |
| SAT DLP `n=5` | complete in 22 ms, scalar 8 | `[k]P = Q` | correctness |
| SAT DLP `n=9` | complete in 1.07 s, scalar 112 | `[k]P = Q` | correctness |
| Pair vs SAT, `n=9`, `w=2`, 16 planted triples | pair 16/16, median 0.362 ms; SAT 14/16 lifted, median 82.5 ms | pair hits verified; 2 SAT models spurious | engineering |

Unbounded SAT is about **228×** slower than pair lookup on the same
planted triples, and two SAT models failed to lift. The e2e DLP oracle
uses a tighter 0.25 s / 2000-conflict budget; it still recovered the
two toy logs because yield at these sizes is high. Neither fact moves
the `n=131` product-law floor. Ratio to the floor remains 1.

## 4. How to run

On the G7e developer host:

```
make -C ecc2k130 indexcalc-cuda ARCH='-gencode arch=compute_120,code=sm_120'
make -C ecc2k130 test-indexcalc-pairs
make -C ecc2k130 test-indexcalc-cuda
python3 ecc2k130/benchmarks/indexcalc-g7e/run.py
python3 -m venv ~/ic-venv
~/ic-venv/bin/pip install -r ecc2k130/benchmarks/indexcalc-g7e/requirements-sat.txt
~/ic-venv/bin/python ecc2k130/codegen/testdecomp.py
~/ic-venv/bin/python ecc2k130/benchmarks/indexcalc-g7e/run_sat.py
```

`--skip-gpu` on `run.py` still recovers the degree-5 and degree-9
logs. `--from-raw` rebuilds `summary.json` from a frozen
`raw-gpu.json`. The CUDA self-test is a packed-add differential against
`Ref` and one planted triple at weight 2. `run_sat.py` writes `sat.json`
and patches `summary.json['sat']`; it does not use the GPU. System
Python on this host is PEP 668, so the venv is required.

## 5. Classification

| change | class | why |
|---|---|---|
| Pair enumeration vs SAT/RR on the same base | engineering | already measured in the RR panel; GPU repeats it |
| Packed ONB affine add on sm_120 | engineering | same algorithm, different device |
| Toy DLP at `n = 5, 9` (pairs and SAT) | correctness | known-answer logs, not an ECC2K-130 result |
| SAT growth 5/9/11 and timeout at 23 | engineering | same Semaev encoding, host CPU; floor flat |
| Unbounded SAT vs pair lookup at `n=9` | engineering | SAT 228× slower; 2/16 spurious; floor flat |
| Quoting pairs/second or SAT seconds as `S` | relabelling | forbidden; wall-clock is a footnote |

Nothing in this thread is an advance against the floor.

## 6. What this is not

Not a discrete logarithm of the Certicom challenge. Not an F₄ or
Riemann–Roch solver. SAT *was* run on this host; it does not beat pair
enumeration and it does not enter the product unit. Not a pair *table*:
96 GiB holds a weight-2 table and not a weight-9 table, and weight 2
has expected yield `2^{-90}` per target. Not a claim that G7e capacity
was previously unavailable for rho — the rho client already runs here;
this is the first IC oracle that uses the same field code on the same
chip.
