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
discrete logs recovered by the same algorithm on the CPU.

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

## 4. How to run

On the G7e developer host:

```
make -C ecc2k130 indexcalc-cuda ARCH='-gencode arch=compute_120,code=sm_120'
make -C ecc2k130 test-indexcalc-pairs
make -C ecc2k130 test-indexcalc-cuda
python3 ecc2k130/benchmarks/indexcalc-g7e/run.py
```

`--skip-gpu` on `run.py` still recovers the degree-5 and degree-9
logs. `--from-raw` rebuilds `summary.json` from a frozen
`raw-gpu.json`. The CUDA self-test is a packed-add differential against
`Ref` and one planted triple at weight 2.

## 5. Classification

| change | class | why |
|---|---|---|
| Pair enumeration vs SAT/RR on the same base | engineering | already measured in the RR panel; GPU repeats it |
| Packed ONB affine add on sm_120 | engineering | same algorithm, different device |
| Toy DLP at `n = 5, 9` | correctness | known-answer logs, not an ECC2K-130 result |
| Quoting pairs/second as `S` | relabelling | forbidden; wall-clock is a footnote |

Nothing in this thread is an advance against the floor.

## 6. What this is not

Not a discrete logarithm of the Certicom challenge. Not a SAT, F₄, or
Riemann–Roch solver. Not a pair *table*: 96 GiB holds a weight-2 table
and not a weight-9 table, and weight 2 has expected yield `2^{-90}` per
target. Not a claim that G7e capacity was previously unavailable for
rho — the rho client already runs here; this is the first IC oracle
that uses the same field code on the same chip.
