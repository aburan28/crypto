#set document(title: "Small-x factor base on NIST P-192: exact two-point relation yield")
#set page(paper: "a4", margin: (x: 2cm, y: 1.8cm), numbering: "1 / 1")
#set text(size: 9.5pt)
#set par(justify: true, leading: 0.6em)
#set heading(numbering: "1.")
#show heading.where(level: 1): set text(size: 12pt)
#show table: set text(size: 8.5pt)
#show raw.where(block: true): set text(size: 8pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[Small-x factor base on NIST P-192: exact two-point relation yield]

  #text(size: 9pt)[Study `research/prime_small_x_p192_20261014` · run 2026-10-14 · seed 7 · `ic prime` · `src/cryptanalysis/prime_small_x.rs`]
]

= Question and answer

*Question.* What is the probability, per trial, that a random point of NIST P-192 is a signed sum of two points of a small-x factor base, and how does the measurement compare with the estimates given before it was run?

*Answer.* With a base of B abscissae (2B points) the probability is *exactly 2B²/(N − 1)* on P-192, where N ≈ 2¹⁹² is the group order: 2⁻¹⁷⁵·⁰ at B = 256, 2⁻¹⁷¹·⁰ at 1024, 2⁻¹⁶⁹·⁰ at 2048, 2⁻¹⁶⁷·⁰ at 4096 and 2⁻¹⁶⁵·⁰ at 8192, the largest base the tool builds. This is a count, not an estimate: every one of the up to 67 million pair sums was computed and none coincide, so no collision removes any. The expected number of trials per relation is the reciprocal, 2¹⁶⁵ or more, and 2¹⁷⁸ trials to collect B + 1 relations at B = 8192. Pollard rho costs about 2⁹⁶ group operations. Monte Carlo on P-192 finds nothing, as it must, because the yield is 2¹⁴⁶ times smaller than the 2M-trial run can resolve. The method is checked on toy primes, where the count equals the Monte Carlo rate to within 0.5% and the detector fires on every planted decomposition.

*Against the estimates.* The earlier estimate for a 2¹³-point base was a yield near 2⁻¹⁶⁶ and about 2⁻¹⁴⁴ expected collisions. The measured yield for 8192 points (B = 4096) is 2⁻¹⁶⁷·⁰, a factor of two lower than estimated. The estimate counted unordered pairs of points (2²⁵), but a sum and its negative have the same abscissa, so those pairs give 2²⁴ distinct abscissae. The collision estimate stands: about 2⁻¹⁴⁵ expected, 0 observed, which cannot discriminate. The estimates for a hypothetical 2⁴⁰-point base and for three-point sums were not measured.

= Method

A small-x base takes the points with the smallest abscissae x = 0, 1, 2, … that lie on y² = x³ − 3x + b, one sign of each, and uses both signs. Let the base have B abscissae. A signed sum of two base points has abscissa x(a + b) or x(a − b) for a pair of distinct base abscissae, or x(2a). There are B(B − 1)/2 pairs with two sums each and B doublings, which is B² abscissae if no two coincide. A target T is a signed sum of two base points exactly when x(T) is one of them, and then so is −T. The number of covered targets is therefore twice the number D of distinct sum abscissae, and a uniformly random group element decomposes in one trial with probability exactly

#raw("yield = 2·D / (N − 1)        and        D ≤ B², with equality when no two sums collide", block: true)

`ic prime` computes D from every pair sum. The abscissa of a ± b is λ² − x₁ − x₂ with λ = (y₁ ∓ y₂)/(x₁ − x₂), and one inversion is shared by a whole row (Montgomery's trick), about seven field multiplications per pair. The keys are exact 192-bit field elements, sorted and deduplicated, so D is exact. P-192 uses the repository's native three-limb Montgomery field; toy curves use a 64-bit field. Because the count is exact, no target is sampled.

Three checks run with every measurement. *Monte Carlo*: random group elements (a random abscissa with a square right-hand side) are looked up in the key set, and the first hits are recovered as real decompositions by subtracting base points. *Planted controls*: a random signed sum of two base points is built by the group law and must be detected by the same lookup and recovered. *Toy curves*, where hits are common: y² = x³ − 3x + b over a prime p ≡ 3 (mod 4) just below 2ᵇⁱᵗˢ, with b chosen so the cubic has no root (no rational 2-torsion, so every abscissa carries two points) and the order counted exactly by walking the field.

= Results

== NIST P-192

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (right, right, right, right, right, right, right, right, right),
    table.hline(),
    [*B*], [*points*], [*pair sums*], [*collisions*], [*exact yield*], [*trials per relation*], [*B + 1 relations*], [*pair sums, s*], [*MC µs/trial*],
    table.hline(stroke: 0.5pt),
    [256], [512], [2¹⁶], [0], [2⁻¹⁷⁵·⁰⁰], [2¹⁷⁵·⁰⁰], [2¹⁸³·⁰¹], [0.01], [9.18],
    [1024], [2048], [2²⁰], [0], [2⁻¹⁷¹·⁰⁰], [2¹⁷¹·⁰⁰], [2¹⁸¹·⁰⁰], [0.10], [9.17],
    [2048], [4096], [2²²], [0], [2⁻¹⁶⁹·⁰⁰], [2¹⁶⁹·⁰⁰], [2¹⁸⁰·⁰⁰], [0.39], [9.44],
    [4096], [8192], [2²⁴], [0], [2⁻¹⁶⁷·⁰⁰], [2¹⁶⁷·⁰⁰], [2¹⁷⁹·⁰⁰], [1.70], [10.11],
    [8192], [16384], [2²⁶], [0], [2⁻¹⁶⁵·⁰⁰], [2¹⁶⁵·⁰⁰], [2¹⁷⁸·⁰⁰], [6.01], [9.59],
    table.hline(),
  ),
  caption: [P-192 small-x bases (largest abscissa x = 508 to 16535). Exact yield is 2D/(N − 1), identical to the collision-free model 2B²/(N − 1) in every row. "B + 1 relations" is log₂((B + 1)/yield). Monte Carlo: 2 000 000 trials, 0 hits in each run; wall microseconds per trial on 4 threads, about 4× that on one thread. Planted controls: 32 of 32 detected and 32 of 32 recovered in every run. Source: `data/p192_*.json`.],
)

On P-192 the Monte Carlo is a null result: 0 hits in 2 000 000 trials bounds the rate above by 1.84 × 10⁻⁶ (exact one-sided 95% bound), and the exact yield is 2⁻¹⁶⁵ at best. It cannot test the count; the toy runs do. What it does show is the cost of a trial at this size, about 38 µs on one thread for a random-point draw plus a lookup, which excludes the cost of deciding a target by solving the Semaev S₃ system.

== Method check on toy primes

#figure(
  image("graphs/fig1_yield.svg", width: 100%),
  caption: [Left: Monte Carlo rate against the exact count (dots, with 95% intervals) and against the collision-free model (crosses), for 17 toy runs with 3 000 000 trials each. Right: P-192, base size against log₂ trials. Run 2026-10-14, seed 7. Source: `data/*.json`; plot `graphs/plot_prime.py`.],
)

In all 17 toy runs the exact count lies inside the 95% interval of the Monte Carlo rate, and the ratio of the two is between 0.993 and 1.002 (the spread is what the number of hits predicts: 5800 hits resolves 1.3%). The collision-free model is not close where collisions are common: at yield 0.87 (B = 4096, N ≈ 2²⁴) the model gives 2.0, and at yield 0.12 it is 5 to 6% high. That is the regime the exact count exists for. On P-192 the two coincide to every digit, because with B ≤ 8192 the expected number of collisions among the sums is below 2⁻¹⁴⁰. In every toy run, 32 of 32 planted decompositions were detected and recovered, and the first 8 Monte Carlo hits were recovered as real decompositions.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (right, right, right, right, right, right, right, right),
    table.hline(),
    [*N*], [*B*], [*collisions*], [*exact*], [*model*], [*MC rate*], [*hits*], [*MC / exact*],
    table.hline(stroke: 0.5pt),
    [65173], [16], [0], [0.007856], [0.007856], [0.007852], [23 557], [0.9995],
    [65173], [32], [17], [0.03090], [0.03142], [0.03085], [92 563], [0.9984],
    [65173], [64], [198], [0.1196], [0.1257], [0.1197], [358 976], [1.0003],
    [65173], [128], [3384], [0.3989], [0.5028], [0.3988], [1 196 271], [0.9995],
    [65173], [256], [37 316], [0.8660], [2.011], [0.8659], [2 597 822], [0.9999],
    [1 047 267], [32], [3], [0.001950], [0.001956], [0.001940], [5821], [0.9951],
    [1 047 267], [128], [298], [0.03072], [0.03129], [0.03077], [92 301], [1.0015],
    [1 047 267], [512], [55 745], [0.3942], [0.5006], [0.3943], [1 182 896], [1.0003],
    [1 047 267], [1024], [594 877], [0.8664], [2.003], [0.8663], [2 598 783], [0.9998],
    [16 777 897], [128], [15], [0.001951], [0.001953], [0.001938], [5814], [0.9932],
    [16 777 897], [1024], [63 063], [0.1175], [0.1250], [0.1172], [351 602], [0.9976],
    [16 777 897], [2048], [895 288], [0.3933], [0.5000], [0.3928], [1 178 290], [0.9987],
    [16 777 897], [4096], [9 523 333], [0.8647], [2.000], [0.8649], [2 594 613], [1.0002],
    table.hline(),
  ),
  caption: [Thirteen of the 17 toy runs; all 17 are in `data/`. Yields are probabilities per trial; the model is 2B²/(N − 1). Source: `data/toy*.json`.],
)

= Checks

- *Brute force.* On a 14-bit toy curve, every sum of two signed base points by the full group law reaches exactly twice as many distinct points as there are distinct sum abscissae (unit test).
- *Independent arithmetic.* On P-192, every batched sum abscissa of a ± b for the first six base points is found among the abscissae of the native Jacobian addition of the same points (unit test).
- *Monte Carlo.* On toy primes, 17 of 17 exact counts lie inside the 95% interval.
- *Controls.* Planted decompositions were detected and recovered in all 22 runs.
- *Counts.* The toy group orders equal the number of points found by walking the field, plus the point at infinity (unit test). The P-192 order is the standard prime N.
- *Inspection.* `ic inspect --curve secp192r1` (also `ic secp192r1`) passes the screens for the order, the modulus, nonsingularity and the generator, and the generator's order. Its cardinality check is reported as not implemented for prime curves.

= Limits

+ *Two-point sums only.* A relation here is a signed sum of two small-x points. Three-point sums, which are what a prime-field attack would use, are not measured, and neither is the 2⁴⁰-point base of the earlier estimate.
+ *Yield, not solving.* The per-trial cost of deciding a target (a Semaev S₃ root find over the base) is not part of this measurement, and no linear algebra or descent runs. A relation yield is not a discrete logarithm, and the 2¹⁷⁸ trials in the last column are a count of trials, not a time.
+ *One base family.* Only the small-x base exists. The tool says so. Other families (Kummer-tower, Dickson, structured bases) are not offered for prime fields.
+ *Size.* The tool builds at most 8192 abscissae because the pair sums are held in memory (67 million 192-bit keys at the maximum). At B = 8192 the yield is 2⁻¹⁶⁵, and the base would need about 2⁹⁶ abscissae for 2B²/N to reach one.
+ *Toy curves are not prime order.* They have odd, not necessarily prime, order and no rational 2-torsion. They validate the counting method and the Monte Carlo, not any property of P-192.
+ *The greedy gain.* The coverage-driven selector of the binary-field study chooses among Frobenius orbits; P-192 has none, and with zero collisions at every base size here a coverage-driven choice has nothing to improve: any base of the same size has the same yield, 2B²/(N − 1), provided its sums are distinct.

= Reproduction

```sh
cargo run --release --bin ic -- --json prime --size 4096 --trials 2000000 --controls 32 --seed 7
cargo run --release --bin ic -- --json prime --curve toy:20 --size 512 --trials 3000000 --controls 32 --seed 7
cargo run --release --bin ic -- inspect --curve secp192r1
python3 research/prime_small_x_p192_20261014/graphs/plot_prime.py
```

`data/COMMANDS.txt` lists every run. Files: `data/*.json` (one `ic prime` report per run, 22 runs), `graphs/plot_prime.py`, `graphs/fig1_yield.{svg,png}`.
