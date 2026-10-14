#set document(title: "Small-x factor base on the SECG curves: exact two-point relation yield")
#set page(paper: "a4", margin: (x: 2cm, y: 1.8cm), numbering: "1 / 1")
#set text(size: 9.5pt)
#set par(justify: true, leading: 0.6em)
#set heading(numbering: "1.")
#show heading.where(level: 1): set text(size: 12pt)
#show table: set text(size: 8pt)
#show raw.where(block: true): set text(size: 8pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[Small-x factor base on the SECG curves: exact two-point relation yield]

  #text(size: 9pt)[Study `research/secg_small_x_20261014` · run 2026-10-14 · seed 7 · `ic prime` · `src/cryptanalysis/{secg_standards,binary_small_x,prime_small_x}.rs`]
]

= Question and answer

*Question.* For the SECG curves `sect113r1`, `sect113r2`, `sect163r1`, `sect163r2` and the other SECG curves up to about 240 bits (among them the 160-bit `secp160k1/r1/r2`), what is the probability, per trial, that a random element of the prime-order subgroup is a signed sum of two points of a small-x factor base?

*Answer.* With a base of B abscissae (2B points) the probability is *2B²/(n − 1)*, where n is the order of the prime-order subgroup. It is a count, not an estimate: every pair sum was computed and, on 18 of the 21 curves, none coincide. At B = 4096 the yield is 2⁻⁸⁷·⁰ on `sect113r1` and `sect113r2`, 2⁻¹³⁷·⁰ on `sect163r1` and `sect163r2`, and 2⁻¹³⁵·⁰ on the three 160-bit prime curves. The trials needed to collect B + 1 relations are about n/(2B): 2⁹⁹·⁰ on `sect113r1` (B = 4096), 2⁹⁸·⁰ at B = 8192, the largest base the tool builds; and 2¹⁴⁹·⁰ on `sect163r2`, 2¹⁴⁷·⁰ on `secp160r1`. Pollard rho costs about √n group operations, 2⁵⁶·⁰ for `sect113r1`, 2⁸¹·⁰ for `sect163r1/r2` and 2⁸⁰·⁰ for the 160-bit primes. The smallest gap to rho is 2⁴¹·⁹ (`secp112r2`, B = 4096) and the gap on `sect113r1` is 2⁴³·⁰ at B = 4096 and 2⁴²·⁰ at B = 8192, so two-point small-x sums lose to rho on every curve measured: the trials for a full system are n/(2B) against √n, so they could reach rho only for a base of about √n/2 abscissae, as large as the rho computation itself. Monte Carlo finds no hit in 66 million random subgroup elements over the 66 SECG runs, as the yields require. The method is checked on toy binary curves where hits are common: the exact count lies inside the 95% interval of the Monte Carlo rate in 14 of 15 runs, and every planted decomposition was detected and recovered in all 81 runs.

*"61 or 63, 160 bit".* The request was read as the 163-bit curves `sect163r1` and `sect163r2` (`sect163k1` is also measured) and, as the alternative reading, the 160-bit curves `secp160k1`, `secp160r1`, `secp160r2`. All are in the tables.

= Method

The base is the points with the smallest abscissae x = 1, 2, 3, … (binary fields: the polynomial whose coefficient bits are the binary digits of x; prime fields: x = 0, 1, 2, …) that lie on the curve, one sign of each, used with both signs. Let the base have B abscissae. A signed sum of two base points has abscissa x(a + b) or x(a − b) for a pair of distinct base abscissae, or x(2a). That is B(B − 1)/2 pairs with two sums each and B doublings, B² abscissae if no two coincide. A target T is a signed sum exactly when x(T) is one of them, and then so is −T. The covered targets are twice the number D of distinct sum abscissae, and a uniform element of the subgroup decomposes in one trial with probability exactly

#raw("yield = 2·D / (n − 1)        and        D ≤ B², with equality when no two sums collide", block: true)

`ic prime` computes D from every pair sum with one batched inversion per row. Binary curves y² + xy = x³ + ax² + b use λ = (y₁ + y₂)/(x₁ + x₂) and x(a + b) = λ² + λ + x₁ + x₂ + a, with −b = (x₂, x₂ + y₂) giving x(a − b), and x(2a) = x² + b/x². Prime curves use the short Weierstrass formulas. Keys are exact field elements, sorted and deduplicated.

*Cofactor.* Every SECG binary curve has cofactor h = 2 or 4, as do `secp112r2` and `secp128r2`. A relation is a relation in the subgroup ⟨G⟩ = [h]E of prime order n, so each small-x point P enters the base as [h]P, and targets are [h] of a uniform curve point. The yield is over ⟨G⟩. The base abscissae are therefore not small after the multiplication. Keeping the raw small-x points and counting only the pairs that sum into ⟨G⟩ would give about 1/h of this yield. The tool does not measure that variant.

*Fields.* Binary: a four-limb F₂ᵐ with carry-less multiplication (PCLMULQDQ where available) and Itoh–Tsujii inversion, checked against the repository's arbitrary-width field. Prime: three-limb Montgomery arithmetic for any odd prime of at most 192 bits with p ≡ 3 (mod 4); P-192 keeps its native field. The constants come from the repository's standards database, `docs/curves/standards/parameters.json`, and the SECG binary constants equal the hand-typed ones in `src/binary_ecc` (unit test).

Checks run with every measurement: *Monte Carlo* on uniform subgroup elements, with the first hits recovered as real decompositions; *planted controls*, random signed sums of two base points, which must be detected and recovered; and, on *toy binary curves* y² + xy = x³ + x² + b over F₂¹³, F₂¹⁷ and F₂¹⁹ with group order 2·n, n prime and counted by walking the field, a comparison of the exact count with the Monte Carlo rate where hits are common.

= Results

== The SECG curves at B = 4096

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right, right, right, right),
    table.hline(),
    [*curve*], [*field bits*], [*h*], [*log₂ n*], [*exact yield*], [*collisions*], [*B + 1 relations*], [*rho*], [*pair sums, s*], [*MC µs/trial*], [*controls*],
    table.hline(stroke: 0.5pt),
    [`secp112r1`], [112], [1], [111.8], [2⁻⁸⁶·⁸], [0], [2⁹⁸·⁸], [2⁵⁵·⁹], [2.89], [4.9], [32/32],
    [`secp112r2`], [112], [4], [109.8], [2⁻⁸⁴·⁸], [0], [2⁹⁶·⁸], [2⁵⁴·⁹], [2.71], [7.7], [32/32],
    [`sect113r1`], [113], [2], [112.0], [2⁻⁸⁷·⁰], [0], [2⁹⁹·⁰], [2⁵⁶·⁰], [2.72], [14.9], [32/32],
    [`sect113r2`], [113], [2], [112.0], [2⁻⁸⁷·⁰], [0], [2⁹⁹·⁰], [2⁵⁶·⁰], [2.87], [15.6], [32/32],
    [`secp128r1`], [128], [1], [128.0], [2⁻¹⁰³·⁰], [0], [2¹¹⁵·⁰], [2⁶⁴·⁰], [3.32], [6.1], [32/32],
    [`secp128r2`], [128], [4], [126.0], [2⁻¹⁰¹·⁰], [0], [2¹¹³·⁰], [2⁶³·⁰], [2.95], [9.8], [32/32],
    [`sect131r1`], [131], [2], [130.0], [2⁻¹⁰⁵·⁰], [0], [2¹¹⁷·⁰], [2⁶⁵·⁰], [3.32], [27.2], [32/32],
    [`sect131r2`], [131], [2], [130.0], [2⁻¹⁰⁵·⁰], [0], [2¹¹⁷·⁰], [2⁶⁵·⁰], [3.30], [26.8], [32/32],
    [`secp160k1`], [160], [1], [160.0], [2⁻¹³⁵·⁰], [0], [2¹⁴⁷·⁰], [2⁸⁰·⁰], [2.81], [8.2], [32/32],
    [`secp160r1`], [160], [1], [160.0], [2⁻¹³⁵·⁰], [0], [2¹⁴⁷·⁰], [2⁸⁰·⁰], [2.68], [8.0], [32/32],
    [`secp160r2`], [160], [1], [160.0], [2⁻¹³⁵·⁰], [0], [2¹⁴⁷·⁰], [2⁸⁰·⁰], [2.63], [7.9], [32/32],
    [`sect163k1`], [163], [2], [162.0], [2⁻¹³⁷·⁰], [12], [2¹⁴⁹·⁰], [2⁸¹·⁰], [3.29], [30.3], [32/32],
    [`sect163r1`], [163], [2], [162.0], [2⁻¹³⁷·⁰], [0], [2¹⁴⁹·⁰], [2⁸¹·⁰], [3.34], [30.1], [32/32],
    [`sect163r2`], [163], [2], [162.0], [2⁻¹³⁷·⁰], [0], [2¹⁴⁹·⁰], [2⁸¹·⁰], [9.28], [32.9], [32/32],
    [`secp192k1`], [192], [1], [192.0], [2⁻¹⁶⁷·⁰], [0], [2¹⁷⁹·⁰], [2⁹⁶·⁰], [2.41], [9.4], [32/32],
    [`secp192r1`], [192], [1], [192.0], [2⁻¹⁶⁷·⁰], [0], [2¹⁷⁹·⁰], [2⁹⁶·⁰], [2.91], [9.1], [32/32],
    [`sect193r1`], [193], [2], [192.0], [2⁻¹⁶⁷·⁰], [0], [2¹⁷⁹·⁰], [2⁹⁶·⁰], [8.01], [29.9], [32/32],
    [`sect193r2`], [193], [2], [192.0], [2⁻¹⁶⁷·⁰], [0], [2¹⁷⁹·⁰], [2⁹⁶·⁰], [6.81], [30.2], [32/32],
    [`sect233k1`], [233], [4], [231.0], [2⁻²⁰⁶·⁰], [20], [2²¹⁸·⁰], [2¹¹⁵·⁵], [8.21], [41.3], [32/32],
    [`sect233r1`], [233], [2], [232.0], [2⁻²⁰⁷·⁰], [0], [2²¹⁹·⁰], [2¹¹⁶·⁰], [7.33], [34.1], [32/32],
    [`sect239k1`], [239], [4], [237.0], [2⁻²¹²·⁰], [18], [2²²⁴·⁰], [2¹¹⁸·⁵], [7.83], [53.4], [32/32],
    table.hline(),
  ),
  caption: [Two-point small-x bases of 4096 abscissae (8192 points), 2²⁴ pair sums. "B + 1 relations" is log₂((B + 1)/yield), "rho" is log₂ √n. Monte Carlo: 1 000 000 trials, 0 hits in every run; wall microseconds per trial on 4 threads. Planted controls detected and recovered. Source: `data/sec*_4096.json`.],
)

#figure(
  image("graphs/fig1_secg.svg", width: 100%),
  caption: [Left: trials for B + 1 relations at B = 4096 (bars) against rho (dashes) on every measured curve; blue is a binary field, orange a prime field. Right: toy binary curves, Monte Carlo rate against the exact count (dots, 95% intervals) and the collision-free model (crosses). Run 2026-10-14, seed 7. Source: `data/*.json`; plot `graphs/plot_secg.py`.],
)

The yield follows 2B²/(n − 1) closely, so a curve's row is set by its subgroup order and the base size, not by its field or family: `sect113r1` and `sect113r2` (n ≈ 2¹¹²) match each other, `sect163r1/r2` (n ≈ 2¹⁶²) match each other, and `secp192r1`, `secp192k1`, `sect193r1/r2` (n ≈ 2¹⁹²) give the same 2⁻¹⁶⁷·⁰ as the earlier P-192 study. Three Koblitz curves are the exception: `sect163k1` has 12 colliding sums at B = 4096 (2 at B = 256 and 1024), `sect233k1` 20 and `sect239k1` 18. Out of 2²⁴ sums that is a relative loss of order 10⁻⁶ in the yield, which the tables show as no change at the printed precision. The cause was not investigated; the collisions are present already at B = 256 on `sect233k1` and `sect239k1` (12), so they come from a few early abscissae and not from the size of the base.

== Growth with the base

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right),
    table.hline(),
    [*curve*], [*B*], [*pair sums*], [*collisions*], [*exact yield*], [*B + 1 relations*], [*pair sums, s*], [*MC hits*],
    table.hline(stroke: 0.5pt),
    [`sect113r1`], [256], [2¹⁶], [0], [2⁻⁹⁵·⁰], [2¹⁰³·⁰], [0.06], [0],
    [`sect113r1`], [1024], [2²⁰], [0], [2⁻⁹¹·⁰], [2¹⁰¹·⁰], [0.18], [0],
    [`sect113r1`], [4096], [2²⁴], [0], [2⁻⁸⁷·⁰], [2⁹⁹·⁰], [2.72], [0],
    [`sect113r1`], [8192], [2²⁶], [0], [2⁻⁸⁵·⁰], [2⁹⁸·⁰], [16.78], [0],
    [`sect163r2`], [256], [2¹⁶], [0], [2⁻¹⁴⁵·⁰], [2¹⁵³·⁰], [0.02], [0],
    [`sect163r2`], [1024], [2²⁰], [0], [2⁻¹⁴¹·⁰], [2¹⁵¹·⁰], [0.26], [0],
    [`sect163r2`], [4096], [2²⁴], [0], [2⁻¹³⁷·⁰], [2¹⁴⁹·⁰], [9.28], [0],
    [`sect163r2`], [8192], [2²⁶], [0], [2⁻¹³⁵·⁰], [2¹⁴⁸·⁰], [33.70], [0],
    [`secp160r1`], [256], [2¹⁶], [0], [2⁻¹⁴³·⁰], [2¹⁵¹·⁰], [0.02], [0],
    [`secp160r1`], [1024], [2²⁰], [0], [2⁻¹³⁹·⁰], [2¹⁴⁹·⁰], [0.08], [0],
    [`secp160r1`], [4096], [2²⁴], [0], [2⁻¹³⁵·⁰], [2¹⁴⁷·⁰], [2.68], [0],
    [`secp160r1`], [8192], [2²⁶], [0], [2⁻¹³³·⁰], [2¹⁴⁶·⁰], [12.25], [0],
    table.hline(),
  ),
  caption: [Three curves from B = 256 to B = 8192. The yield grows as B², so each doubling of B gains 2 bits of yield and costs 2 bits of pair-sum work, while the trials for a full system fall by one bit. Source: `data/{sect113r1,sect163r2,secp160r1}_*.json`. Timings are one run each on a shared machine; the `sect163r2` B = 4096 value (9.3 s) is about 3× the neighbouring `sect163r1` and `sect163k1` rows and is not a property of the curve.],
)

== Method check on toy binary curves

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (right, right, right, right, right, right, right, right, right),
    table.hline(),
    [*n*], [*B*], [*collisions*], [*exact*], [*model*], [*MC rate*], [*hits*], [*MC / exact*], [*in 95% interval*],
    table.hline(stroke: 0.5pt),
    [4133], [16], [14], [0.11713], [0.12391], [0.11736], [352089], [1.0020], [yes],
    [4133], [32], [197], [0.40029], [0.49564], [0.39992], [1199759], [0.9991], [yes],
    [4133], [64], [2265], [0.88625], [1.9826], [0.8863], [2658891], [1.0000], [yes],
    [4133], [8], [0], [0.030978], [0.030978], [0.030947], [92841], [0.9990], [yes],
    [65579], [128], [3448], [0.39452], [0.49968], [0.39451], [1183524], [1.0000], [yes],
    [65579], [16], [0], [0.0078075], [0.0078075], [0.0078197], [23459], [1.0016], [yes],
    [65579], [256], [37123], [0.86654], [1.9987], [0.86629], [2598871], [0.9997], [yes],
    [65579], [32], [12], [0.030864], [0.03123], [0.030827], [92481], [0.9988], [yes],
    [65579], [64], [243], [0.11751], [0.12492], [0.11748], [352433], [0.9997], [yes],
    [261983], [1024], [917625], [0.99969], [8.0049], [0.99969], [2999056], [1.0000], [yes],
    [261983], [128], [976], [0.11763], [0.12508], [0.11748], [352443], [0.9988], [yes],
    [261983], [256], [13799], [0.39497], [0.50031], [0.39536], [1186070], [1.0010], [yes],
    [261983], [32], [0], [0.0078173], [0.0078173], [0.007788], [23364], [0.9962], [yes],
    [261983], [512], [149023], [0.86358], [2.0012], [0.86339], [2590181], [0.9998], [yes],
    [261983], [64], [33], [0.031017], [0.031269], [0.030804], [92411], [0.9931], [NO],
    table.hline(),
  ),
  caption: [All 15 toy binary runs, 3 000 000 trials each, cofactor 2. Yields are probabilities per trial; the model is 2B²/(n − 1). Source: `data/toyb*.json`.],
)

The exact count lies inside the 95% interval in 14 of 15 runs. The ratio of Monte Carlo rate to exact count is between 0.9931 and 1.0020. The run outside the interval (n = 261983, B = 64, ratio 0.9931) is 2.1 standard errors from the count (92 411 hits, relative standard error 0.33%); one miss in 15 runs is what 95% intervals predict, and the other 14 runs are inside. Where collisions are common the collision-free model is not close: at yield 0.8863 (n = 4133, B = 64) the model is 1.98. In every toy run the 32 planted decompositions were detected and recovered, and the first 8 Monte Carlo hits were recovered as real decompositions. This validates the counting, the cofactor map and the Monte Carlo on binary curves of cofactor 2; it does not exercise cofactor 4.

= Checks

- *Field.* The four-limb binary field agrees with the repository's `F2mElement` on random products, squares and inverses on `sect113r1` and `sect163r2`, and the carry-less kernels (hardware and portable) agree (unit test).
- *Group law.* 3·G on `sect113r1` equals the repository's `scalar_mul` (unit test). n·G = O and (n − 1)·G = −G for the published generators of `sect113r1`, `sect113r2`, `sect163r1`, `sect163r2` (binary) and for `secp112r1`, `secp128r1`, `secp160k1/r1/r2`, `secp192k1/r1` (prime) (unit tests).
- *Constants.* The database entries for `sect113r1`, `sect113r2` and `sect163r2` equal the hand-typed constructors in `src/binary_ecc` (a, b, Gx, Gy, n, polynomial, cofactor).
- *Base in the subgroup.* Every base point of `sect113r1`, `sect163r2` and of `secp112r2` (h = 4) satisfies n·P = O (unit tests).
- *Brute force.* On a 13-bit toy binary curve, every sum of two signed base points by the full group law reaches exactly twice as many distinct points as there are distinct sum abscissae (unit test).
- *Prime backend.* Montgomery arithmetic at 128, 160 and 192 bits agrees with big-integer arithmetic on random operands, and the generic P-192 curve gives the same yield, largest abscissa and distinct sums as the native P-192 field (unit tests).
- *Inspection.* `ic inspect --curve sect113r1` (and the other SECG names in `ic list`) passes the screens for the order, the modulus or reduction polynomial, nonsingularity, that the generator is on the curve and that [r]G = O.

= Limits

+ *Two-point sums only.* Three-point sums are not measured.
+ *Yield, not solving.* The cost of deciding a target (a Semaev S₃ solve) is not part of this measurement, and no linear algebra or descent runs. The "B + 1 relations" column is a count of trials, not a time.
+ *Cofactor convention.* The base is [h] of the small-x points, so it is in ⟨G⟩ but its abscissae are not small. The variant that keeps raw small-x points is not measured; the 1/h figure is an argument, not a measurement. The cofactor-4 path (`secp112r2`, `secp128r2`, `sect233k1`, `sect239k1`) is checked structurally (n·P = O, controls recovered) but has no toy statistical check.
+ *One base family.* Only the small-x base exists. The coverage-driven greedy selector of the binary Koblitz study chooses among Frobenius orbits and does not apply here.
+ *Size.* At most 8192 abscissae, bounded by the memory of the pair-sum set. 8192 was run on `sect113r1`, `sect163r2` and `secp160r1` only; the other curves stop at 4096.
+ *Coverage of the SECG list.* `secp224` and larger primes need more than three limbs and `sect283` and larger binary fields more than four; they are named and inspectable, but `ic prime` refuses them. `secp521r1` is not in `ic list` because the inspector reads at most 512-bit primes.
+ *Toy curves are not SECG curves.* They validate the method, not any property of a SECG curve.
+ *Timings are single runs on a shared 4-vCPU machine* and carry no more than their order of magnitude.

= Reproduction

```sh
cargo run --release --bin ic -- --json prime --curve sect113r1 --size 4096 --trials 1000000 --controls 32 --seed 7
cargo run --release --bin ic -- --json prime --curve sect163r2 --size 4096 --trials 1000000 --controls 32 --seed 7
cargo run --release --bin ic -- --json prime --curve secp160r1 --size 4096 --trials 1000000 --controls 32 --seed 7
cargo run --release --bin ic -- --json prime --curve toy-binary:17 --size 128 --trials 3000000 --controls 32 --seed 7
cargo run --release --bin ic -- inspect --curve sect163r2
python3 research/secg_small_x_20261014/graphs/plot_secg.py
```

`data/COMMANDS.txt` lists every run. Files: `data/*.json` (one `ic prime` report per run, 81 runs), `graphs/plot_secg.py`, `graphs/make_tables.py`, `graphs/fig1_secg.{svg,png}`.
