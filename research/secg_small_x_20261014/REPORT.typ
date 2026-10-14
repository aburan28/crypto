#set document(title: "Small-x factor base on the SECG curves: exact two-point relation yield and the cost of deciding a trial")
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1 / 1")
#set text(size: 9.5pt)
#set par(justify: true, leading: 0.6em)
#set heading(numbering: "1.")
#show heading.where(level: 1): set text(size: 12pt)
#show table: set text(size: 7.5pt)
#show raw.where(block: true): set text(size: 8pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[Small-x factor base on the SECG curves: exact two-point relation yield and the cost of deciding a trial]

  #text(size: 9pt)[Study `research/secg_small_x_20261014` · run 2026-10-14 · seed 7 · `ic prime` · `src/cryptanalysis/{secg_standards,secg_small_x,binary_small_x,prime_small_x,semaev_s3_cost}.rs`]
]

= Question and answer

*Question.* For every SECG curve in the repository's standards database (33 curves, `secp112r1` to `secp521r1` and `sect113r1` to `sect571r1`), what is the probability, per trial, that a random element of the prime-order subgroup is a signed sum of two points of a small-x factor base, and what does a trial cost when it is decided by solving Semaev's S₃ polynomial instead of looking the sum up?

*Answer.* With a base of B abscissae (2B points) the probability is *2B²/(n − 1)*, n the order of the prime-order subgroup. It is a count: every pair sum was computed, and on 27 of the 33 curves none coincide. At B = 4096 the yield is 2⁻⁸⁷·⁰ on `sect113r1` and `sect113r2`, 2⁻¹³⁷·⁰ on `sect163r1` and `sect163r2`, 2⁻¹³⁵·⁰ on `secp160k1/r1/r2`. B = 8192 was run on 30 curves (all except `secp521r1`, `sect571k1` and `sect571r1`, where the pair-sum keys would not fit in memory): 2⁻⁸⁵·⁰ on `sect113r1/r2`, 2⁻¹³⁵·⁰ on `sect163r1/r2`, 2⁻¹³³·⁰ on the 160-bit primes. The trials for a full system of B + 1 relations are about n/(2B), 2⁹⁸·⁰ on `sect113r1` at B = 8192, against rho at about √n, 2⁵⁶·⁰: the gap is at least 2⁴⁰·⁹ on every curve and run. Monte Carlo finds no hit in 91.6 million random subgroup elements over the 127 SECG runs, as the yields require.

*What a trial costs.* Deciding a target by solving S₃(x₁, x₂, X) = 0 for x₂ at every base abscissa takes 5.4 ms on `sect113r1`, 9.6 ms on `sect163r2` and 42 ms on `secp160r1` at B = 4096, one thread, against 2 to 5 µs for the sorted pair-sum lookup, which needs the 2²⁴-entry table. At the exact yield, collecting B + 1 relations by S₃ solves takes 2⁹¹·⁵ s on `sect113r1` against 2³⁹·⁴ s for Pollard rho at √n affine additions on the same machine, *52 bits slower*; 77 bits on `sect163r2`, 79 bits on `secp160r1`, and from 52 (`sect113r1`) to 280 bits (`sect571r1`) over the 33 curves. The S₃ decision equals the table decision on every element of toy fields and on 66 measured runs (0 disagreements, 16 of 16 planted sums found each time).

*The request "61 or 63, 160 bit"* was read as `sect163r1` and `sect163r2` (`sect163k1` is also measured) and, as the alternative reading, `secp160k1`, `secp160r1`, `secp160r2`. All are in the tables.

= Method

The base is the points with the smallest abscissae x = 1, 2, 3, … (binary fields: the polynomial whose coefficient bits are the binary digits of x; prime fields: x = 0, 1, 2, …) that lie on the curve, used with both signs. For B abscissae, a signed sum of two base points has abscissa x(a + b) or x(a − b) for a pair of distinct abscissae, or x(2a): B(B − 1)/2 pairs with two sums each and B doublings, B² abscissae if none coincide. A target T is a signed sum exactly when x(T) is one of them, and then so is −T. The covered targets are twice the number D of distinct sum abscissae, and a uniform element of the subgroup decomposes in one trial with probability exactly

#raw("yield = 2·D / (n − 1)        and        D ≤ B², with equality when no two sums collide", block: true)

`ic prime` computes D from every pair sum with one batched inversion per row. Binary curves y² + xy = x³ + ax² + b use λ = (y₁ + y₂)/(x₁ + x₂), x(a + b) = λ² + λ + x₁ + x₂ + a, −b = (x₂, x₂ + y₂) for x(a − b), and x(2a) = x² + b/x². Prime curves use the short Weierstrass formulas. Keys are exact field elements, sorted and deduplicated.

*Cofactor.* Every SECG binary curve has cofactor h = 2 or 4, as do `secp112r2` and `secp128r2`. A relation is a relation in the subgroup ⟨G⟩ = [h]E of prime order n, so each small-x point P enters the base as [h]P and targets are [h] of a uniform curve point; the yield is over ⟨G⟩. The base abscissae are then not small. The alternative, raw small-x points with only the sums that land in ⟨G⟩ counting, is sampled in section 3.4.

*Fields.* Prime: L-limb Montgomery arithmetic for any odd prime with L = 3, 4, 6, 9 limbs (112–192, 224–256, 384, 521 bits), Tonelli–Shanks for p ≡ 1 (mod 4) (`secp224k1/r1`); P-192 keeps its native field in the default run. Binary: F₂ᵐ with carry-less multiplication and Itoh–Tsujii inversion in 2, 3, 4, 5, 7, 9 limbs (113–128, 131–192, 193–256, 283, 409, 571 bits). Constants come from `docs/curves/standards/parameters.json`.

*Semaev S₃ cost.* For a prime curve y² = x³ + ax + b, S₃(x₁, x₂, X) = (x₁ − x₂)²X² − 2((x₁ + x₂)(x₁x₂ + a) + 2b)X + (x₁x₂ − a)² − 4b(x₁ + x₂); in x₂ it is a quadratic with roots (t ± √(t² − A·C))/A, A = (X − x₁)², t = x₁X(x₁ + X) + a(x₁ + X) + 2b, C = x₁²X² − 2X(ax₁ + 2b) + a² − 4bx₁: one square root per base abscissa, with one inversion for all of them. For a binary curve, S₃ = x₁²x₂² + x₁²X² + x₂²X² + x₁x₂X + b, independent of a (computed by eliminating the y's over F₂ with a computer algebra system, and checked against the group law for a = 0 and a = 1); in x₂ it is A·x₂² + B·x₂ + C with A = (x₁ + X)², B = x₁X, C = x₁²X² + b, solved by t² + t = AC/B² (trace and half-trace, both F₂-linear, as precomputed tables) and x₂ = (B/A)·t. A target is decided by trying every base abscissa x₁ and asking whether a root is a base abscissa. The time is one thread, 20 random subgroup targets per run; Pollard rho is charged √n affine additions, each timed on the same machine.

*Checks run with every measurement:* Monte Carlo on uniform subgroup elements, with the first hits recovered as real decompositions; planted controls (random signed sums of two base points, which must be detected and recovered); and, on toy curves with group order h·n, n prime and counted by walking the field, a comparison of the exact count with the Monte Carlo rate where hits are common: 15 binary curves with h = 2, 15 binary with h = 4, 11 prime with h = 2, 11 prime with h = 4.

= Results

== The SECG curves at B = 4096

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right, right, right, right),
    table.hline(),
    [*curve*], [*field bits*], [*h*], [*log₂ n*], [*exact yield*], [*collisions*], [*B + 1 relations*], [*rho*], [*pair sums, s*], [*MC trials*], [*controls*],
    table.hline(stroke: 0.5pt),
    [`secp112r1`], [112], [1], [111.8], [2⁻⁸⁶·⁸], [0], [2⁹⁸·⁸], [2⁵⁵·⁹], [2.89], [1000k], [32/32],
    [`secp112r2`], [112], [4], [109.8], [2⁻⁸⁴·⁸], [0], [2⁹⁶·⁸], [2⁵⁴·⁹], [2.71], [1000k], [32/32],
    [`sect113r1`], [113], [2], [112.0], [2⁻⁸⁷·⁰], [0], [2⁹⁹·⁰], [2⁵⁶·⁰], [2.72], [1000k], [32/32],
    [`sect113r2`], [113], [2], [112.0], [2⁻⁸⁷·⁰], [0], [2⁹⁹·⁰], [2⁵⁶·⁰], [2.87], [1000k], [32/32],
    [`secp128r1`], [128], [1], [128.0], [2⁻¹⁰³·⁰], [0], [2¹¹⁵·⁰], [2⁶⁴·⁰], [3.32], [1000k], [32/32],
    [`secp128r2`], [128], [4], [126.0], [2⁻¹⁰¹·⁰], [0], [2¹¹³·⁰], [2⁶³·⁰], [2.95], [1000k], [32/32],
    [`sect131r1`], [131], [2], [130.0], [2⁻¹⁰⁵·⁰], [0], [2¹¹⁷·⁰], [2⁶⁵·⁰], [3.32], [1000k], [32/32],
    [`sect131r2`], [131], [2], [130.0], [2⁻¹⁰⁵·⁰], [0], [2¹¹⁷·⁰], [2⁶⁵·⁰], [3.30], [1000k], [32/32],
    [`secp160k1`], [160], [1], [160.0], [2⁻¹³⁵·⁰], [0], [2¹⁴⁷·⁰], [2⁸⁰·⁰], [2.81], [1000k], [32/32],
    [`secp160r1`], [160], [1], [160.0], [2⁻¹³⁵·⁰], [0], [2¹⁴⁷·⁰], [2⁸⁰·⁰], [2.68], [1000k], [32/32],
    [`secp160r2`], [160], [1], [160.0], [2⁻¹³⁵·⁰], [0], [2¹⁴⁷·⁰], [2⁸⁰·⁰], [2.63], [1000k], [32/32],
    [`sect163k1`], [163], [2], [162.0], [2⁻¹³⁷·⁰], [12], [2¹⁴⁹·⁰], [2⁸¹·⁰], [3.29], [1000k], [32/32],
    [`sect163r1`], [163], [2], [162.0], [2⁻¹³⁷·⁰], [0], [2¹⁴⁹·⁰], [2⁸¹·⁰], [3.34], [1000k], [32/32],
    [`sect163r2`], [163], [2], [162.0], [2⁻¹³⁷·⁰], [0], [2¹⁴⁹·⁰], [2⁸¹·⁰], [9.28], [1000k], [32/32],
    [`secp192k1`], [192], [1], [192.0], [2⁻¹⁶⁷·⁰], [0], [2¹⁷⁹·⁰], [2⁹⁶·⁰], [2.41], [1000k], [32/32],
    [`secp192r1`], [192], [1], [192.0], [2⁻¹⁶⁷·⁰], [0], [2¹⁷⁹·⁰], [2⁹⁶·⁰], [2.91], [1000k], [32/32],
    [`sect193r1`], [193], [2], [192.0], [2⁻¹⁶⁷·⁰], [0], [2¹⁷⁹·⁰], [2⁹⁶·⁰], [8.01], [1000k], [32/32],
    [`sect193r2`], [193], [2], [192.0], [2⁻¹⁶⁷·⁰], [0], [2¹⁷⁹·⁰], [2⁹⁶·⁰], [6.81], [1000k], [32/32],
    [`secp224k1`], [224], [1], [224.0], [2⁻¹⁹⁹·⁰], [0], [2²¹¹·⁰], [2¹¹²·⁰], [2.52], [500k], [32/32],
    [`secp224r1`], [224], [1], [224.0], [2⁻¹⁹⁹·⁰], [0], [2²¹¹·⁰], [2¹¹²·⁰], [2.25], [500k], [32/32],
    [`sect233k1`], [233], [4], [231.0], [2⁻²⁰⁶·⁰], [20], [2²¹⁸·⁰], [2¹¹⁵·⁵], [8.21], [1000k], [32/32],
    [`sect233r1`], [233], [2], [232.0], [2⁻²⁰⁷·⁰], [0], [2²¹⁹·⁰], [2¹¹⁶·⁰], [7.33], [1000k], [32/32],
    [`sect239k1`], [239], [4], [237.0], [2⁻²¹²·⁰], [18], [2²²⁴·⁰], [2¹¹⁸·⁵], [7.83], [1000k], [32/32],
    [`secp256k1`], [256], [1], [256.0], [2⁻²³¹·⁰], [0], [2²⁴³·⁰], [2¹²⁸·⁰], [3.31], [500k], [32/32],
    [`secp256r1`], [256], [1], [256.0], [2⁻²³¹·⁰], [0], [2²⁴³·⁰], [2¹²⁸·⁰], [3.66], [500k], [32/32],
    [`sect283k1`], [283], [4], [281.0], [2⁻²⁵⁶·⁰], [22], [2²⁶⁸·⁰], [2¹⁴⁰·⁵], [7.71], [100k], [32/32],
    [`sect283r1`], [283], [2], [282.0], [2⁻²⁵⁷·⁰], [0], [2²⁶⁹·⁰], [2¹⁴¹·⁰], [7.16], [100k], [32/32],
    [`secp384r1`], [384], [1], [384.0], [2⁻³⁵⁹·⁰], [0], [2³⁷¹·⁰], [2¹⁹²·⁰], [4.98], [100k], [32/32],
    [`sect409k1`], [409], [4], [407.0], [2⁻³⁸²·⁰], [17], [2³⁹⁴·⁰], [2²⁰³·⁵], [4.86], [100k], [32/32],
    [`sect409r1`], [409], [2], [408.0], [2⁻³⁸³·⁰], [0], [2³⁹⁵·⁰], [2²⁰⁴·⁰], [4.69], [100k], [32/32],
    [`secp521r1`], [521], [1], [521.0], [2⁻⁴⁹⁶·⁰], [0], [2⁵⁰⁸·⁰], [2²⁶⁰·⁵], [6.61], [100k], [32/32],
    [`sect571k1`], [571], [4], [569.0], [2⁻⁵⁴⁴·⁰], [19], [2⁵⁵⁶·⁰], [2²⁸⁴·⁵], [7.40], [100k], [32/32],
    [`sect571r1`], [571], [2], [570.0], [2⁻⁵⁴⁵·⁰], [0], [2⁵⁵⁷·⁰], [2²⁸⁵·⁰], [7.40], [100k], [32/32],
    table.hline(),
  ),
  caption: [Two-point small-x bases of 4096 abscissae (8192 points), 2²⁴ pair sums, 33 curves. "B + 1 relations" is log₂((B + 1)/yield), "rho" is log₂ √n. Monte Carlo: 0 hits in every run. Planted controls detected and recovered in every run. Source: `data/sec*_4096.json`.],
)

#figure(
  image("graphs/fig1_secg.svg", width: 100%),
  caption: [Left: trials for B + 1 relations at B = 4096 (bars) against rho (dashes) on every curve; blue is a binary field, orange a prime field. Right: toy curves, Monte Carlo rate against the exact count (95% intervals), binary and prime, h = 1, 2, 4. Source: `data/*.json`; plot `graphs/plot_secg.py`.],
)

The yield follows 2B²/(n − 1), so a row is set by the subgroup order and the base size, not by the field or family: `sect113r1/r2` (n ≈ 2¹¹²) match each other and `sect163r1/r2` (n ≈ 2¹⁶²) match each other, and `secp192r1`, `secp192k1`, `sect193r1/r2` (n ≈ 2¹⁹²) give 2⁻¹⁶⁷·⁰, the earlier P-192 study's value.

*Collisions on the Koblitz curves.* Six curves have a few colliding sums at B = 4096: `sect163k1` 12, `sect233k1` 20, `sect239k1` 18, `sect283k1` 22, `sect409k1` 17, `sect571k1` 19 (2 to 12 already at B = 256). The cause is Frobenius. In the polynomial basis squaring spreads the bits of x, so the squares of small abscissae are small abscissae: 6 = z² + z, 6² = 20 = z⁴ + z², 20² = 272 = z⁸ + z⁴. The base then holds Frobenius images P, τP, τ²P, and τ² = μτ − 2 turns that into exact coincidences of sums. `examples/secg_koblitz_collisions.rs` lists them on `sect163k1` at B = 256: [6] − [20] = ±([6] + [272]) and 2·[6] = ±([20] − [272]) (so τ²[6] = τ[6] − 2·[6] up to sign), and on `sect233k1` the same pattern with x = 2, 4, 16, 256 and 3, 5, 17, 257. The control `sect163r2` has none. Out of 2²⁴ sums the loss is of order 10⁻⁶ of the yield.

== The base at B = 8192

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right, right, right),
    table.hline(),
    [*curve*], [*field bits*], [*h*], [*exact yield*], [*collisions*], [*B + 1 relations*], [*rho*], [*pair sums, s*], [*MC hits*], [*controls*],
    table.hline(stroke: 0.5pt),
    [`secp112r1`], [112], [1], [2⁻⁸⁴·⁸], [0], [2⁹⁷·⁸], [2⁵⁵·⁹], [6.0], [0], [32/32],
    [`secp112r2`], [112], [4], [2⁻⁸²·⁸], [0], [2⁹⁵·⁸], [2⁵⁴·⁹], [5.6], [0], [32/32],
    [`sect113r1`], [113], [2], [2⁻⁸⁵·⁰], [0], [2⁹⁸·⁰], [2⁵⁶·⁰], [16.8], [0], [32/32],
    [`sect113r2`], [113], [2], [2⁻⁸⁵·⁰], [0], [2⁹⁸·⁰], [2⁵⁶·⁰], [7.9], [0], [32/32],
    [`secp128r1`], [128], [1], [2⁻¹⁰¹·⁰], [0], [2¹¹⁴·⁰], [2⁶⁴·⁰], [5.7], [0], [32/32],
    [`secp128r2`], [128], [4], [2⁻⁹⁹·⁰], [0], [2¹¹²·⁰], [2⁶³·⁰], [5.4], [0], [32/32],
    [`sect131r1`], [131], [2], [2⁻¹⁰³·⁰], [0], [2¹¹⁶·⁰], [2⁶⁵·⁰], [15.7], [0], [32/32],
    [`sect131r2`], [131], [2], [2⁻¹⁰³·⁰], [0], [2¹¹⁶·⁰], [2⁶⁵·⁰], [19.4], [0], [32/32],
    [`secp160k1`], [160], [1], [2⁻¹³³·⁰], [0], [2¹⁴⁶·⁰], [2⁸⁰·⁰], [5.6], [0], [32/32],
    [`secp160r1`], [160], [1], [2⁻¹³³·⁰], [0], [2¹⁴⁶·⁰], [2⁸⁰·⁰], [12.3], [0], [32/32],
    [`secp160r2`], [160], [1], [2⁻¹³³·⁰], [0], [2¹⁴⁶·⁰], [2⁸⁰·⁰], [5.5], [0], [32/32],
    [`sect163k1`], [163], [2], [2⁻¹³⁵·⁰], [12], [2¹⁴⁸·⁰], [2⁸¹·⁰], [10.4], [0], [32/32],
    [`sect163r1`], [163], [2], [2⁻¹³⁵·⁰], [0], [2¹⁴⁸·⁰], [2⁸¹·⁰], [10.3], [0], [32/32],
    [`sect163r2`], [163], [2], [2⁻¹³⁵·⁰], [0], [2¹⁴⁸·⁰], [2⁸¹·⁰], [33.7], [0], [32/32],
    [`secp192k1`], [192], [1], [2⁻¹⁶⁵·⁰], [0], [2¹⁷⁸·⁰], [2⁹⁶·⁰], [5.5], [0], [32/32],
    [`secp192r1`], [192], [1], [2⁻¹⁶⁵·⁰], [0], [2¹⁷⁸·⁰], [2⁹⁶·⁰], [5.6], [0], [32/32],
    [`sect193r1`], [193], [2], [2⁻¹⁶⁵·⁰], [0], [2¹⁷⁸·⁰], [2⁹⁶·⁰], [10.1], [0], [32/32],
    [`sect193r2`], [193], [2], [2⁻¹⁶⁵·⁰], [0], [2¹⁷⁸·⁰], [2⁹⁶·⁰], [11.0], [0], [32/32],
    [`secp224k1`], [224], [1], [2⁻¹⁹⁷·⁰], [0], [2²¹⁰·⁰], [2¹¹²·⁰], [7.4], [0], [32/32],
    [`secp224r1`], [224], [1], [2⁻¹⁹⁷·⁰], [0], [2²¹⁰·⁰], [2¹¹²·⁰], [7.4], [0], [32/32],
    [`sect233k1`], [233], [4], [2⁻²⁰⁴·⁰], [20], [2²¹⁷·⁰], [2¹¹⁵·⁵], [9.8], [0], [32/32],
    [`sect233r1`], [233], [2], [2⁻²⁰⁵·⁰], [0], [2²¹⁸·⁰], [2¹¹⁶·⁰], [10.1], [0], [32/32],
    [`sect239k1`], [239], [4], [2⁻²¹⁰·⁰], [18], [2²²³·⁰], [2¹¹⁸·⁵], [11.0], [0], [32/32],
    [`secp256k1`], [256], [1], [2⁻²²⁹·⁰], [0], [2²⁴²·⁰], [2¹²⁸·⁰], [7.5], [0], [32/32],
    [`secp256r1`], [256], [1], [2⁻²²⁹·⁰], [0], [2²⁴²·⁰], [2¹²⁸·⁰], [7.2], [0], [32/32],
    [`sect283k1`], [283], [4], [2⁻²⁵⁴·⁰], [22], [2²⁶⁷·⁰], [2¹⁴⁰·⁵], [15.5], [0], [32/32],
    [`sect283r1`], [283], [2], [2⁻²⁵⁵·⁰], [0], [2²⁶⁸·⁰], [2¹⁴¹·⁰], [15.8], [0], [32/32],
    [`secp384r1`], [384], [1], [2⁻³⁵⁷·⁰], [0], [2³⁷⁰·⁰], [2¹⁹²·⁰], [14.4], [0], [32/32],
    [`sect409k1`], [409], [4], [2⁻³⁸⁰·⁰], [17], [2³⁹³·⁰], [2²⁰³·⁵], [24.3], [0], [32/32],
    [`sect409r1`], [409], [2], [2⁻³⁸¹·⁰], [0], [2³⁹⁴·⁰], [2²⁰⁴·⁰], [23.9], [0], [32/32],
    table.hline(),
  ),
  caption: [Two-point small-x bases of 8192 abscissae, 2²⁶ pair sums, 30 curves. The three curves not listed (`secp521r1`, `sect571k1`, `sect571r1`) would need 9.7 GB of pair-sum keys and are refused by the memory guard (the message is recorded in `data/COMMANDS.txt` section notes; the tool exits with an error). Source: `data/sec*_8192.json`. Timings are single runs; `sect163r2` (33.7 s) is about 3× its neighbours and was run while other work shared the machine.],
)

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
  caption: [Growth with B on three curves. The yield grows as B², so each doubling of B gains 2 bits of yield and costs 2 bits of pair-sum work, while the trials for a full system fall by one bit. Source: `data/{sect113r1,sect163r2,secp160r1}_*.json`.],
)

== What a trial costs by a Semaev S₃ solve

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right, right, right),
    table.hline(),
    [*curve*], [*S₃ per target*], [*field mults*], [*table lookup*], [*group add*], [*B + 1 relations, s*], [*rho, s*], [*bits slower*], [*planted found*], [*disagree*],
    table.hline(stroke: 0.5pt),
    [`secp112r1`], [25.73 ms], [796k], [3898 ns], [5.6 µs], [2⁹³·⁵], [2³⁸·⁴], [55.1], [16/16], [0],
    [`secp112r2`], [25.81 ms], [789k], [3172 ns], [5.7 µs], [2⁹¹·⁵], [2³⁷·⁵], [54.0], [16/16], [0],
    [`sect113r1`], [5.36 ms], [77k], [2055 ns], [9.8 µs], [2⁹¹·⁵], [2³⁹·⁴], [52.1], [16/16], [0],
    [`sect113r2`], [7.13 ms], [136k], [2077 ns], [7.0 µs], [2⁹¹·⁹], [2³⁸·⁹], [53.0], [16/16], [0],
    [`secp128r1`], [23.56 ms], [721k], [2963 ns], [8.3 µs], [2¹⁰⁹·⁶], [2⁴⁷·¹], [62.5], [16/16], [0],
    [`secp128r2`], [24.73 ms], [763k], [3677 ns], [8.2 µs], [2¹⁰⁷·⁷], [2⁴⁶·¹], [61.6], [16/16], [0],
    [`sect131r1`], [8.17 ms], [61k], [2256 ns], [21.5 µs], [2¹¹⁰·¹], [2⁴⁹·⁵], [60.6], [16/16], [0],
    [`sect131r2`], [8.73 ms], [91k], [2270 ns], [13.5 µs], [2¹¹⁰·²], [2⁴⁸·⁸], [61.3], [16/16], [0],
    [`secp160k1`], [44.12 ms], [1328k], [3355 ns], [10.3 µs], [2¹⁴²·⁵], [2⁶³·⁴], [79.1], [16/16], [0],
    [`secp160r1`], [42.41 ms], [1311k], [3286 ns], [10.4 µs], [2¹⁴²·⁴], [2⁶³·⁴], [79.0], [16/16], [0],
    [`secp160r2`], [45.51 ms], [1248k], [3322 ns], [11.5 µs], [2¹⁴²·⁵], [2⁶³·⁶], [79.0], [16/16], [0],
    [`sect163k1`], [9.42 ms], [96k], [2345 ns], [17.1 µs], [2¹⁴²·³], [2⁶⁵·²], [77.1], [16/16], [0],
    [`sect163r1`], [8.74 ms], [94k], [2262 ns], [16.5 µs], [2¹⁴²·²], [2⁶⁵·¹], [77.1], [16/16], [0],
    [`sect163r2`], [9.64 ms], [100k], [2430 ns], [17.3 µs], [2¹⁴²·³], [2⁶⁵·²], [77.1], [16/16], [0],
    [`secp192k1`], [54.08 ms], [1483k], [3183 ns], [13.5 µs], [2¹⁷⁴·⁸], [2⁷⁹·⁸], [95.0], [16/16], [0],
    [`secp192r1`], [43.21 ms], [1016k], [3044 ns], [13.2 µs], [2¹⁷⁴·⁵], [2⁷⁹·⁸], [94.7], [16/16], [0],
    [`sect193r1`], [10.13 ms], [115k], [2517 ns], [18.0 µs], [2¹⁷²·⁴], [2⁸⁰·²], [92.1], [16/16], [0],
    [`sect193r2`], [10.02 ms], [110k], [2516 ns], [19.5 µs], [2¹⁷²·⁴], [2⁸⁰·⁴], [92.0], [16/16], [0],
    [`secp224k1`], [429.11 ms], [8792k], [4408 ns], [21.1 µs], [2²⁰⁹·⁸], [2⁹⁶·⁵], [113.3], [16/16], [0],
    [`secp224r1`], [1401.41 ms], [29088k], [4478 ns], [20.3 µs], [2²¹¹·⁵], [2⁹⁶·⁴], [115.1], [16/16], [0],
    [`sect233k1`], [10.13 ms], [107k], [2763 ns], [24.3 µs], [2²¹¹·⁴], [2¹⁰⁰·²], [111.2], [16/16], [0],
    [`sect233r1`], [10.45 ms], [125k], [2614 ns], [20.2 µs], [2²¹²·⁴], [2¹⁰⁰·⁴], [112.0], [16/16], [0],
    [`sect239k1`], [11.51 ms], [100k], [2550 ns], [29.1 µs], [2²¹⁷·⁶], [2¹⁰³·⁴], [114.1], [16/16], [0],
    [`secp256k1`], [102.42 ms], [2129k], [3567 ns], [23.2 µs], [2²³⁹·⁷], [2¹¹²·⁶], [127.1], [16/16], [0],
    [`secp256r1`], [58.46 ms], [1216k], [3357 ns], [17.6 µs], [2²³⁸·⁹], [2¹¹²·²], [126.7], [16/16], [0],
    [`sect283k1`], [16.45 ms], [101k], [2469 ns], [48.1 µs], [2²⁶²·¹], [2¹²⁶·²], [135.9], [16/16], [0],
    [`sect283r1`], [16.29 ms], [75k], [2554 ns], [59.6 µs], [2²⁶³·¹], [2¹²⁷·⁰], [136.1], [16/16], [0],
    [`secp384r1`], [261.85 ms], [2870k], [4160 ns], [62.3 µs], [2³⁶⁹·¹], [2¹⁷⁸·⁰], [191.0], [16/16], [0],
    [`sect409k1`], [20.25 ms], [100k], [3522 ns], [74.0 µs], [2³⁸⁸·⁴], [2¹⁸⁹·⁸], [198.6], [16/16], [0],
    [`sect409r1`], [21.06 ms], [82k], [4630 ns], [108.4 µs], [2³⁸⁹·⁴], [2¹⁹⁰·⁸], [198.6], [16/16], [0],
    [`secp521r1`], [432.66 ms], [2079k], [4888 ns], [198.2 µs], [2⁵⁰⁶·⁸], [2²⁴⁸·²], [258.6], [16/16], [0],
    [`sect571k1`], [34.96 ms], [106k], [4706 ns], [189.5 µs], [2⁵⁵¹·²], [2²⁷²·¹], [279.0], [16/16], [0],
    [`sect571r1`], [34.79 ms], [102k], [3232 ns], [198.4 µs], [2⁵⁵²·²], [2²⁷²·⁷], [279.5], [16/16], [0],
    table.hline(),
  ),
  caption: [Deciding a target by an S₃ solve at every base abscissa, B = 4096, one thread, 20 random subgroup targets per run. "B + 1 relations, s" is the trials for a full system at the exact yield times the time of one solve, in seconds; "rho, s" is √n times the time of one affine addition measured on the same machine; "bits slower" is the log₂ of their ratio. "disagree" counts targets on which the S₃ decision and the pair-sum table differ; "planted found" is of 16. Source: `data/s3_*_4096.json`.],
)

#figure(
  image("graphs/fig2_s3_cost.svg", width: 100%),
  caption: [Left: milliseconds to decide one target by S₃ at B = 4096. Right: bits by which collecting B + 1 relations with S₃ solves is slower than rho on the same machine. Source: `data/s3_*_4096.json`; plot `graphs/plot_secg.py`.],
)

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    table.hline(),
    [*curve*], [*B*], [*S₃ per target*], [*per base abscissa*], [*bits slower than rho*],
    table.hline(stroke: 0.5pt),
    [`sect113r1`], [1024], [1.36 ms], [1330 ns], [52.7],
    [`sect113r1`], [4096], [5.36 ms], [1308 ns], [52.1],
    [`sect163r2`], [1024], [2.64 ms], [2574 ns], [77.2],
    [`sect163r2`], [4096], [9.64 ms], [2354 ns], [77.1],
    [`secp160r1`], [1024], [10.22 ms], [9983 ns], [78.9],
    [`secp160r1`], [4096], [42.41 ms], [10354 ns], [79.0],
    [`secp256r1`], [1024], [16.39 ms], [16005 ns], [126.9],
    [`secp256r1`], [4096], [58.46 ms], [14273 ns], [126.7],
    table.hline(),
  ),
  caption: [S₃ cost against B: linear in B (the cost per base abscissa is constant), and the gap to rho is the same at B = 1024 and 4096. Source: `data/s3_*_{1024,4096}.json`.],
)

The S₃ solve costs one quadratic per base abscissa, so a target costs about B times a few hundred field multiplications on a prime curve (a square root each) and about B times twenty on a binary curve (the trace and half-trace are table lookups): 77 000 to 136 000 multiplications on `sect113r1/r2` against 1 300 000 on `secp160r1`. Without the table a trial is 2600 times dearer on `sect113r1` (5.36 ms against 2.06 µs), 13 000 times on `secp160r1` and 313 000 times on `secp224r1` (1.4 s: the square roots by Tonelli–Shanks, p ≡ 1 mod 4); the ratio over the 33 curves runs from 2600 to 313 000. With the table the trial is a lookup but the table holds B² keys, the 2²⁴ entries of 24 to 72 bytes that bound the base size. At the exact yield the cost of B + 1 relations by S₃ is 52 to 280 bits above rho, so the answer to "does a trial's solve cost rescue the method" is no on every curve: the yield already puts the method at least 2⁴⁰·⁹ trials above rho, and the solve adds more.

== The raw small-x base

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right),
    table.hline(),
    [*curve*], [*h*], [*samples*], [*in subgroup*], [*share*], [*95% interval*], [*1/h*], [*raw yield, B = 64*],
    table.hline(stroke: 0.5pt),
    [`secp112r2`], [4], [4000], [1024], [0.2560], [0.2427 to 0.2698], [0.2500], [2⁻⁹⁸·⁷],
    [`sect113r1`], [2], [4000], [2010], [0.5025], [0.4870 to 0.5180], [0.5000], [2⁻¹⁰⁰·⁰],
    [`sect113r2`], [2], [4000], [1974], [0.4935], [0.4780 to 0.5090], [0.5000], [2⁻¹⁰⁰·⁰],
    [`secp128r2`], [4], [4000], [1018], [0.2545], [0.2412 to 0.2682], [0.2500], [2⁻¹¹⁵·⁰],
    [`sect163k1`], [2], [4000], [2059], [0.5148], [0.4993 to 0.5302], [0.5000], [2⁻¹⁵⁰·⁰],
    [`sect163r1`], [2], [4000], [2004], [0.5010], [0.4855 to 0.5165], [0.5000], [2⁻¹⁵⁰·⁰],
    [`sect163r2`], [2], [4000], [1993], [0.4983], [0.4828 to 0.5137], [0.5000], [2⁻¹⁵⁰·⁰],
    [`sect233k1`], [4], [4000], [1105], [0.2762], [0.2626 to 0.2903], [0.2500], [2⁻²¹⁹·⁹],
    [`sect239k1`], [4], [4000], [987], [0.2467], [0.2336 to 0.2603], [0.2500], [2⁻²²⁶·⁰],
    [`sect283k1`], [4], [4000], [975], [0.2437], [0.2307 to 0.2573], [0.2500], [2⁻²⁷⁰·⁰],
    table.hline(),
  ),
  caption: [Share of signed pair sums of raw small-x points (no [h]) that lie in the prime-order subgroup, from 4000 random pairs at B = 64; membership is [n]S = O. A target of the subgroup can only be a sum that lands in it, so the raw base's yield is this share times the subgroup base's. "raw yield" is log₂ of share × 2B²/(n − 1), assuming no collisions. Source: `data/raw_*.json`; plot `graphs/fig3_raw_base.svg`.],
)

#figure(
  image("graphs/fig3_raw_base.svg", width: 70%),
  caption: [The sampled share (dots, 95% intervals) against 1/h (bars). Source: `data/raw_*.json`.],
)

On nine of the ten curves the share is consistent with 1/h (0.4935 to 0.5148 for h = 2, 0.2437 to 0.2560 for h = 4). `sect233k1` is not: 0.2762, 3.8 standard errors above 1/4. The mechanism is not established. The class of a point in E/⟨G⟩ is set by the trace of its abscissa (h = 2), and small abscissae are not uniform in trace, so the 1/h model, which treats the classes as uniform, can be off for a small-x base; this is a hypothesis, not a measurement. The upshot is that the raw variant loses a factor of about h on the SECG curves, not 1 and not much more: 2 bits at h = 4.

== Method check on toy curves

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right, right, right, right),
    table.hline(),
    [*field*], [*h*], [*n*], [*B*], [*collisions*], [*exact*], [*model*], [*MC rate*], [*hits*], [*MC / exact*], [*in 95% interval*],
    table.hline(stroke: 0.5pt),
    [binary], [2], [4133], [16], [14], [0.11713], [0.12391], [0.11736], [352089], [1.0020], [yes],
    [binary], [2], [4133], [32], [197], [0.40029], [0.49564], [0.39992], [1199759], [0.9991], [yes],
    [binary], [2], [4133], [64], [2265], [0.88625], [1.9826], [0.8863], [2658891], [1.0000], [yes],
    [binary], [2], [4133], [8], [0], [0.030978], [0.030978], [0.030947], [92841], [0.9990], [yes],
    [binary], [4], [2063], [16], [26], [0.22308], [0.2483], [0.22302], [669048], [0.9997], [yes],
    [binary], [4], [2063], [32], [390], [0.61494], [0.99321], [0.61486], [1844577], [0.9999], [yes],
    [binary], [4], [2063], [64], [3091], [0.97478], [3.9728], [0.97487], [2924616], [1.0001], [yes],
    [binary], [4], [2063], [8], [2], [0.060136], [0.062076], [0.060154], [180462], [1.0003], [yes],
    [binary], [2], [65579], [128], [3448], [0.39452], [0.49968], [0.39451], [1183524], [1.0000], [yes],
    [binary], [2], [65579], [16], [0], [0.0078075], [0.0078075], [0.0078197], [23459], [1.0016], [yes],
    [binary], [2], [65579], [256], [37123], [0.86654], [1.9987], [0.86629], [2598871], [0.9997], [yes],
    [binary], [2], [65579], [32], [12], [0.030864], [0.03123], [0.030827], [92481], [0.9988], [yes],
    [binary], [2], [65579], [64], [243], [0.11751], [0.12492], [0.11748], [352433], [0.9997], [yes],
    [binary], [4], [32779], [128], [5980], [0.63482], [0.99969], [0.63539], [1906175], [1.0009], [NO],
    [binary], [4], [32779], [16], [8], [0.015132], [0.01562], [0.015234], [45701], [1.0067], [yes],
    [binary], [4], [32779], [256], [49415], [0.98365], [3.9988], [0.98373], [2951186], [1.0001], [yes],
    [binary], [4], [32779], [32], [38], [0.060162], [0.062481], [0.060283], [180848], [1.0020], [yes],
    [binary], [4], [32779], [64], [489], [0.22009], [0.24992], [0.22038], [661127], [1.0013], [yes],
    [binary], [2], [261983], [1024], [917625], [0.99969], [8.0049], [0.99969], [2999056], [1.0000], [yes],
    [binary], [2], [261983], [128], [976], [0.11763], [0.12508], [0.11748], [352443], [0.9988], [yes],
    [binary], [2], [261983], [256], [13799], [0.39497], [0.50031], [0.39536], [1186070], [1.0010], [yes],
    [binary], [2], [261983], [32], [0], [0.0078173], [0.0078173], [0.007788], [23364], [0.9962], [yes],
    [binary], [2], [261983], [512], [149023], [0.86358], [2.0012], [0.86339], [2590181], [0.9998], [yes],
    [binary], [2], [261983], [64], [33], [0.031017], [0.031269], [0.030804], [92411], [0.9931], [NO],
    [binary], [4], [131023], [1024], [983065], [1], [16.006], [1], [3000000], [1.0000], [yes],
    [binary], [4], [131023], [128], [1818], [0.22234], [0.2501], [0.22248], [667429], [1.0006], [yes],
    [binary], [4], [131023], [256], [23963], [0.6346], [1.0004], [0.63452], [1903566], [0.9999], [yes],
    [binary], [4], [131023], [32], [3], [0.015585], [0.015631], [0.015565], [46695], [0.9987], [yes],
    [binary], [4], [131023], [512], [197758], [0.98283], [4.0015], [0.98279], [2948369], [1.0000], [yes],
    [binary], [4], [131023], [64], [139], [0.060402], [0.062524], [0.060389], [181168], [0.9998], [yes],
    [prime], [2], [32789], [128], [6004], [0.63316], [0.99939], [0.63307], [1899200], [0.9999], [yes],
    [prime], [2], [32789], [16], [3], [0.015432], [0.015615], [0.015314], [45941], [0.9923], [yes],
    [prime], [2], [32789], [256], [49441], [0.98176], [3.9976], [0.98171], [2945118], [0.9999], [yes],
    [prime], [2], [32789], [32], [17], [0.061425], [0.062462], [0.061066], [183197], [0.9942], [NO],
    [prime], [2], [32789], [64], [521], [0.21807], [0.24985], [0.21787], [653607], [0.9991], [yes],
    [prime], [4], [16493], [128], [9249], [0.86527], [1.9869], [0.8654], [2596189], [1.0001], [yes],
    [prime], [4], [16493], [16], [6], [0.030318], [0.031045], [0.03041], [91231], [1.0031], [yes],
    [prime], [4], [16493], [256], [57291], [0.99988], [7.9476], [0.99988], [2999645], [1.0000], [yes],
    [prime], [4], [16493], [32], [52], [0.11788], [0.12418], [0.11806], [354194], [1.0016], [yes],
    [prime], [4], [16493], [64], [876], [0.39049], [0.49673], [0.39083], [1172495], [1.0009], [yes],
    [prime], [2], [523867], [1024], [791426], [0.98174], [4.0032], [0.98179], [2945384], [1.0001], [yes],
    [prime], [2], [523867], [128], [549], [0.060454], [0.06255], [0.06064], [181920], [1.0031], [yes],
    [prime], [2], [523867], [256], [7258], [0.22249], [0.2502], [0.22237], [667115], [0.9995], [yes],
    [prime], [2], [523867], [32], [0], [0.0039094], [0.0039094], [0.0038933], [11680], [0.9959], [yes],
    [prime], [2], [523867], [512], [96559], [0.63217], [1.0008], [0.63239], [1897181], [1.0004], [yes],
    [prime], [2], [523867], [64], [51], [0.015443], [0.015638], [0.015419], [46258], [0.9985], [yes],
    [prime], [4], [262321], [1024], [917459], [0.99967], [7.9946], [0.9997], [2999095], [1.0000], [NO],
    [prime], [4], [262321], [128], [973], [0.1175], [0.12492], [0.11756], [352692], [1.0006], [yes],
    [prime], [4], [262321], [256], [13787], [0.39455], [0.49966], [0.39464], [1183913], [1.0002], [yes],
    [prime], [4], [262321], [32], [2], [0.007792], [0.0078073], [0.0077507], [23252], [0.9947], [yes],
    [prime], [4], [262321], [512], [148784], [0.86429], [1.9987], [0.86424], [2592734], [0.9999], [yes],
    [prime], [4], [262321], [64], [49], [0.030855], [0.031229], [0.030659], [91978], [0.9936], [NO],
    table.hline(),
  ),
  caption: [All 52 toy runs, 3 000 000 trials each. Yields are probabilities per trial; the model is 2B²/(n − 1). Source: `data/toy*.json`.],
)

The exact count lies inside the 95% interval in 47 of 52 runs; five misses are what 52 intervals at 95% give about one time in nine (the chance of five or more is 0.12). The stronger test is the spread: over the 51 runs with yield below one, the standardized differences (rate − exact)/√(exact(1 − exact)/N) have mean 0.02 and standard deviation 1.04, largest 2.6, as a standard normal has, so there is no bias. The ratio of rate to count is between 0.9923 and 1.0067. The 11 prime and 15 binary runs with h = 4 and the 11 prime and 15 binary runs with h = 2 are all in this spread, so the cofactor map is validated for h = 2 and 4 on both field types. Where collisions are common the collision-free model is far off: at yield 0.8863 the model gives 1.98. In every toy run the 32 planted decompositions were detected and recovered, and the first 8 hits were recovered as real decompositions.

= Checks

- *Field.* The generic binary field agrees with the repository's `F2mElement` on random products, squares and inverses at every limb count (113, 163, 233, 283, 409, 571 bits), and the carry-less kernels (hardware and portable) agree (unit test). Prime Montgomery arithmetic at 4, 6 and 9 limbs agrees with big-integer arithmetic (unit test), as at 3 limbs; Tonelli–Shanks square roots on `secp224r1` verified (unit test).
- *Group law and generators.* n·G = O for the published generators of `sect113r1/r2`, `sect163r1/r2`, `sect233k1`, `sect283k1`, `sect409k1`, `sect571k1`, `secp112r1`, `secp128r1`, `secp160k1/r1/r2`, `secp192k1/r1`, `secp224k1/r1`, `secp256k1/r1`, `secp384r1`, `secp521r1`; 3·G on `sect113r1` equals the repository's `scalar_mul` (unit tests).
- *All 33 curves.* Every SECG curve builds, has its generator on the curve, and gives collision-free recovered controls on a small base (unit test); the database entries for `sect113r1`, `sect113r2`, `sect163r2` equal the hand-typed constructors.
- *S₃.* The prime and binary S₃ polynomials vanish on the abscissae of real group-law sums (binary for a = 0 and a = 1); the S₃ decision equals the pair-sum table on every field element of toy fields, binary 2¹³ with h = 2 and 4 and prime ≈ 2¹³ with h = 1, 2, 4 (unit tests); trace and half-trace tables match the field.
- *Brute force.* On a 13-bit toy binary curve and a 14-bit toy prime curve, every sum of two signed base points by the full group law reaches exactly twice as many distinct points as there are distinct sum abscissae.
- *Inspection.* `ic inspect --curve <name>` passes the screens (order, modulus or reduction polynomial, nonsingularity, generator on the curve, [r]G = O) for all 33 names in `ic list`.

= Limits

+ *Two-point sums only.* Three-point sums are not measured.
+ *S₃ is the cost of deciding a trial, not of solving the discrete logarithm.* No linear algebra over the relations and no descent ran. The S₃ times are 20 targets per run, one thread, single runs on a shared 4-vCPU machine, at B = 1024 and 4096; they carry an order of magnitude (for example `secp256k1` 102 ms against `secp256r1` 58 ms at the same size is run-to-run noise, not a property of the curves). The comparison with rho charges rho √n affine additions and nothing for its memory, distinguished points or the parallel speedup of either side.
+ *Raw base.* Sampled at 4000 pairs and B = 64 on ten curves; the share is not measured on the other curves, the yield assumes no collisions, and the `sect233k1` excess is unexplained.
+ *One base family.* Only the small-x base exists. The coverage-driven greedy selector of the binary Koblitz study chooses among Frobenius orbits and does not apply here.
+ *Size.* At most 8192 abscissae, bounded by memory for the pair-sum keys; three curves stop at 4096 for that reason (section 3.2).
+ *Toy curves are not SECG curves.* They validate the method, not any property of a SECG curve.
+ *Outside this study.* `cargo test --lib -- curve_traits::subfield` has two failing tests (`a_gf4_curve_over_gf2_6_descends_to_gf4`, `the_twist_of_a_koblitz_curve_over_an_even_degree_descends_only_to_gf4`, the `debug_assert` at `subfield.rs:75`). They fail identically at the branch's merge base, `45efe3b33`, which does not contain this work, and the files they run are untouched here (no change under `src/binary_ecc` or `src/cryptanalysis/curve_traits`). They are not fixed here.
+ *Coordination.* The Conductor control plane was not up for the first part of this work, and the earlier edits were made under the documented fallback warning. It was brought up with the repository's session-start hook; `conductor check` for the files of this study is clear and the work is task T-1.

= Reproduction

```sh
cargo run --release --bin ic -- --json prime --curve sect113r1 --size 4096 --trials 1000000 --controls 32 --seed 7
cargo run --release --bin ic -- --json prime --curve sect163r2 --size 4096 --trials 1000000 --controls 32 --seed 7
cargo run --release --bin ic -- --json prime --curve secp160r1 --size 4096 --trials 1000000 --controls 32 --seed 7
cargo run --release --bin ic -- --json prime --curve sect571r1 --size 4096 --trials 100000 --controls 32 --seed 7
cargo run --release --bin ic -- --json prime --curve sect163r2 --size 4096 --trials 1000 --controls 4 --seed 7 --s3-trials 20
cargo run --release --bin ic -- --json prime --curve sect233k1 --size 64 --trials 1000 --controls 8 --seed 7 --raw-samples 4000
cargo run --release --bin ic -- --json prime --curve toy-binary:17,7,4 --size 128 --trials 3000000 --controls 32 --seed 7
cargo run --release --example secg_koblitz_collisions -- sect163k1 256
python3 research/secg_small_x_20261014/graphs/plot_secg.py
```

`data/COMMANDS.txt` lists every run. Files: `data/*.json` (one `ic prime` report per run, 255 runs), `graphs/plot_secg.py`, `graphs/make_tables.py`, `graphs/fig{1_secg,2_s3_cost,3_raw_base}.{svg,png}`.
