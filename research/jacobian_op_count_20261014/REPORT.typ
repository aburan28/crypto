#set page(paper: "a4", margin: (x: 1.8cm, y: 1.5cm))
#set text(size: 9.5pt)
#set par(justify: true)
#show heading: set block(above: 1em, below: 0.5em)

= Counted Cantor against the cost model's 10 g² multiplications

_2026-10-14. Frozen protocol `PROTOCOL.md` (sha256 `c27e54a3c7bbbb33…`). Native instrument `examples/jacobian_op_count.rs`. All 10400 counted results agree with the repository's group law (0 mismatches)._

== Result (exact computation)

Mean counts for one Cantor addition of two random degree-`g` divisors over `F_1009`, 200 pairs per genus, with the total weighting an inversion as 8 multiplications. `ratio` is `log2` of the total over the model's `10 g²`; `local e` is the exponent between this genus and the previous measured one.

#table(
  columns: (auto, auto, auto, auto, auto, auto, auto),
  align: (right, right, right, right, right, right, right),
  stroke: 0.4pt + luma(180),
  inset: 3pt,
  [*g*], [*mults*], [*inversions*], [*total*], [*10 g²*], [*ratio (bits)*], [*local e*],
  [2], [124.6], [5.00], [164.7], [40], [2.041], [—],
  [4], [505.5], [8.01], [569.6], [160], [1.832], [1.807],
  [6], [1266.4], [11.06], [1354.9], [360], [1.912], [2.022],
  [8], [2508.6], [14.04], [2620.9], [640], [2.034], [2.294],
  [10], [4448.1], [17.18], [4585.6], [1000], [2.197], [2.507],
  [14], [10657.0], [23.32], [10843.5], [1960], [2.468], [2.624],
  [16], [15132.6], [26.36], [15343.5], [2560], [2.583], [2.600],
  [24], [46026.4], [38.63], [46335.5], [5760], [3.008], [2.766],
  [32], [101876.9], [50.70], [102282.5], [10240], [3.320], [2.752],
)

*Verdict (fixed wording).* The model is below textbook Cantor by 2.041 bits at genus 2 and 3.320 bits at genus 32 (addition, `p = 1009`; doubling 3.210 bits and `p = 65521` 3.231 bits at genus 32). The ratio is outside the one-bit band at every measured genus. The fitted exponent is 2.347 (95% bootstrap interval 2.346 to 2.349) for addition, 2.349 (2.348 to 2.350) for doubling and 2.319 (2.319 to 2.320) for addition at `p = 65521`: the exponent 2 is outside each interval. The intervals cover divisor sampling only; the local exponents rise from about 1.8 to about 2.75, so a single power law is an approximation.

#figure(image("ratio_by_genus.svg", width: 90%), caption: [Counted cost over `10 g²` by genus, addition and doubling, two primes. The grey band is the one-bit decision rule; dotted lines mark genus 6, 10 and 14, where most extension routes sit. Source: `summary.csv`.])

== Could the gap move a verdict? (exact count on the committed JSON)

Of the 848 `formula_advantage` verdicts, 14 lie within 1 bit of the 8-bit margin, 23 within 2 bits and 32 within 3 bits; 50 more are `formula_marginal`. If the measured cost replaced `10 g²` and the *whole* route cost scaled with it, 24 of the 848 would drop below the margin. That is an upper bound on fragility: routes priced by a different solver do not scale entirely with the group operation.

#figure(image("verdict_fragility.svg", width: 90%), caption: [Advantage verdicts close to the margin, and the conservative count that the measured cost would push below it. Source: `summary.csv`.])

== Claims and scope

- *Exact computation:* the counts per genus, the fit, and the verdict counts above.
- *Not shown:* that the model is wrong. Textbook Cantor with extended Euclid is not the cheapest algorithm (explicit formulas exist for `g = 2, 3`), so the ratio is a statement about this algorithm, not about the best one. The model's unit may also count an optimised operation.
- *Not shown:* anything about `F_p` of cryptographic size; the primes are `1009` and `65521`, and the two primes differ by at most 0.09 bits in the ratio at any measured genus, so no dependence on `p` is visible here.
- *Dominant caveat:* 846 of the 848 advantage verdicts are extension-field `cover_descent` routes marked "proposed", whose genus-`k` cover over `GF(p)` is not established (the route's own note). Whether those covers exist matters more to those verdicts than any constant in this report; see the odd-characteristic descent study.
