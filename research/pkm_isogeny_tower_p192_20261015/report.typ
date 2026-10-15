#set page(paper: "a4", margin: 2cm, numbering: "1")
#set text(size: 10pt)
#set par(justify: true)
#show heading: set block(above: 1.2em, below: 0.6em)

#align(center)[
  #text(size: 15pt, weight: "bold")[The PKM isogeny-tower auxiliary curve at the NIST P-192 prime]\
  #text(size: 11pt)[A 2-Sylow halving search: cost per draw, success probability, and a supersingular shortcut]\
  #v(0.3em)
  2026-10-15 · branch `claude/isogeny-tower-p192` · protocol `fc6317e51` (committed before data) · run `20261015T2010Z`
]

= Summary

- *Built.* `two_sylow_search.rs` and the `pkm-aux-search` binary. Each draw is a curve with full rational 2-torsion over $p = 2^192 - 2^64 - 1$. Its 2-Sylow is $ZZ\/2 times ZZ\/2^b$, and x-only point halving computes $b$ with no point counting. A point of order $2^t$ exists iff $b >= t$. Found points are certified over `BigUint` and accepted by the unmodified checker `backlog_check_aux2t`.
- *Per draw (measured).* On the $X_1(2,8)$ family (`z2xz8`, $b >= 3$ forced), the median is *10.06 µs*, against 30.83 µs for the baseline (0.326×). On `z2xz2` it is 13.89 µs (0.561× of 24.74). All four pre-registered speed claims pass.
- *Success probability (measured).* $2^28$ `z2xz8` draws and $2^25$ `z2xz2` draws match $P(b >= t) = 2^(3-t)$ and $3 dot 2^(-t)$ inside 99.9% Poisson bands at every tested $t$. The deepest ordinary curve has $b = 38$.
- *$t = 48$ (derived).* $2^45$ draws × 10.06 µs = *98,352 core-hours ≈ 11.2 core-years*.
- *Shortcut (measured).* $y^2 = x^3 + x$ is supersingular over this prime, with $\#E = p + 1$ and $v_2(p+1) = 64$. One draw gives a point of order $2^48$ (or $2^64$), and the checker accepts it. Whether PKM's relation phase accepts a supersingular $E'$ is untested.

= Method

For $p equiv 3 space (mod 4)$ and $E: y^2 = (x - e_1)(x - e_2)(x - e_3)$, the Weil pairing rules out $ZZ\/4 times ZZ\/4$, so $E(FF_p)[2^infinity] ≅ ZZ\/2 times ZZ\/2^b$. A point $P$ lies in $2E$ iff all $x - e_i$ are squares. The halves of $plus.minus P$ have $x_H = x + r_1 r_2 + r_1 r_3 + r_2 r_3$ with $r_i = plus.minus sqrt(x - e_i)$, $x_H - e_i = (r_i + r_j)(r_i + r_k)$ and $y_H = plus.minus (r_1+r_2)(r_1+r_3)(r_2+r_3)$. Negating $r_n$ gives $H + T_n$. The descent map $delta : E\/2E -> {plus.minus 1}^2$ reduces the question "can some half be halved again?" to one symbol $chi(x_H - e_k)$. The `fast` arm reads that symbol from the exponentiation the next level needs anyway, $i_k = (x_H - e_k)^((p-3)\/4)$ (189 squarings). A second exponentiation $u$ decides between $H$ and $H + T_n$ and yields the next roots without inversions.

#figure(image("figures/halving_chain.svg", width: 100%), caption: [One draw of the `fast` arm. Measured per `z2xz8` draw (C4 counters): 4.0 exponentiations, 1.0 Legendre symbol, 1.0 halving.])

= Cost per draw (C1, C2)

#figure(
  table(columns: 7, align: (left, left, right, right, right, right, right),
    [family], [arm], [median µs], [min], [max], [Legendre/draw], [exps/draw],
    [z2xz8], [baseline], [30.83], [29.30], [32.05], [8.49], [6.00],
    [z2xz8], [fast-chi], [12.48], [12.15], [13.30], [4.00], [2.00],
    [z2xz8], [*fast*], [*10.06*], [9.78], [10.19], [1.00], [4.00],
    [z2xz2], [baseline], [24.74], [24.00], [25.72], [9.15], [4.49],
    [z2xz2], [fast-chi], [16.81], [16.70], [17.63], [5.24], [2.99],
    [z2xz2], [fast], [13.89], [13.52], [14.96], [3.00], [3.93],
  ),
  caption: [C1: 65,536 draws × 7 interleaved repetitions, one core (`taskset -c 2`) of a shared 4-core host (load 2.3–3.5). Pre-registered: S1 fast/baseline 0.326 ≤ 0.40 ✓; S2 fast/fast-chi 0.806 ≤ 0.85 ✓; S3 10.06 ≤ 12 µs ✓; S4 (`z2xz2`) 0.561 ≤ 0.65 ✓. No min–max spreads overlap.],
)

#figure(image("figures/per_draw_cost.svg", width: 85%), caption: [C1 medians with min–max bars.])

C2 (ns per call, median of 11): mul 10.6; dedicated squaring 5.1; $(p-3)\/4$ chain 1,808; `sqrt` 2,185; binary Jacobi 1,793; posdivsteps Jacobi 1,952; inversion 2,153; drawing $m$ 146. *Negative result:* the exact posdivsteps Jacobi (libsecp256k1-style) is not faster than the binary one at 192 bits. Interleaving two exponentiation chains was slower and was dropped before the protocol.

= Success probability (C3, C4)

#figure(image("figures/draws_vs_t.svg", width: 100%), caption: [Measured draws per success $1\/hat(P)(b >= t)$ with 95% intervals (C3: $N = 2^25$ `z2xz2`; C4: $N = 2^28$ `z2xz8`). Dashed: the pre-registered models, extrapolated to $t = 48$. Dotted: the card's earlier hypotheses. Stars: the supersingular curve (C6).])

#figure(image("figures/obs_vs_model.svg", width: 85%), caption: [Observed / model count of draws with $b >= t$, with the model's 99.9% Poisson band; only $t$ with model count ≥ 20 is tested. H1 and H2: consistent at every tested $t$. For `z2xz8` the ratio drifts above 1 for $t >= 21$ (up to 1.31 at $t = 25$) while staying inside its band. The counts are cumulative, so neighbouring ratios are not independent. This is recorded as an observation.])

#table(columns: 6, align: right,
  [$t$], [count $b >= t$], [model], [ratio], [mean draws to success], [model],
  [16], [32,894], [32,768], [1.004], [8,161 ± 45], [8,192],
  [20], [2,103], [2,048], [1.027], [127,621 ± 2,772], [131,072],
  [24], [147], [128], [1.148], [1,822,411 ± 157,208], [2,097,152],
  [28], [12], [8], [1.5], [(22.6 ± 5.2)·10⁶], [33.6·10⁶],
)
`z2xz2` (C3) gives $P(b = 1)$: 8,387,239 observed against 8,388,608 predicted. The deepest `z2xz8` draw has $b = 38$ (seed 196253408), a 1-in-128 upper-tail event under H2.

= Curves, checker and tower (C5–C7)

`backlog_check_aux2t --expect-t t` accepts the ordinary artifacts at $t = 16, 24, 28, 38$ and rejects $t+1$. It also accepts the supersingular artifacts at $t = 48$ and $64$. The tower $V = x(R' + ⟨ Q ⟩)$, built at P-192 on the $b = 38$ curve ($T = 8, 12, 16$) and on the supersingular curve ($T = 8, 16$), has $|V| = 2^T$ distinct abscissae. $L(x) = c$ holds on all of $V$ in the field backend and again over `BigUint`. None of 4096 random abscissae outside $V$ hits $c$. About half of $V$ lifts to P-192 (32,609 of 65,536 at $T = 16$).

= The $t = 48$ cost, and what it does not decide

$2^45$ `z2xz8` draws × 10.06 µs = 98,352 core-hours (95,617–99,624 over the C1 spread), about 12 core-hours per $2^32$-draw shard. `z2xz2` would need 362,005. The card's earlier 39,094 core-hours assumed $2^46$ draws at 2 µs. $ZZ\/2 times ZZ\/8$ buys 3 bits, not 4, because $v_2(p-1) = 1$ caps the non-cyclic part at $ZZ\/2$. The figure is derived: it extrapolates H2 seventeen bits past the deepest tested $t$ on one prime.

The supersingular curve removes the search for P-192, and for any prime with large $v_2(p+1)$, if a supersingular $E'$ is admissible for PKM. This study does not test that, nor whether a long tower helps PKM at all (encyclopedia §8 lead 1).

#v(0.5em)
#text(size: 8.5pt)[Data: `data/` (C1–C7 raw JSON), `derived.json` (from `examples/two_sylow_study_stats.rs`), figures from `figures/plot.py` (layout only). Code: `src/cryptanalysis/two_sylow_search.rs`, `src/bin/pkm_aux_search.rs`, `examples/{two_sylow_bench,p192_isogeny_tower,two_sylow_study_stats}.rs`.]
