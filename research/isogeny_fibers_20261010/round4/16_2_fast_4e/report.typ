#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8.5pt)

#align(center)[
  #text(size: 14pt, weight: "bold")[Lead 16.2: a fast 4E test as a subgroup check]\
  Round 4 of `RESEARCH_ISOGENY_FIBERS_ECDLP.md` (§14.2, §16.2) · 2026-10-11 · worktree commit `440c0e163` + untracked round-4 files\
  Code: `examples/r4_16_2_fast_4e.rs` · runs started 2026-10-11T19:50:38Z · x86_64, single thread, contended host
]

= Question and verdict

On the a = 0 Koblitz curve $y^2 + x y = x^3 + 1$ over $FF_(2^n)$, $n$ odd, with $\#E = 4r$, does the closed form of Prop. 14.2 decide membership in $G = 4E$ exactly, and what does it cost against $[r]P = O$?
The falsifier is any disagreement.

*Verdict: confirmed (measured).*
- 0 errors on every point for $n = 5, 7, dots, 23$ (exhaustive).
- 0 errors on $10^6$ sampled points each at $n = 131$ (the ECC2K-130 curve) and $n = 233$ (sect233k1), across all four cosets of $4E$.
- On members of $G$ the test costs 1 M, 3 traces, and either $n-1$ S or $ceil(n\/8)$ table lookups. τ-NAF $[r]P$ costs 346 M + 609 S ($n=131$) and 618 M + 1085 S ($n=233$).
- *No novelty is claimed.* A two-level halving validation for cofactor 4 already appears in the literature (§5).

= The test (derived)

Doubling $(u, v)$ gives $lambda = u + v\/u$, $x = lambda^2 + lambda$ and $y = u^2 + (lambda+1)x$.
- $P$ has a rational half iff $"Tr"(x) = 0$. For odd $n$, $lambda = H(x)$ or $H(x)+1$, where $H$ is the half-trace.
- A half $Q$ lies in $2E$ iff $"Tr"(u) = 0$. Since $"Tr"(u) = "Tr"(u^2) = "Tr"(y) + "Tr"(lambda x) + "Tr"(x)$, this is $"Tr"(y) + "Tr"(x H(x))$ when $"Tr"(x)=0$.
- The two choices of $lambda$ give $Q$ and $Q+T$. They agree because $T = 2 dot (1,0) in 2E$.

#block(inset: 6pt, stroke: 0.5pt + gray)[
  $P in 4E <==> P = O "or" ["Tr"(x) = 0 "and" "Tr"(y) + "Tr"(x dot H(x)) = 0]$.
  The traced elements are $x$, $y$ and the product $x dot H(x)$.
]

- $T = (0,1)$ is rejected automatically, since $"Tr"(1) = n mod 2 = 1$.
- $(1, 0)$ has order 4 and no half, because $"Tr"(1) = 1$. So the 2-part is exactly $ZZ\/4$ for every odd $n$, and $4E = G$ when $r$ is prime.
- Preconditions: $P$ is on the curve, $n$ is odd and $a = 0$.

= Measured results

#figure(
  image("costs.svg", width: 100%),
  caption: [Counted field operations (left, middle) and contended wall time (right) per check, for members of $G$.
    Counts come from 249 896 ($n=131$) and 250 583 ($n=233$) points. Wall time is the median of 5 reps over 2048 points.
    Sources: `sample_n131.json`, `sample_n233.json`, `plot_data.json`.],
)

#table(
  columns: 7,
  align: (right, left, right, right, right, right, right),
  table.header([n], [true coset], [in G], [accept], [reject], [bits (Tr x, Tr y + Tr xH)], [errors]),
  [131], [0], [yes], [249 897], [0], [(0,0)], [0],
  [131], [1], [no], [0], [249 504], [(1,0)], [0],
  [131], [2], [no], [0], [250 049], [(0,1)], [0],
  [131], [3], [no], [0], [250 559], [(1,1)], [0],
  [233], [0], [yes], [250 584], [0], [(0,0)], [0],
  [233], [1], [no], [0], [249 533], [(1,1)], [0],
  [233], [2], [no], [0], [249 993], [(0,1)], [0],
  [233], [3], [no], [0], [249 899], [(1,0)], [0],
)

Exhaustive runs found 0 false accepts and 0 false rejects at every odd $n$ from 5 to 23; the largest is $n = 23$ with 8 383 412 points.

Further checks:
- Binary and τ-NAF $[r]P$ agree on every point.
- The library's `is_valid_public_point` agrees on 2000 of 2000 points per $n$.
- $r(131)$ equals the literature ECC2K-130 order, cross-checked by search only. $r(233)$ equals the library's sect233k1 order.

*Observation (not derived):* on these samples the bit pair also separates cosets 1 and 3, with the assignment swapped between $n = 131$ and $n = 233$.

#table(
  columns: 6,
  align: (right, left, right, right, right, right),
  table.header([n], [method], [M], [S], [traces / lookups], [contended wall, median]),
  [131], [$[r]P$ binary LD], [661], [816], [0], [232 µs],
  [131], [$[r]P$ τ-NAF LD], [346], [609], [0], [160 µs],
  [131], [closed form, $H$ by squarings], [1], [130], [3 / 0], [22.8 µs],
  [131], [closed form, $H$ by table], [1], [0], [3 / 17], [0.124 µs],
  [233], [$[r]P$ binary LD], [1255], [1506], [0], [165 µs],
  [233], [$[r]P$ τ-NAF LD], [618], [1085], [0], [98 µs],
  [233], [closed form, $H$ by squarings], [1], [232], [3 / 0], [11.8 µs],
  [233], [closed form, $H$ by table], [1], [0], [3 / 30], [0.122 µs],
)

The counts match the derived formulas: binary $3(L-1) + 8(w-2) + 2$ M; τ-NAF $8(z-2) + 2$ M; the closed form $1$ M + $(n-1)$ S.
Multiplication ratio, τ-NAF against the closed form: 346× and 618×.
Wall-time ratio, τ-NAF against the table variant: about 1300× and 800×. These are contended, with untuned squaring, so they are order-of-magnitude only.

= Limits

- The exhaustive check stops at $n = 23$. A run to $n = 25$ was stopped to respect the per-run time limit.
- Not constant-time. Odd $n$ and $a = 0$ only.
- There is no `identity.certificate/v1` record for the identity.
- `cargo clippy -D warnings` fails on 28 errors in `src/` at this commit. The example itself has no lints. `rustfmt --check` passes.

= Literature status (search snippets only; PDFs blocked)

- Halvability $"Tr"(x) = "Tr"(a)$: Knudsen 1999; Schroeppel; Fong–Hankerson–López–Menezes.
- US patent 7,257,709 / US 2004/0114760 (priority 2002): validation of binary-curve keys with cofactor 2 or 4 by testing for a half via traces.
- King, PKC 2004 (LNCS 2947).
- Eagle–Galbraith–Ong, ePrint 2009/086. For $\#E = 4 ell$: halve, then check $"Tr"(x_Q) = "Tr"(a^2)$.
- Pornin, ePrint 2022/1164: membership by halvings on twisted Edwards curves.
- ePrint 2026/749 (Divide-and-Pair).

The closed form here only inlines the half's trace, so no square root or second half-trace is needed. No novelty is claimed.
