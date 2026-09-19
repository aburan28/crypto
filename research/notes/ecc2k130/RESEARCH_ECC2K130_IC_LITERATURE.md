# What the literature has, and has not, for index calculus on ECC2K-130

**Method:** automated multi-source survey, 103 agents, 21 sources fetched,
105 candidate claims extracted, top 25 put through 3-vote adversarial
verification. **12 confirmed, 13 refuted.** The refutation rate is the
useful part: most of what died were over-strong claims on the *negative*
side, which is why the conclusion below is "stalemate", not "settled".

**Target:** `E : y² + xy = x³ + 1` over `F_{2^131}`, scored against
parallel Pollard rho with the `⟨−1⟩ × ⟨π⟩` speed-up, **2^60.8** group
operations on the prime-order subgroup.

**Excluded as already built here** (the survey was told to report on
these only to say how a new idea differs): Gaudry/Diem subspace factor
bases and symmetrised `S₄`; Frobenius-invariant *and* Frobenius-union
factor bases with the GGMP collapse; quasi-subfield polynomials; GHS and
hyperelliptic covers; F4/XL and CDCL SAT on the descended systems; the
`u = 1/(x+1)` 2-torsion frame.

## Bottom line

**No published approach outside the exclusion set gives a numerically
explicit route that beats 2^60.8 at `n = 131`.** For the two families
examined in depth the reasons are now sharp rather than hand-waved, and
one of them is permanent. Three of the six research items came back
*empty*, which is the actionable part — see [§5](#5-the-open-ground).

## 1. The first-fall-degree controversy is a stalemate, not a negative

Every subexponential small-characteristic ECDLP claim is gated on the
first fall degree assumption, `D_reg = D_ff + o(1)`. The survey confirms
the assumption is **positively suspect**, not merely unproven:
Huang–Kosters–Yeo (CRYPTO 2015) give a reductio — chain `m = O(n)`
copies of `S₃` with x-coordinates in an `O(n/m)`-dimensional subspace and
you get a degree-≤3 Weil descent system with `D_ff ≤ 5` (usually 2), so
the assumption would make ECDLP polynomial-time. Their conclusion:

> "Consequently, we have a polynomial-time algorithm to solve ECDLP,
> which is highly improbable. We conclude that the first fall degree
> assumption is unlikely to hold for this system as well."

Kousidis–Wiemers (JMC 2019) independently downgrade Petit–Quisquater's
Assumption 2 to an open *question* and replicate Kosters–Yeo's anomalous
rise in `D_reg` at `m = 2`. Galbraith–Gaudry (DCC 2016) §10.2 record
that there is **"no consensus whether there is a subexponential
algorithm for ECDLP in characteristic 2."** Nothing in 2020–2026
rehabilitates it.

**But the refutation side is heuristic too**, which is what makes this a
stalemate worth caring about:

- HKY's conclusion is explicitly scoped — "unlikely to hold **for this
  system**" — and their own §5.3 says the degree of regularity of a Weil
  descent system "tends to grow more slowly than a random system… The
  slower it grows, the better algorithms there will be for ECDLP". The
  authors do not claim the line is dead.
- Kosters–Yeo put "prove" (P=NP) in scare quotes; Theorem 3.5i is
  conditional on Assumption 3.3, Remark 3.6 notes the unconditional
  variant fails for `p = 2`, and Remark 3.7 says outright: "it does not
  suggest that ECDLP itself is a hard problem."
- Rigorously bounding last fall degrees of Weil descent systems is
  labelled an open problem in 2015 and is **still open in 2026**.

*Consequence for this repository:* neither a rigorous upper nor lower
bound on the cost of solving these systems exists, so **our own measured
`D_reg` growth at small `n` is as good as the published state of the
art.** The FFD ladder in `research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md` is not
duplicating a settled question.

## 2. The number that actually decides it

Kousidis–Wiemers Theorem 3.2 proves **unconditionally** that
`D_ff ≤ m² − m + 1` for the Weil descent of `S_{m+1}` over `F_{2^n}`
(`n' ≥ m ≥ 3`), improving Petit–Quisquater's `m² + 1`, and it is
empirically sharp (7 at `m = 3`, 13 at `m = 4`). The sharpened
complexity is `O(2^{c·log(n)·(n^{2/3} − n^{1/3} + 1)})` with
`c = 2ω/3`.

| | value |
|---|---|
| turning point vs generic `O(2^{n/2})` | `n ≈ 1250` (PQ's was ≈ 2000) |
| the formula at `n = 131` | **2^86.0** |
| generic `2^{n/2}` | 2^65.5 |
| rho reference | **2^60.8** |

So even **granting the suspect assumption** and using the best proven
bound, the algebraic line is about **2^25 short** at our size. The
verifier recomputed the exponent independently with `ω = log₂ 7`; only
the log₁₀ convention reproduces both of the paper's stated turning
points (`n = 1250 → 615.9` vs 625; PQ at `n = 2000 → 986.9` vs 1000), and
**the gap exceeds 2^20 under every reading**, far beyond any suppressed
constant.

One logical point worth pinning, because it is easy to get backwards:
**a `D_ff` upper bound gives no `D_reg` upper bound at all** — `D_ff`
bounds `D_reg` from *below*. Theorem 3.2 is a theorem; the complexity
that follows from it is assumption-gated.

Note also the optimal `m` for that formula is `≈ n^{1/3} ≈ 5.1` at
`n = 131`. Our harness reaches `m = 3` and strains at `m = 4`. That gap
— between the `m` we can run and the `m` the (conditional) theory needs
— is the real practical barrier, and it is the metric
`research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md` already tracks.

## 3. Rigorous theory certifies nothing useful here

The only general rigorous bound on the solving degree of a Weil
restriction is Caminata–Ceria–Gorla, `sd(Weil(F)) ≤ n·reg(F^h) − n + 1`,
which at `n = 131` gives **263** (reg 3) or **394** (reg 4). Adding the
Boolean field equations does not improve it (their Cor. 4.12 gives the
identical number). It is numerically vacuous: for `m = 2` there are 262
Boolean unknowns, so a bound of 263 exceeds the maximum possible monomial
degree and is formally content-free. The paper contains no ECDLP
instance, no field size and no operation count — a full-text search
returned zero hits for "ECDLP", "subexponential", "Petit", "Koblitz".

## 4. Two things to stop doing

**Do not cite Yokoyama et al. (JMC 14, 2020) as closing this question.**
It looks like a lower bound that ends the argument. It is not: it is
conditional on unproven Assumptions 7 and 11 plus unproved Conjecture 10
("We leave a theoretical proof of the conjecture for our future work"),
its title and abstract restrict it to Semaev's ***naive*** index calculus,
its experiments are over **prime fields only** (B = 10–25 bits, validating
assumptions rather than costing an attack), and its bound is parametric
and never instantiated at a cryptographic size. Its `cited_by_count` is
0 — neither independently validated nor refuted.

Worse for our purposes, its §5.4 carves out exactly two exceptions —
Diem's subspace factor bases and HKPYY quasi-subfield polynomials —
**both of which are already in our exclusion set**. The literature's own
pointer to "index calculus can beat rho" leads straight back to what we
have already built.

**Cross off Joux–Vitse cover-and-decomposition permanently.** It is
structurally inapplicable, not merely inefficient: the algorithm is
defined on a tower `F_{q^d}/F_q/F_p` and needs `q` a strict power of `p`,
i.e. a **composite** extension degree. With `n = 131` prime there is no
intermediate field, so neither the cover stage nor the decomposition
stage has a target. Galbraith's ellipticnews post (2012-02-06):
"For prime values of `n`… one cannot combine the two methods, and so the
paper does not discuss these cases." The obvious evasion was tested and
fails: embedding `E(F_2^131)` into `F_2^{131m}` collapses both towers —
with `q = 2^131` the cover step is trivial and decomposition reverts to
the original problem over a larger field, and with `d = 131` the GHS
cover has genus `≈ 2^130`, which is exactly the constraint
`research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md` derived independently. No 2011–2026
work removes the composite-degree condition.

Its "more relations than unknowns" trick is worth one line of caution:
the `(ng+2)`-point variant has an **asymptotically worse** exponent,
`2 − 2/(ng+2)` against `2 − 2/(ng)`, and is faster only on measured
constants (~2.5× at `log₂ p ≈ 27`). Any transplant of that idea must be
argued on measured constants, never on asymptotics. The stronger claim —
that its precompute-the-Gröbner-basis-once structure is the key
transferable idea — was **voted down 0-3** in verification and should be
treated as unverified.

## 5. The open ground

**Research items 3, 4 and 5 produced zero surviving claims.** The survey
flags this honestly as absence of evidence in the collected corpus, not
evidence of absence — none of the 12 confirmed and none of the 13
refuted claims touch any of:

- **Crossbred, BDD/DDR, hybrid XL, or modern MQ solvers** applied to
  binary ECDLP decomposition systems;
- **index-calculus variants at prime extension degree** needing neither a
  subfield nor an invariant subspace — Semaev's later proposals, Nagao's
  decomposition as a standalone method, Diem's higher-dimensional
  constructions, extension-field Jacobian decomposition;
- **Koblitz structure beyond Frobenius** for index calculus rather than
  rho or scalar multiplication — the τ-adic expansion, CM by
  `(1 ± √−7)/2`, the class group of `Z[τ]`.

The hard structural constraints (`131` prime, `2` primitive mod `131`,
genera `1, 2^129, 2^130` only) pre-filter much of item 4. **They do not
touch item 3 at all** — the choice of MQ/hybrid solver is orthogonal to
the factor-base construction — and they do not obviously touch item 5.

Item 3 is the shortest path, because half of it is already here:
`src/cryptanalysis/crossbred.rs` implements Joux–Vitse Crossbred and is
already tested against Koblitz decomposition systems at `m = 2` and
`m = 3`. The survey found **no publication combining the two**. That is
the one place where "novel idea not yet seen" and "already half-built"
coincide.

## Honest limits

- This is a literature survey, not a result. It says what has been
  published, at 21 sources; it cannot prove nothing exists.
- Items 3–5 are *under-searched*, and the null there is weaker evidence
  than the confirmed findings in §§1–4.
- Every number in §2 is the published authors', recomputed for
  convention but not re-derived from scratch.
- Nothing here threatens a deployed curve, and nothing here claims to.

## Sources

- M.-D. Huang, M. Kosters, S. L. Yeo, *On the last fall degree of zero-dimensional Weil descent systems*, CRYPTO 2015 — [eprint 2015/573](https://eprint.iacr.org/2015/573)
- M.-D. Huang, M. Kosters, C. Petit, S. L. Yeo, Y. Yun — [arXiv:1505.02532](https://arxiv.org/pdf/1505.02532)
- M. Kosters, S. L. Yeo — [arXiv:1503.08001](https://arxiv.org/abs/1503.08001)
- S. Kousidis, A. Wiemers, *On the first fall degree of summation polynomials*, J. Math. Cryptol. 2019 — [arXiv:1906.05594](https://arxiv.org/abs/1906.05594)
- S. Galbraith, P. Gaudry, *Recent progress on the elliptic curve discrete logarithm problem*, DCC 2016 — [eprint 2015/1022](https://eprint.iacr.org/2015/1022)
- A. Caminata, M. Ceria, E. Gorla, *Solving degree of Weil restrictions* — [arXiv:2112.10506](https://arxiv.org/abs/2112.10506)
- K. Yokoyama, M. Yasuda, Y. Takahashi, J. Kogure, J. Math. Cryptol. 14 (2020) 460–485 — [doi:10.1515/jmc-2019-0029](https://www.degruyterbrill.com/document/doi/10.1515/jmc-2019-0029/html)
- A. Joux, V. Vitse, *Cover and decomposition index calculus…* — [eprint 2011/020](https://eprint.iacr.org/2011/020)
