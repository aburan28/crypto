# Protocol: the vertical ℓ-step in JMV self-reducibility, measured and reformulated

Frozen 2026-10-08, after the pilot in `README.md` and before any
instrument below is run.  **Structural and asymptotic study, toy sizes.**
`S`, end-to-end cost and ECDLP speedup are **unset** and no outcome
changes them.  The subject is the cost and the algebraic shape of one
vertical ℓ-isogeny between endomorphism-ring levels of an ordinary
isogeny class over `F_p`.

## Known at freezing

- The verdict table and the ℓ² obstruction lemma of `README.md`.
- Pilot: 13 of 14 attempted depth-1 instances, `ℓ ≤ 19`, pass the DLP
  transfer, the `F_{p^r}` torsion field and the `π − λ` eigenline checks.
- Registered P-256 has `f_π = 1` and no vertical step; every instance
  below is CM-constructed.

## Shared definitions

- Depth `d = v_ℓ(f_π)`; level `k` counts from the surface, `k = 0`
  surface, `k = d` floor; `O_k = Z + ℓ^k O_K` locally at `ℓ`.
- For a curve at level `k < d`, `θ_k = (π − λ_k)/ℓ^{d−k}` with `λ_k ∈ Z`
  chosen so that `θ_k ∈ O_k`; `λ_k ≡ t/2 (mod ℓ^{d−k})` and the choice is
  unique modulo `ℓ^{d−k}`.
- `r = ord(t/2 mod ℓ)`, the degree of the field over which the surface
  ℓ-torsion is rational.
- Cost unit: wall seconds in Sage 10.10 on one core, reported with the
  machine; and, where a count is available, `F_p` multiplications.
- CM construction: `p = (t² − f² D_K)/4` prime with `D_K ∈ {−8, −11, −20}`,
  `f = ℓ` or `ℓ²`, `t` odd; the surface curve from the Hilbert class
  polynomial of `D_K` with the twist of trace `t`.

## J1: depth-2 volcanoes, intermediate level (conservative)

Instrument: extend `jmv_floor_toy.sage` with `f = ℓ²`, `ℓ ∈ {3, 5, 7}`,
two instances per `ℓ`, seed `20261008`; also rerun the aborted pilot
instance `(ℓ, p) = (17, 659)`.  At level 1 all `ℓ + 1` edges are rational;
classify each codomain by its own edge count (floor: 1; surface: `ℓ + 1`
with `h` horizontal) and identify the ascending edge.  Evaluate `θ_1` on
`E[ℓ]` by lifting: `θ_1(P) = (π − λ_1)(P′)` for `P′ ∈ E[ℓ²]` with `ℓP′ = P`.

- **J1.1 (eigenline).**  The ascending kernel at level 1 is the unique
  eigenline of `θ_1` on `E[ℓ]`.  *Falsified by one instance.*
- **J1.2 (ℓ² obstruction at level 1).**  Among `a + bθ_1`, `0 ≤ a, b <
  ℓ²`, every element killing the ascending kernel has norm `≡ 0 (mod
  ℓ²)`, and the set of such `(a, b)` has exactly `ℓ³` elements (index `ℓ`
  in `O_1/ℓ²O_1`).  *Falsified by one element of norm `≢ 0`.*
- **J1.3 (pilot completion).**  `(17, 659)` passes the three pilot checks
  with 600 draws.

## J3: cost of the torsion route in `r` and `ℓ` (conservative)

Instrument: new Sage script, CM classes with `f = ℓ`, `ℓ ∈ {101, 211, 401,
809, 1601, 3203}`, `D_K = −8`.  For each `ℓ` choose `t` to realise `r ∈
{1, 2, 3, 4, 6}` and one `t` with `r ≥ (ℓ − 1)/4`, `p` of 48 to 64 bits.
Time: (a) find a point of order `ℓ` in `E₀(F_{p^r})` by cofactor
multiplication; (b) the descending isogeny from that point with
`algorithm="velusqrt"`; (c) baseline, the kernel polynomial by factoring
the ℓ-division polynomial over `F_p`, for `ℓ ≤ 401` only.  Three timings
each, median reported.  Output `results/j3.jsonl`.

- **J3.1 (exponent in `r`).**  At each `ℓ`, the fitted exponent of (a)+(b)
  in `r` over `r ∈ {1, 2, 3, 4, 6}` lies in `[1.5, 2.5]`.
- **J3.2 (exponent in `ℓ`).**  At `r = 2`, the fitted exponent of (b) in
  `ℓ` lies in `[0.4, 1.2]`.
- **J3.3 (crossover).**  For every `ℓ ≥ 101` and `r ≤ 4`, (a)+(b) is below
  (c).  *Falsified by one `(ℓ, r)` where the division polynomial wins.*
- **J3.4 (large `r`).**  At the `t` with `r ≥ (ℓ − 1)/4`, (a)+(b) exceeds
  the `r = 2` time by at least `(r/2)^{1.5}`.

Decision rule: J3.3 and J3.4 pass: the vertical step at a large prime is
cheap exactly on the sparse set of classes with small `r`, mirroring C4 of
PR #1565, and the JMV gap is reduced to the large-`r` case.  J3.1 fails
low: `F_{p^r}` arithmetic is not the bottleneck and the model is revised
before J4.

## J4: replicate the `q^{1/4}` vertical step of [Gal24] (representation-changing)

Instrument: not yet written.  Spec: implement the algorithm of ePrint
2024/924 for the flat-volcano case at toy sizes, `q ∈ {2²⁰, 2²⁴, 2²⁸, 2³²}`,
CM classes with `f = ℓ`, `ℓ` the largest prime below `q^{1/2}/8`, three
instances per `q`.  The paper's sections that define the algorithm are
read and cited in the instrument header before any code; this protocol
does not assume knowledge of the method beyond the abstract.

- **J4.1 (exponent).**  Fitted exponent of wall time in `q` lies in `[0.2,
  0.35]`.
- **J4.2 (correctness).**  Every output is replayed by the pilot's checks:
  the codomain has the right edge count and the DLP transfer pulls back.
- **J4.3 (comparison).**  At `q = 2³²` the replicated step is faster than
  J3's torsion route at the same instance's `r`.  *Not a falsification
  of [Gal24] if it fails; a statement about this implementation.*

## J5: the invertible-ideal reformulation (representation-changing, then speculative)

Instrument: extend the pilot at depth 1, `ℓ ∈ {5, 7, 11}`, `D_K` chosen so
that each of inert, split and ramified `ℓ` occurs.  For each floor curve
`E₁′ ≠ E₁` under `E₀`, form `ψ = φ̂₁′ ∘ φ₁`, compute a generator `K` of its
kernel over the field where it lives, the scalar `μ (mod ℓ²)` with `π(K) =
μK`, and the kernel ideal `(ℓ², π − μ)` as a binary quadratic form of
discriminant `D_π`.

- **J5.1 (invertibility pattern).**  Exactly `ℓ − h` of the `ℓ` cyclic
  order-`ℓ²` groups over `G` give primitive forms, `h = 1 + (D_K/ℓ)`; the
  primitive ones are exactly the floor-to-floor composites and the
  non-primitive ones the composites landing on surface curves.
  *Falsified by one miscounted instance.*
- **J5.2 (class-group evaluation).**  For each primitive form, its class
  in `Cl(O′)`, written as a product of small split primes of `O′` by
  brute force in the toy class group, acts on `E₁` through small-degree
  horizontal isogenies and lands on `E₁′`.  This is the toy stand-in for
  [PR23].  *Falsified by one wrong landing.*
- **J5.3 (obstruction for floor-to-floor maps; proof obligation).**
  Every element of `Hom(E₁, E₁′)` killing `G` has degree `≡ 0 (mod ℓ²)`.
  Toy check: among the `ψ`-type composites and their compositions with
  small horizontal steps, none has degree `ℓ · d` with `gcd(ℓ, d) = 1`.
  The proof, if the check passes, goes through `Hom(E₁, E₁′) · I(G) =
  ℓ · Hom(E₁, E₁′) O_K`; it is written before J5.4.
- **J5.4 (speculative, no prediction).**  The gcd problem: from efficient
  representations of `ψ, ψ″` with `ker ψ ∩ ker ψ″ = G`, obtain `φ_G`.
  Registered only as a problem statement; a dimension-2 attempt from
  `2^e`-torsion images over the curve's own small extension is admissible
  only after a written proof obligation naming the diamond and its
  degrees.

Decision rule: J5.1 and J5.2 pass: the vertical step is equivalent to the
gcd problem J5.4 modulo polynomial-time work, and that becomes the
registered open problem of this line.  J5.1 fails: the derivation in
`README.md` is wrong in the counting and is corrected before anything
else is built on it.

## Stop conditions and inadmissible moves

Wall-time caps: J1 30 minutes; J3 10 minutes per `(ℓ, r)` timing, 4 hours
total; J4 2 hours per `q`; J5 1 hour.  A size that does not finish is
reported as partial.

Inadmissible: changing the exponent windows after a fit; dropping
instances; choosing `t` after seeing `r`-dependence except as the
registered design says; reading any result as a statement about P-256,
whose class has no vertical step, or about `S`.
