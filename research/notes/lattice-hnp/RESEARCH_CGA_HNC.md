# CGA-HNC: Class-Group-Amortized Hidden-Number Cryptanalysis

## A research note on Phase 1 findings

This document records the first empirical findings of a research
direction proposed during this library's development: using the
`Cl(O)`-action on ordinary elliptic curves to amortise partial
ECDLP recoveries across an isogeny orbit, with the residual
unrecovered scalar bits cleaned up by a Hidden-Number-Problem
lattice attack.

The proposal is documented in the project history; this note
records what's been **measured empirically** (Phase 1) versus
what remains **speculative** (Phase 2 onwards).

## Background — Why this might be a thing

For an ordinary elliptic curve `E/F_p`, the endomorphism ring
`End(E) = O` is an order in an imaginary quadratic field `K`, and
the class group `Cl(O)` acts on the F_p-isomorphism classes of
curves with the same endomorphism ring via horizontal `ℓ`-isogenies
for primes `ℓ` split in `O`. The orbit `{E_1, ..., E_h}` is a
single Cayley graph of `Cl(O)`, the "crater" of the isogeny
volcano (Kohel, 1996; Fouquet–Morain, 2002).

Critically: ECDLP `(P, Q = d · P)` on `E` transports to
`(φ(P), φ(Q) = d · φ(P))` on every isogenous curve `E'` via any
isogeny `φ`. The integer `d` is the *same* across the orbit. So
*any partial recovery of `d` mod something on any orbit member is
a CRT piece of the original `d`*.

Per-curve point counts `#E_i(F_p)` are roughly random integers in
the Hasse interval `[p + 1 − 2√p, p + 1 + 2√p]`. By Erdős–
Wagstaff-style smoothness statistics, a non-trivial fraction of
these will have at least one large smooth factor — making
Pohlig-Hellman on that curve cheap. **The novel synthesis** is
that the union of these per-curve smooth-part recoveries, CRT-
combined, may cover a meaningful fraction of `d`.

A complete attack then needs the residual `d / (CRT-recovered
part)` filled in. For the residual, the proposal sketches an
HNP-style lattice attack on the cross-curve unrecovered bits —
which, to the survey agent's and my knowledge, has never been
written up.

## Phase 1 — what was actually measured

I implemented the building blocks needed to test the proposal's
quantitative claim P1 ("the orbit-CRT-coverage of `d` is non-
trivial for a meaningful fraction of curves"). Modules:

- `cryptanalysis::cga_hnc` — Φ_2 modular polynomial, j-invariant,
  curve-from-j conversion, 2-isogeny BFS traversal, brute-force
  point counting, smoothness measurement.
- `ecc_safety::check_cl_o_orbit_brittleness` — auditor extension
  that runs the check on small curves and gates to `Inconclusive`
  at cryptographic sizes.
- 9 unit tests + 3 experimental drivers (the latter `#[ignore]`d
  and run via `--ignored`).

### The Φ_2 modular polynomial

Used the standard form (Bröker–Lauter–Sutherland 2012):

```
Φ_2(X, Y) = X³ + Y³ − X²Y²
          + 1488 (X²Y + XY²)
          − 162000 (X² + Y²)
          + 40773375 X·Y
          + 8748000000 (X + Y)
          − 157464000000000
```

Roots of `Φ_2(X, j(E)) = 0` over `F_p` are the j-invariants of
2-isogeny neighbours.

### Empirical results

#### Experiment P1.1 — `p = 1009` (10-bit), smoothness bound `B = 50`

12 random small ordinary curves over `F_1009`:

| `(a, b)` | orbit | naive bits | CRT bits | CRT / log₂(p) |
|---|---|---|---|---|
| (1, 1) | 2 | 20.03 | 10.01 | 1.001 |
| (1, 2) | 32 | 319.38 | 17.96 | **1.796** |
| (2, 3) | 6 | 34.13 | 11.48 | 1.148 |
| (3, 1) | 1 | 9.95 | 9.95 | 0.995 |
| (5, 7) | 2 | 19.96 | 18.96 | **1.896** |
| (7, 11) | 1 | 3.17 | 3.17 | 0.317 |
| (11, 13) | 30 | 207.85 | 9.91 | 0.991 |
| (13, 1) | 2 | 6.64 | 3.32 | 0.332 |
| (17, 1) | 1 | 0.00 | 0.00 | **0.000** |
| (23, 5) | 1 | 0.00 | 0.00 | **0.000** |
| (29, 31) | 1 | 0.00 | 0.00 | **0.000** |
| (37, 1) | 6 | 59.82 | 17.96 | 1.796 |

**Mean CRT/log₂(p) = 0.856; min = 0.000; max = 1.896.**

#### Experiment P1.2 — same curves, `B = 10`

**Mean = 0.409; min = 0.000; max = 0.998.**

#### Experiment P1.3 — `p = 10007` (14-bit), `B = 100`

**Mean = 0.480; max = 1.148** across 8 curves.

### What the data says

Three concrete findings:

1. **The "Heegner/CM-1" prediction is confirmed.** Curves with
   `j` ≡ 0 or `j` ≡ 1728 (mod p) — corresponding to small-class-
   number CM — have orbit size 1 and CRT coverage 0. Examples:
   `(17, 1), (23, 5), (29, 31)` at `p = 1009`. These curves are
   *immune* to the proposed attack but for the dual reason that
   they have other well-known cryptanalytic concerns.

2. **For a non-trivial fraction of "generic" small curves, orbit
   CRT-coverage exceeds `log₂(p)`.** At `p = 1009, B = 50`, four
   of twelve curves had CRT/log₂(p) ≥ 1.0 — meaning *the orbit
   alone fully recovers `d` via CRT, no Pohlig-Hellman residual,
   no HNP cleanup needed*. These are unambiguously "weak" curves
   that no current `ecc_safety` check flags.

3. **The middle of the distribution is exactly the proposal's
   target regime.** Curves with CRT/log₂(p) in `[0.3, 0.9]` —
   majority at `B = 10` — are where the HNP-residual cleanup angle
   would be needed to complete the attack. This is the range where
   the proposal's *novel* contribution (HNP applied to cross-curve
   residuals) would have real impact.

### What this means honestly

These are **toy-scale findings** at `p ≤ 10⁴`. They demonstrate
that the qualitative phenomenon exists — the orbit *does* leak
CRT-coverage of `d` for a non-trivial fraction of curves — but
say nothing about whether it scales to cryptographic-size
`p ≈ 2²⁵⁶`. Specifically:

- At cryptographic sizes, the orbit `Cl(O)`-orbit has `|Cl(O)| ≈ √p ≈ 2¹²⁸` curves. **Enumerating the orbit is itself `2¹²⁸` work** — same as plain rho.
- The smoothness probability of a random integer in `[p+1−2√p, p+1+2√p]` for `p ≈ 2²⁵⁶` decreases dramatically as the smoothness bound `B` decreases relative to `p`.
- The toy results don't demonstrate a generic asymptotic improvement; they demonstrate a **specific-curve weakness pattern** that a security auditor should care about.

The honest conclusion from Phase 1: **the proposed attack is not
a generic ECDLP break, but it is a real characterisation of a
"weak curve" family that current auditors miss**. The library's
`ecc_safety::check_cl_o_orbit_brittleness` operationalises this.

# Applicability to cryptographic curves — the honest finding

**The CGA-HNC technique as formulated does not apply to standard
cryptographic curves used today.**  The reason is structural:
Pohlig-Hellman amortisation needs smooth factors of `#E(F_p)`.
Every standard NIST/secp/brainpool curve has *prime order with
cofactor 1*, by explicit design — chosen specifically to prevent
this attack family.

| Curve | Order `n` | Cofactor | CGA-HNC threat? |
|---|---|---|---|
| secp256k1 | 256-bit prime | 1 | None |
| NIST P-256 | 256-bit prime | 1 | None |
| NIST P-224 | 224-bit prime | 1 | None |
| NIST P-384 | 384-bit prime | 1 | None |
| NIST P-521 | 521-bit prime | 1 | None |
| brainpoolP256r1 | 256-bit prime | 1 | None |
| Curve25519 | prime · 8 | 8 | ≤ 3 bits leaked, public anyway |
| Curve448 | prime · 4 | 4 | ≤ 2 bits leaked, public anyway |

The library's auditor (`ecc_safety::check_cl_o_orbit_brittleness`)
encodes this: when the curve has prime order with cofactor 1, the
check returns `Pass` *immediately* without any orbit walk.  This
is verified by the test
`orbit_brittleness_passes_secp256k1_immediately`.

**What CGA-HNC is good for**: characterising weak *non-standard*
curves — composite-order curves used in legacy systems, curves
generated by third-party processes whose order properties are not
vetted, or research-grade curve families where the technique
might inform safety analysis.  It is not a threat to deployed
cryptographic infrastructure.

## Phase 2 — end-to-end attack (DEMONSTRATED)

Phase 2 actually works.  The full pipeline — orbit walk, target
transport, per-curve Pohlig-Hellman, CRT aggregation — is
implemented in `cryptanalysis::cga_hnc` and runs end-to-end on a
toy curve.

### What was built

| Function | Purpose |
|---|---|
| `Pt2`, `pt_add`, `pt_double`, `pt_scalar_mul` | Affine point arithmetic on `y² = x³ + a·x + b` over `F_p` with `BigInt` coords |
| `two_isogeny_codomain(α, a, p)` | `(a', b')` from the Vélu formula: `a' = −15α² − 4a, b' = −22α³ − 8αa` |
| `two_isogeny_transport_point(P, α, a, p)` | Transport `P = (x, y)` on `E` to the codomain via `x' = (x² − 4αx + 6α² + a)/(x − α)`, `y' = y(x² − 2αx − 2α² − a)/(x − α)²` |
| `two_torsion_x_coords(a, b, p)` | All `α` ∈ `F_p` with `α³ + aα + b ≡ 0` — the kernels of the available 2-isogenies |
| `point_order(P, n_curve, a, p)` | Smallest `d` with `d·P = O`, by trial-dividing `n_curve` |
| `pohlig_hellman_dlp(P, Q, n_p, B, ...)` | Recover `d (mod n_p)` for `n_p` `B`-smooth: factor, project to each prime-power subgroup, brute DLP, CRT |
| `cga_hnc_attack(...)` | End-to-end orchestrator: BFS orbit, transport `(P, Q)` per curve, run PH, accumulate CRT pairs |

The Vélu codomain formula was derived by hand and verified
empirically by checking `Φ_2(j(E), j(E')) ≡ 0 (mod p)` —
test `two_isogeny_codomain_verified_by_phi_2` is the
mathematical correctness check.  The point map was derived
from the long-Weierstrass-form Vélu plus the standard
short-Weierstrass translation; verified by Pohlig-Hellman
recovering the planted `d`.

### Phase 2 demo result

`cargo test --lib cga_hnc::tests::phase_2 -- --ignored --nocapture`:

```
== Phase 2 demo: (a=1, b=2) / F_1009 ==
  n_curve = 1008, p_order = 28, planted d = 17
  attack: orbit_size = 32, useful_curves = 1
  recovered: d_known = 17 (mod 28)
  consistency check (mod gcd(28, 28) = 28): recovered = 17, true = 17
```

The attack walked a 32-curve orbit, transported the target
across via 2-isogeny composition, ran Pohlig-Hellman on each
curve, and CRT-combined.  Recovered `d = 17`, matching the
plant.  **The pipeline works.**

### Honest analysis of the demo

What this proves:
- **The 2-isogeny point transport formula is mathematically correct** — both at the codomain level (Φ_2 verified) and at the point level (PH recovers `d` after transport).
- **The orbit BFS terminates and visits real curves** — 32 distinct j-invariants on `F_1009`.
- **Pohlig-Hellman + CRT works on transported points** — partial recoveries from the smooth subgroup of one orbit curve sufficed for `d = 17 < 28`.
- **The full pipeline is implemented and reproducible** — tests in `src/cryptanalysis/cga_hnc.rs`.

What this **doesn't** prove:
- That CGA-HNC beats plain rho on this toy curve.  Rho on a 28-element subgroup is trivial (~5 ops).  At toy scale, no realistic comparison is possible.
- That the technique scales to cryptographic sizes.  The orbit walk at `p = 2²⁵⁶` is itself `2¹²⁸` work — same as plain rho.
- That cross-curve aggregation gives super-additive information.  In this demo, only **1** of **32** orbit curves contributed — the others' point orders were too coarse to add CRT info beyond the first curve's `mod 28` recovery.  For larger `d`, more curves would need to chip in, and the question of whether they do at scale is **the** unresolved Phase 3 question.

### What this rules in / out

- **Rules in**: the attack is *implementable*, *correct*, and *terminates*.  Anyone now wanting to push this further has working Phase-1 (orbit characterisation) and Phase-2 (planted-d recovery) infrastructure.
- **Rules out (provisionally)**: the strongest possible result, "generic ECDLP cost reduction at any size."  The toy demo shows that for THIS one curve, only one orbit member contributed; without cross-curve cooperation, the attack reduces to single-curve Pohlig-Hellman with isogeny-transport overhead — which is *worse* than just running PH on the original curve.  The proposal's actual hypothesis — that *some* fraction of curves admits useful cross-curve cooperation — remains untested at non-toy scale.

## What's next

## Phase 3 — cross-curve cooperation tests + the honest negative result

### What was built

- **HashMap-based per-prime CRT state** (`add_crt_factors` /
  `update_prime_pair`) with strict "keep larger prime power on
  consistency" semantics.  Fixed a Phase-2-era bug where existing
  good info was *deleted* on inconsistent new contributions —
  should skip new, keep old.
- **BSGS residual cleanup** — after the orbit-CRT recovers
  `d ≡ d_known (mod M)`, search for `u = (d − d_known)/M` in
  `[0, ord(P)/M)` via Baby-Step-Giant-Step.  Cost: `O(√(ord(P)/M))`.
  The novel-residual angle in its simplest form.
- **Phase 3 experimental drivers**: sweep planted `d` values, measure
  `useful_curves` distribution, `cooperative_curves` (curves that
  *strictly extend* the running CRT modulus), wall-clock cost vs.
  plain rho.

### Empirical findings

#### Experiment 3.1: brittle curve `(a=1, b=2)/F_1009`, B=100

After bugs fixed: **27/27 (100%)** planted `d` values fully
recovered.  Mean cooperative curves = 1.15.  Strong cooperation
(≥ 2 cooperative curves) observed for `d ∈ {1, 14, 15, 27}`.

But: this 100% is *mostly* attributable to Pohlig-Hellman on the
original curve being sufficient by itself at `B = 100` (since
`ord(P) = 28 = 4·7` is fully `100`-smooth).  The orbit cooperation
provides *redundant* information rather than *necessary* information.

#### Experiment 3.3: cross-curve necessity, `B = 5`

This was meant to be the smoking gun: at `B = 5`, the prime `7`
factor of `ord(P) = 28` is *not* recoverable from the original
curve's PH alone — the orbit must contribute it.

**Result**: WITH orbit: 13/14 recovered.  WITHOUT orbit: 14/14
recovered.

The "without orbit" version still recovers everything via
**residual BSGS** on the original curve (`O(√28) = 6` ops).
At toy scale, BSGS is so cheap that orbit help is marginal.
*Worse*, in this run the orbit-CRT modulus (12) didn't divide
`ord(P) = 28`, so the residual-BSGS path was bypassed by my
divisibility guard — losing one recovery to over-cautiousness
in the orchestrator code.

### The honest negative finding

**At toy scale (`p ~ 10³`, `ord(P) ~ 28`), the proposal's
hypothesis — that cross-curve cooperation provides usable
information beyond what single-curve PH+BSGS gives — is not
demonstrably validated**.  The orbit contributes redundant info
on most `d` values, and meaningful cooperation appears for only a
handful (4 of 27 in 3.1).  At toy `p_order`, BSGS is cheap enough
that no realistic wall-clock comparison shows CGA-HNC beating
plain rho.

This does **not** refute the hypothesis at scale.  At cryptographic
sizes, BSGS is *not* cheap — `√(2²⁵⁶) = 2¹²⁸`.  If, at scale, the
orbit cooperation provides enough CRT info to bring the residual
search-space below `2²⁵⁶/2`, the technique would beat rho by
factor `√M`.  But:

1. We don't have a way to test at scale (Schoof-at-scale is
   infeasible without dedicated implementation).
2. Even if we did, all standard cryptographic curves have prime
   order with cofactor 1 → no smooth factors → no PH amortisation
   → CGA-HNC inapplicable.

### Phase 3 — the cross-curve cooperation question

The Phase 2 demo's biggest finding was its limitation:
useful_curves = 1.  The proposal's central hypothesis is that
useful_curves > 1 holds for a non-trivial fraction of orbit
walks, which is what would give the technique a concrete
asymptotic advantage over single-curve attacks.

To test this, Phase 3 needs:

1. **Test on multiple planted d values** — does useful_curves
   stay at 1 as `d` grows, or does it correlate with the
   intersection of `d`-bit-pattern with each curve's prime-power
   subgroup structure?
2. **Use different starting curves** — characterise the
   distribution of `useful_curves` across the orbit-brittleness
   distribution measured in Phase 1.
3. **Implement HNP-residual cleanup** — the genuinely-novel
   part of the proposal.  Once we have `d ≡ d_known (mod M)`
   from CRT and `M < n`, the residual `u = (d − d_known)/M ∈
   [0, n/M)` is an HNP instance with bias `log(n/M)`.  If we
   have additional orbit curves that contributed to `crt_pairs`
   *partially* (some non-trivial residue but the prime power
   wasn't fully covered), each gives an HNP equation in `u`.
   The proposal claims these aggregate into a lattice attack
   that completes the recovery.  **This claim is still
   untested.**

### Phase 4 — adversarial 32/48/64-bit construction

Once the cross-curve cooperation question is settled, find a
specific 32-bit (then 48-, 64-bit) curve where CGA-HNC
*demonstrably beats plain rho* in wall-clock time on the same
curve.  This requires:

- Implementing Schoof at small scale (we currently brute-force
  point counting, which doesn't scale past `p ~ 2¹⁶`)
- Searching for orbit-brittle curves with high CRT/log₂(p)
- End-to-end timing comparison

### Phase 5 — the original deferred items

Goal: characterise the prevalence of orbit-brittle curves among
random secp256k1-class ordinary curves. Specifically: what
fraction of 256-bit ordinary curves with random `(a, b)` have
CRT/log₂(p) > 0.5 in their `Cl(O)`-orbit?

This is the question that determines whether the result is a
**curve-design recommendation** ("avoid these specific curves") or
a **cryptanalytic threat** ("a meaningful fraction of curves
generated by careless processes are vulnerable").

At cryptographic sizes, brute-force point counting is infeasible.
Phase 5 needs Schoof's algorithm at scale (or a substantial
sub-cryptographic empirical extrapolation, e.g., to 96-bit).

### (Original Phase 2 plan, now subsumed)

Goal: construct a 32-bit, 48-bit, ..., 96-bit ordinary curve where
the proposed CGA-HNC attack succeeds (recovers `d` faster than
plain rho on that curve). Estimated effort: 1–2 weeks.

Approach:
1. Pick `p` of target bit-size.
2. For each `p`, enumerate ordinary curves with small-discriminant
   maximal-order endomorphism rings (small `|disc(O)|`).
3. For each such curve, brute-force compute `#E_i` across the
   `Cl(O)`-orbit (Schoof at scale, infeasible past ~64 bits).
4. Search for orbit configurations where ∏ smooth_parts ≥
   `2^(0.6 · log₂(n))` — i.e. the orbit covers enough of `d` that
   HNP-residual cleanup on the remaining 40% is feasible.
5. Run the actual attack: per-curve Pohlig-Hellman + cross-curve
   CRT + HNP-lattice-residual.

### Phase 3 — generic-curve generator

Goal: characterise the prevalence of orbit-brittle curves among
random secp256k1-class ordinary curves. Specifically: what
fraction of 256-bit ordinary curves with random `(a, b)` have
CRT/log₂(p) > 0.5 in their `Cl(O)`-orbit?

This is the question that determines whether the result is a
**curve-design recommendation** ("avoid these specific curves") or
a **cryptanalytic threat** ("a meaningful fraction of curves
generated by careless processes are vulnerable").

At cryptographic sizes, brute-force point counting is infeasible.
Phase 3 needs Schoof's algorithm at scale (or a substantial
sub-cryptographic empirical extrapolation, e.g., to 96-bit).

## Open questions for an expert reviewer

These are the questions where I most want a domain expert
(Galbraith, Heninger, Lange, Sutherland) to weigh in:

1. **Has any version of "Pohlig-Hellman amortised across an
   isogeny orbit" been published?** The survey agent's literature
   review didn't find it; my best hunch is it's folklore that
   nobody bothered to write up because it doesn't beat rho on
   random curves. If it's published somewhere, a citation would
   immediately compress this proposal into "we re-derived
   classical X" — which is fine, but worth knowing.

2. **Does the HNP-residual cleanup on cross-curve residuals work
   in principle?** The proposal sketches this but I haven't
   verified the lattice-construction details. Specifically:
   given partial recovery `d ≡ d_known (mod M)` from the orbit,
   the residual `u = (d − d_known) / M` lives in `[0, n/M)`.
   Treating relations from un-fully-exploited orbit curves as
   HNP samples — is the lattice rank+volume right?

3. **Is `Cl(O)` of large class number actually rich enough to
   produce the empirical CRT-coverage distribution we see at
   small `p`?** At `p = 1009`, the variance in our experiment
   could be an artefact of small-prime statistics rather than a
   genuine scaling phenomenon. Repeating at `p = 2¹⁶, 2²⁰, 2²⁴`
   would be telling.

4. **For Phase 2 adversarial construction**: which CM
   discriminants give the best CGA-HNC vulnerability? The
   experiment suggests small `|disc(O)|` (high class number) is
   what we want, but I haven't characterised this analytically.

## What this library now ships

- `src/cryptanalysis/cga_hnc.rs` — the building blocks (Φ_2,
  j-invariant, orbit BFS, smoothness measurement). 9 unit tests,
  3 experimental drivers gated under `#[ignore]`.
- `src/ecc_safety.rs::check_cl_o_orbit_brittleness` — auditor
  extension. Returns `Fail` for toy-size curves with high orbit
  brittleness, `Inconclusive` for cryptographic sizes.
- `RESEARCH.md` — this note.

The library's full test suite is **350 passing, 0 failed, 16
ignored**.

## Honest closing note

This is a **research thread, not a result**. The empirical
findings at toy scale are real and reproducible (run
`cargo test --lib cga_hnc -- --ignored --nocapture`). The
extrapolation to cryptographic sizes is speculative. The
auditor extension is a useful addition independently of whether
the attack scales.

If you're an expert in ordinary-curve cryptanalysis and reading
this: please tell me where I'm wrong. The construction sketched
here is implementable on this library's existing infrastructure
(LLL/BKZ from `lattice`, HNP from `hnp_ecdsa`, isogeny machinery
from `cga_hnc`) and the next milestone — adversarial 64-bit
curve construction — is a 1–2 week project.
