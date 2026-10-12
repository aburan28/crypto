# Transfer-aware cross-representation index calculus: what a transfer can and cannot change

**Date:** 2026-10-10.  **Status:** structural note with one measured toy
experiment (`research/transfer_aware_ic_20261010/`).  **Result class:**
boundary; `S`, end-to-end cost and speedup are **unset**.  Nothing here is
a claim about ECC2K-130, P-256 or any deployed curve beyond arithmetic that
is labelled as such.
**Conductor task:** T-251.
**Companions:** `RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md` (isogenous-curve
PDP probe), `RESEARCH_ISOGENY_CONDUCTOR_GAP_20261007.md` and
`EXPERIMENTS_ISOGENY_CONDUCTOR_GAP_20261007.md` (levels of the ECC2K-130
class), `RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md` (§A1, the multiplier
closure), `research/isogeny_class_difficulty_20261008/` (protocols I-1 to
I-5), `research/jmv_self_reducibility_20261008/`,
`research/graph_cycle_ic_20261007/`, `research/cheap_closures_20261007/`,
`research/cover_route_prime_binary_20261008/`,
`research/hyperelliptic_cover_infrastructure_20261006/`.

## 0. The question

The request was to look for relation-yield gains from (a) horizontal
isogenies, (b) curves at lower volcano levels, (c) factor bases spread
over several isogenous curves, (d) loops of isogenies, and more generally
(e) a *transfer-aware* index calculus that searches over sequences of
isogenies, field representations and Jacobian transfers for the
representation in which the DLP instance is cheapest.

Most of (a)–(d) has been measured in this repository already and the
pointers are collected in §2.  What was missing is the statement that
ties the measurements together: **a transfer that is a group homomorphism
cannot change decomposition counts**, so every one of (a)–(d) reduces to
three curve-specific quantities that are already enumerated.  §1 states
and proves that, §2 maps every transfer type onto it, §3 gives the scoring
function a transfer-aware pipeline would have to use and shows which
destinations can score at all, §4 reports the one lever that was still
unmeasured (the multiplier closure, which is what "loops of isogenies"
means at the crater), and §5 says what remains.

## 1. Transport conservation

Let `φ: E → E′` be an `F_q`-rational isogeny (or any group homomorphism of
the rational points), `R ∈ E(F_q)`, and `F′ ⊂ E′(F_q)` a factor base on the
codomain.  Write `F := φ⁻¹(F′) ∩ E(F_q)` for the pulled-back base and
`κ := ker φ ∩ E(F_q)` for the rational part of the kernel.

**Lemma (conservation).**  The `m`-sum decompositions of `φ(R)` over `F′`
and of `R` over `F` are in bijection up to the action of `κ^{m−1}`:

    #{(Q_1,…,Q_m) ∈ F′^m : Σ Q_i = φ(R)}  ·  ∏_i |φ⁻¹(Q_i) ∩ E(F_q)|
        =  |κ|  ·  #{(P_1,…,P_m) ∈ F^m : Σ P_i = R}.

*Proof.*  If `Σ P_i = R` with `P_i ∈ F` then `Σ φ(P_i) = φ(R)` with
`φ(P_i) ∈ F′`.  Conversely, given `Q_i ∈ F′` with `Σ Q_i = φ(R)`, choose
rational preimages `P_i ∈ φ⁻¹(Q_i) ∩ E(F_q)`; then `Σ P_i − R ∈ κ`, and
exactly `|κ|^{m−1}` of the `∏ |φ⁻¹(Q_i) ∩ E(F_q)|` choices land on `R`
itself (translate `P_1` by the kernel element).  ∎

Three consequences, each the closing statement of a thread that was run
empirically:

- **C1 (transported and pulled-back bases are the same experiment).**
  Decomposing on `E′` over `φ(F)` or on `E` over `φ⁻¹(F′)` gives the same
  counts as the native experiment, up to the `|κ|^{m−1}` symmetry.  This is
  why `ic-conductor-leads/LEADS.md` found transported bases "closed", why
  Wave B item 9 of `research/isogeny_volcano_ic_20260924/` demands
  "target-wise decomposition-count conservation" (it is a theorem, not a
  gate), and why the E1 panel of `EXPERIMENTS_ISOGENY_CONDUCTOR_GAP_20261007.md`
  measured level-blind trials per relation.
- **C2 (multi-curve factor bases glued by isogenies are single-curve
  factor bases).**  A base `F_1 ∪ φ⁻¹(F_2) ∪ …` spread over several
  isogenous curves and glued by `log φ(P) = log P` is, on any one member,
  the union of `F_1` with the rational pullbacks of the others.  Its yield
  is that of a base of the same total size; its decomposition problem is
  the mixed one in which some abscissae lie in a subspace and others are
  images under rational functions of degree `deg φ`.  For `deg φ = 2` on
  a binary curve this *is* the `T₂` symmetrisation already in use
  (`x ↦ x + 1/x`); for `deg φ ≥ 3` the kernel is not rational on the
  curves in question, so nothing folds and the mixed system only carries
  higher-degree membership polynomials.  There is no multi-curve yield
  effect to look for beyond the symmetry factor.
- **C3 (sequences compose; loops are multipliers).**  A sequence of
  isogenies and isomorphisms composes to one homomorphism, so the lemma
  applies to the whole route.  A loop `E → … → E` is an endomorphism
  `α ∈ End(E)`, acting on the prime-order subgroup as the scalar `α(μ)`;
  its only use is as a *known-log multiplier* `P ↦ α(P)`.  Loops on the
  floor of the ECC2K-130 class are elements of `Z + 263·O_K`, scalar on
  logs and of degree `≥ 121,046` (`research/ecc2k130_endo_ring_263_20260925/`);
  loops at the crater are `Z[τ]`, which is §4.

**What is not conserved.**  Exactly three things vary between the source
and destination of a homomorphic transfer, and they are the complete list
of model-dependent inputs to an index calculus:

1. the rational kernel symmetry `|κ|^{m−1}`, bounded by the rational
   torsion of the class (`Z/4` for ECC2K-130 and its toy ladder, which is
   the `T₂`/`T₄` folding already used; trivial for prime-order curves such
   as P-256), and shown in `research/cheap_closures_20261007/` to be worth
   at most 2 bits per summand on the Koblitz model;
2. the multiplier set `End(E′)` with its evaluation cost (the crater has
   `Z[τ]` with free `τ`; every other level of the ECC2K-130 class has
   nothing cheap; a prime-field curve of large discriminant has nothing
   but scalars);
3. the factor-base family's on-curve fraction `|F|` and the algebraic
   form of the decomposition polynomial on the destination model.  For
   the binary Semaev `S₃` the curve enters only through the constant `b`,
   so the Macaulay rank profile is a class invariant
   (`RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md` §2, identical ranks across every
   class member measured) and `|F|` moves yield by `≤ 1.3×`; for the
   prime-field enumeration PDP the model does not enter at all.

So "relation yield across horizontal isogenies" and "across lower volcano
levels" are, for a fixed family, the same number up to item 3, and the
measurements agree: `±8 %` in `|F|`, level-blind trials per relation,
identical first-fall degree.  A search over the class for the member with
the largest `|F|` is the only yield lever the class offers, and it is
capped by that `1.3×` while costing the `2n` crater folding whenever the
search leaves `K_0`.

## 2. Transfer catalogue

Every representation change the request names, what it can change under
§1, and where the repository has measured or proved it.

| transfer | homomorphic? | what can change | repository verdict |
|:--|:--|:--|:--|
| horizontal `ℓ`-isogeny, same level | yes | item 3 only (`|F|`, constants of `S₃`) | PDP probe: ranks identical, yield `±8 %` via `|F|`; `K_0` has `h = 1`, so its horizontal edges are loops |
| vertical isogeny to a lower level | yes | items 2 and 3; loses `τ` below the crater | E1 panel at T23: trials/relation 2.41–2.52 crater vs 2.50–2.78 floor; crater advantage is the `2n` folding and nothing else |
| loop of isogenies | yes | item 2: known-log multiplier | `graph_cycle_ic` (negative as an outer algorithm); §4 here as a multiplier |
| field representation (ONB, polynomial, Hamming-weight, cyclotomic-sparse) | identity on the group | changes the *factor-base family*, not the transfer; §1 does not apply, the family's `(γ, C)` does | the subject of the landed rungs; A2 of the directions note is the one untried family |
| model change on the same curve (Edwards, Montgomery, λ-coordinates) | isomorphism | item 1 through the torsion symmetrisation it exposes | `cheap_closures`: pullbacks by `1 ± τ` and `[2]`, `≤ 2` bits |
| base change `E(F_q) ⊂ E(F_{q^k})` | injective | factor base must grow to `≈ q` to use the extension structure | cost `≥ q` against rho's `√r ≈ √q` for the embedded subgroup; never competitive for `k ≥ 2` |
| gluing into a Jacobian over the same field (`(ℓ,ℓ)`-isogeny `E × E′ → Jac C`, conorm `E → Jac C`) | injective on the subgroup | factor base `C(F_q)` of size `q` with cheap decomposition | cost `≥ q^{2−2/g} ≥ q` against `√q`; loses for every genus |
| cover over a subfield (GHS, Galbraith–Hess–Smart isogeny extension) | injective on the subgroup when the subgroup survives the conorm | the one transfer that changes the exponent, because the destination field is smaller | closed for every deployed curve by field structure (`cover_route_prime_binary_20261008`: no proper subfield over `F_p` or prime-degree `F_{2^n}`; the only standardised targets are RFC 2409's two EC2N groups) |
| pairing to `F_{q^k}^×` (MOV/Frey–Rück) | injective | a different index calculus altogether | embedding degree is a class invariant; large for every registered curve |
| twist / `F_{q²}` structure (ECC2K-130: `E[263] ⊂ E(F_{q²})`, `263² \| #K_1(F_q)`) | not a map of the `r`-subgroup | nothing: `⟨G⟩` is trivial in `E/263E` | `RESEARCH_ISOGENY_CONDUCTOR_GAP_20261007.md` V8 |

The cover row is the only one whose destination has a *different*
factor-base family with a better `(γ, C)` after the group-size change,
which is why it is the only transfer with a historical exponent gain, and
why its absence on prime fields and prime-degree binary fields is a
statement about fields, not about search effort.

## 3. The transfer-aware scoring function, and why the search space is empty for deployed curves

A transfer-aware pipeline would pick, among reachable destinations `D`
with group `G_D ⊇ ⟨G⟩`, factor-base family `B_D`, relation yield `γ_D`,
decomposition cost `C_D` and route cost `T_D`, the one minimising

    cost(D) = T_D + |B_D| · C_D / γ_D + LA(|B_D|) + online(D),

and compare it with matched rho on the source, `R ≈ √(πr/4)` folded by
the available automorphisms.  Under §1:

- if `D` is isogenous to the source (any sequence of isogenies,
  isomorphisms and field representations), `γ_D`, `C_D` and `|B_D|` are
  those of the source family up to items 1–3, and `T_D ≥ 0`; the search
  can only lose, except by the bounded `|F|` screen and by finding a
  member with a larger multiplier ring or more rational torsion, both of
  which are class invariants readable from `(q, #E, D_π)` without any
  walk;
- if `D` enlarges the group (`E(F_{q^k})`, a Jacobian over `F_q`, a
  pairing target), `|B_D| ≥ q` for every family with a cheap
  decomposition, so `cost(D) ≥ q > √r` whenever `r ≈ q`;
- if `D` shrinks the field (a cover over a subfield), the exponent can
  drop, and there is no subfield to shrink to for `F_p` or prime-degree
  `F_{2^n}`.

So for every curve in the registered ledger the transfer-aware search has
no destination that can beat its source by more than a constant, and the
constant is bounded by the items in §1.  This is the same conclusion the
cover note reached for covers and the isogeny-class note reached for
generic attacks, now for all homomorphic transfers at once.  It is not a
proof that no non-generic method exists; it is a proof that *representation
search over homomorphic transfers* is not where one would come from.  The
two things in this family that remain mathematically open are the vertical
step at a large prime conductor (`jmv_self_reducibility_20261008`, which
would not change any IC count even if solved) and explicit small-degree
covers of `E/F_{q³}` classes (`cover_route_prime_binary_20261008`, which
does not concern deployed fields).

One observation that is useful rather than negative: a factor-base log
table computed on `E` is a log table for every curve reachable from `E` by
a cheap isogeny, because `log_{φ(G)} φ(P) = log_G P`.  For ECC2K-130 that
is the crater plus the 262 constructible floor curves; for a walked
prime-field neighbourhood it is every visited node.  Precomputation in
this pipeline is therefore a *class-level* asset, exactly as a rho DP
table is (`EXPERIMENTS_ISOGENY_CONDUCTOR_GAP_20261007.md` §0 item 2), and
the amortised-cost accounting of `PLAN_IC_ACCOUNTING_FIXES_20261007.md`
should count the targets it serves across the class, not per curve.

## 4. The one unmeasured lever: multiplier closure

At the crater, loops of isogenies are `Z[τ]`.  The landed compact-orbit
producer uses the multipliers `{±τ^k}`, `2n` per column; §A1 of the
directions note proposes adding the small-norm elements of `Z[τ]`.  Two
derivations, stated before the experiment:

1. **In the 4-sum pair-table pipeline the Frobenius orbit's role in the
   yield model is "cheap known-log multiples", not Frobenius.**  With `K`
   columns and `c` multipliers per column (any multipliers), the summand
   set has `L = 2n·K·c` points, a probe hits with probability `≈ L²/(2r)`,
   relations needed are `K`, and the τ-canonicalised table has `≈ L²/(4n)`
   entries.  Integer multiples `{1, …, c}` and `Z[τ]` elements of norm up
   to about `c` are interchangeable in every one of these quantities;
   `Z[τ]` buys only free evaluation of `α(P)`, which is a precompute cost
   of `K·c` point operations and negligible beside `K·2r/L²` probes.
2. **At equal table memory the closure cannot change the per-target probe
   count.**  `L` is fixed by memory, so `2r/L²` is fixed; the closure
   trades columns for multipliers, `K ↦ K/c`, which divides relations and
   therefore total precompute probes by `c` and the linear algebra by
   `c²`.  The trade stops when the multiples of a single column start to
   coincide additively (`aP + bP = (a+b)P` and its `Z[τ]` analogue): the
   fraction of same-column pairs is about `1/K`, so the loss appears at
   small `K`, not at small `c`.  The limit `K = 1` is a four-list
   meet-in-the-middle on multiples of `G`, i.e. a BSGS with a Frobenius
   table saving, which is the generic time–memory trade the whole
   pipeline sits on.

The frozen protocol and instrument are in
`research/transfer_aware_ic_20261010/` (`PROTOCOL.md`, `multipliers.rs`,
`score.py`).  Three cells, `n = 23, 37, 41` (`r` prime, `≈ 2²¹`, `2^{27.8}`,
`2^{39}`), each with the `frob` arm at `K = Kc` and `tau`/`int` arms at four
divisors of `Kc` down to `K = 1`, all at the same `L` and the same table.
Every arm solves the column logs against hidden scalars and recovers one
hidden target.

### 4.1 Results

Full tables, scoring and files: `research/transfer_aware_ic_20261010/README.md`.
Cells C23, C37 and the amended C37b are complete; C41 and its amended C41b
were relaunched after the first launch was killed with its session and
are reported from their partial logs until the JSON lands.

| cell | `2r/L²` | `frob` probes | `tau`/`int` at `K ≥ 3`, vs `frob` | `tau`/`int` at `K = 1`, vs `frob` | canonical at `K = 1`, vs `frob` |
|:--|--:|--:|:--|:--|:--|
| C23 (`n=23`, `Kc=12`) | 13.8 | 15.1 | 0.91–1.24 | 1.73 / 1.73 | 0.58 / 0.62 |
| C37 (`n=37`, `Kc=36`) | 65.0 | 63.9 | 1.05–1.26 | 1.71 / 1.49 | 0.64 / 0.72 |
| C37b (amended, 128 hits) | 65.0 | 68.0 | 0.88–1.18 | 2.15 / 1.49 | 0.64 / 0.71 |
| C41 (partial, 48 hits) | 31,543 | 35,861 | 0.81–1.01 (`K = 36, 12`) | pending | pending |
| C41b (partial, 160 hits) | 31,543 | 29,143 | 1.19 (`K = 36`) | pending | pending |

- M1 (yield law `2r/L²`): pass in every cell, within 10 %.
- M2 (multiplier arms within `[0.7, 1.5]×` `frob` at `K ≥ 3`): pass in
  every arm.
- M3 (`tau/int ∈ [0.6, 1.6]`): pass; the decisive outcome `tau ≤ 0.67 × int`
  that would have given `Z[τ]` a structural advantage did not occur.
- M4 (coincidence loss at `K = 1`): the table shrinks to `0.58–0.72×` as
  predicted; probes rise `1.49–2.15×`, which misses the frozen `1.5×` bar
  on two `int` arms by 0.01.  Recorded as a failure of the band, with the
  effect present.
- M5/M6: every `tau`/`int` arm reached full rank in `K` or `K + 1` rows
  with no degenerate rows and verified its logs and hidden target.  The
  `frob` arms exposed two protocol errors (rank excess for sparse rows not
  allowed for; C41 registered with fewer hits than columns) and one tail
  event (a column untouched in 64 rows at C37); amended cells with
  `hits ≥ 2K` are reported beside the registered ones.
- M7 (total probes to full rank scale with `K`): pass on 9 of 11 scorable
  arms.

### 4.2 Reading, and the arithmetic at the landed rungs

1. **The Frobenius orbit's yield role is replicable by integers.**  Every
   `int` arm sits on top of its `tau` twin and both sit on `frob`.  In the
   pair-table pipeline, "Koblitz structure" in the *yield* is the
   availability of `2n` cheap known-log multiples per column, and the
   genuinely Frobenius-specific saving is elsewhere: the `2n`
   canonicalisation of the table (the `canonical` column is `L²/(4n)`
   whatever the multipliers, because the multipliers commute with `τ`)
   and the free evaluation of `τ^k P`.
2. **The closure is a precompute trade, not a yield lever.**  Per-target
   probes do not move with `c` at fixed `L` (M2); relations and linear
   algebra move by `c` and `c²` (M7).  At the landed `n = 73` shape
   (`K = 600`, `2nK = 87,600` summands, `5.7·10⁷` probes per relation,
   18,398 s of rank wall), `c = 10` at equal table memory would mean
   `K = 60`, the same `5.7·10⁷` probes per relation, ten times fewer
   relations, a `60 × 60` instead of `600 × 600` solve, and an online
   figure that is unchanged.  Arithmetic, not a measurement.
3. **The trade has a floor set by coincidences, and it is far from the
   landed rungs.**  The loss appears at `K ≤ 2` (`1.2–1.3×` at `K = 2`,
   `1.5–2.2×` at `K = 1`) and is absent at `K ≥ 3`; its mechanism is the
   `≈ 1/K` fraction of same-column pairs whose sums collapse
   (`aP + bP = (a + b)P`).  At `K = 60` it is negligible.  The `K = 1`
   arm is the four-list BSGS on multiples of `G` with a Frobenius table
   saving, which is the generic time–memory trade the pipeline is built
   on; the measured `1.5–2.2×` is the price of that limit.
4. **What this does and does not say about ECC2K-130.**  It says that
   "loops of isogenies at the crater" are a known-log multiplier set, that
   such a set changes the precompute size and not the exponent, and that
   its benefit is bounded by the `K → 1` limit.  It does not move any
   ledger row: no `S`, no end-to-end cost, no online figure changes.

## 5. What remains, by class

- **Closed as yield levers** (this note, with the pointers in §2):
  transported, pulled-back and multi-curve factor bases; horizontal and
  vertical neighbours for a fixed family; loops beyond their use as
  multipliers; model changes beyond the torsion symmetrisation;
  extension-field and same-field Jacobian embeddings; covers on deployed
  fields.
- **Constant-factor engineering option, measured at toy size:** the
  multiplier closure as a precompute trade (§4).  Whether to implement it
  in the compact-orbit producer is a question of whether precompute size
  or linear-algebra width is the binding constraint at the next rung; the
  online figure the ledger reports does not move.
- **Preregistered, cheap, expected null:** protocol I-4
  (`research/isogeny_class_difficulty_20261008/PROTOCOL-I4-ic-yield-across-class.md`)
  with the native `ic` pipeline would confirm §1 item 3 on prime-field
  toys, where the enumeration PDP predicts exact flatness.
- **Open mathematics that would not change an IC count:** the vertical
  step at a large prime conductor; explicit covers of `F_{q³}` classes.
- **Not a transfer question at all, and the actual levers of the Koblitz
  ledger:** higher-arity extraction (`m = 5, 6`), the implicit `S₃` index,
  the equivariant F4 (`research/g2_equivariant_f4_20261008/`), the
  cyclotomic-sparse family (A2).  Those are where the next rung comes
  from, and none of them is reached by changing the curve.
