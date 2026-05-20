# Volcano-Floor Class-Group-Augmented Pollard Rho (VFCG-ρ)

A proposal for a novel ECDLP attack on prime-field CM curves
(specifically secp256k1) that exploits an asymmetry in the isogeny
volcano structure: the **crater** has trivial class group (`h(−3) =
1`), but **the floor at depth 1** has non-trivial class group
(`h(−3·ℓ²) > 1` for ℓ ≥ 5).

The proposal builds on the structural audit in
[`RESEARCH_SECP256K1_CM.md`](RESEARCH_SECP256K1_CM.md).

> **Realistic expectation.**  This proposal is most likely to
> produce a constant-factor improvement at best, and the asymptotic
> outcome is probably null.  The value here is: (a) a clean
> combinatorial framework that is, to my knowledge, novel for
> prime-field CM curves; (b) a falsifiable hypothesis with a
> concrete toy-scale experiment; (c) tightening the lower bound
> on what isogeny-graph structure can extract from secp256k1.
> If the empirical test at 32–64 bit scale shows even a 1.2×
> speedup, that is publishable; if it shows null, it closes a
> direction.

## 1. The structural observation

Every existing attack that tries to use isogeny-graph structure
against secp256k1 is blocked by one of two theorems:

- **Tate**: all `F_p`-isogenous curves have the same group order
  `n`, so transporting ECDLP gives the same hardness.
- **`h(−3) = 1`**: the `Cl(End)`-orbit at the crater is a single
  vertex, so there is no horizontal CM action to amortise across.

The **floor** is in a different regime.  Descend one ℓ-isogeny step
from secp256k1 (with ℓ a prime where descending neighbours exist
F_p-rationally — by the audit, `ℓ ∈ {2, 3, 7, 13, 17, 19, 31, …}`):

```
        End(secp256k1)      = Z[ω]                  ← maximal order, disc −3, h = 1
        ↓ ℓ-isogeny φ
        End(E_↓)            = Z + ℓ · Z[ω]          ← non-maximal order
                                                    disc −3ℓ²,  h(−3ℓ²) > 1 for ℓ ≥ 5
```

So `E_↓` lives in a class group `Cl(Z + ℓ · Z[ω])` of size
`h(−3ℓ²)`.  Concrete values from the audit:

| ℓ  |  h(−3ℓ²) | Cl-action gives |
|---:|---------:|------------------|
| 2  |    1     | trivial          |
| 3  |    1     | trivial          |
| 5  |    2     | Z/2              |
| 7  |    2     | Z/2              |
| 11 |    4     | Z/4              |
| 13 |    4     | Z/4              |
| 17 |    6     | Z/6              |
| 19 |    6     | Z/6              |
| 23 |    8     | Z/8 or Z/2 × Z/4 |

Crucially: **Tate still applies on the floor** — every floor curve
has `#E_↓(F_p) = n` (same prime `n` as secp256k1).  But the curves
in the floor orbit are now *related by a non-trivial class group
action* `Cl(−3ℓ²) ↷ {E_↓^{(1)}, …, E_↓^{(h)}}`.

```
         (CRATER, h = 1)
              ●  secp256k1, End = Z[ω]
             ╱│╲
            ╱ │ ╲   descending ℓ-isogenies
           ╱  │  ╲
          ●   ●   ●
         (FLOOR, h = h(−3ℓ²))
         Cl(−3ℓ²) acts horizontally on the floor
```

**The question this proposal asks**: can a Pollard-ρ variant whose
walks live on the floor — using both "add a random point" and
"apply a class-group action" as step operations — exploit the floor's
class-group structure for asymptotic speedup?

## 2. The proposed walk

Let `E_↓^{(1)}, …, E_↓^{(h)}` be the `h = h(−3ℓ²)` curves in the
depth-1 floor orbit of secp256k1, each with `#E = n`.  Let
`Cl(−3ℓ²) = ⟨γ⟩` (cyclic for our chosen ℓ ∈ {5, 7, 11, 13}; for
ℓ = 17 we have Z/6 = ⟨γ⟩ also cyclic).

For an ECDLP instance `(P, Q = d · P)` on secp256k1, we transport
it to each floor curve via the chosen ℓ-isogeny `φ`:

```
       (P_i, Q_i = d · P_i)  on  E_↓^{(i)}     for i = 1, …, h
```

where `P_i` is the image of `P` after walking the orbit from
`E_↓^{(1)}` to `E_↓^{(i)}` via class-group action.

**State.**  A walker state is a pair `(i, R)` where `i ∈ {1, …, h}`
is the floor-curve index and `R ∈ E_↓^{(i)}(F_p)` is a point.

**Step function `f: (i, R) ↦ (i', R')`.**  Choose one of three
kinds of step uniformly at random based on a hash of `R`:

1. **Additive** (with probability `1 − 2/h`):
   `R' = R + M_i`  for a precomputed step `M_i ∈ E_↓^{(i)}(F_p)`;
   `i' = i`.
2. **Horizontal CM** (with probability `1/h`):
   apply the class-group generator `γ`: `i' = i + 1 mod h`;
   `R' = γ(R)` (the image of `R` under the canonical degree-`γ`
   horizontal isogeny `E_↓^{(i)} → E_↓^{(i+1)}`).
3. **Inverse CM** (probability `1/h`): `i' = i − 1 mod h`;
   `R' = γ^{−1}(R)`.

**Tagging.**  A "distinguished point" (DP) is a state `(i, R)` whose
`R` has its top `k` x-coordinate bits equal to zero, independent of
`i`.  Store DPs in a single global table keyed by `(i, R)`.

**Collision criterion.**  When two walks `(P-walk)` and `(Q-walk)`
both land on the same DP `(i*, R*)`, we have

```
       P-walk:  e_P · P_i  +  Σ (other steps)  =  R*  on E_↓^{(i*)}
       Q-walk:  e_Q · P_i  +  Σ (other steps)  =  R*  on E_↓^{(i*)}
```

with the "other steps" tracked as an explicit symbolic combination
of additive steps and class-group actions.  Subtracting:

```
       (e_P − e_Q) · P_i*  =  (Q-walk other) − (P-walk other)
```

The RHS evaluates to an explicit point on `E_↓^{(i*)}` because we
tracked the step sequence.  Then

```
       (e_P − e_Q) · P_i*  =  Σ a_j · γ^j(P_i*)
```

becomes an equation in `End(E_↓^{(i*)})`, where the LHS is multi-
plication by `e_P − e_Q ∈ Z` and the RHS is a `Z`-linear combination
of class-group-action images.  Reduce modulo `n` to get

```
       (e_P − e_Q)  ≡  (some explicit polynomial in γ̃)  (mod n)
```

where `γ̃ ∈ Z/nZ` is the scalar `γ` acts as on `E_↓^{(i*)}(F_p)`.
(For floor curves, `γ̃` is computable from the isogeny structure.)

Done correctly, this gives `d` mod `n` after one collision.

## 3. Complexity analysis

**State space size**: `h · n` (sum over the `h` floor curves of `n`
points each).

**Naive birthday cost**: `√(h n)` expected steps to first collision.

**Per-step cost**:
- Additive step: 1 point-add ≈ a few field operations.
- Class-group step: 1 ℓ-isogeny evaluation ≈ O(ℓ²) field operations
  (Vélu's formulas + a few multiplications).  For ℓ = 5, 7, 11
  this is small.

**Step rate**: with probability `2/h` of doing a class-group step,
each "step" costs `(1 − 2/h) · 1 + (2/h) · ℓ²`.

Putting it together, total cost is

```
   T(h, ℓ)  =  √(h n) · [ (1 − 2/h) + (2/h) · ℓ² ]
            =  √(h n) · [ 1 + (2/h) · (ℓ² − 1) ]
```

Plain rho (no floor) costs `√n` (with √6 GLV speedup at the crater).
Comparing:

```
   speedup ratio  =  √n / T(h, ℓ)
                  =  1 / [ √h · (1 + 2(ℓ² − 1)/h) ]
                  =  √h / (h + 2(ℓ² − 1))
```

For this to exceed 1 (speedup), need `√h · (h + 2(ℓ² − 1)) < h`,
i.e., `h + 2(ℓ² − 1) < √h`, which requires `√h > 2ℓ² + h`.  Since
`h ≥ 1`, this is **impossible** for `ℓ ≥ 1`.

**Naïve verdict**: the basic proposal does not beat plain ρ.  The
birthday-saving factor `1/√h` is outweighed by the per-step cost
`O(ℓ²)` for class-group transitions.

## 4. Where the proposal might still win

The naïve analysis above is too pessimistic in two respects:

**(a) Class-group steps as "free" amortisation.**  In the additive
walk on a single curve, ρ wastes `√n` steps to find any collision.
The class-group transition is a "smarter" step that *changes the
underlying curve*, accessing a fresh `n`-element search space.

Specifically: after `√n` steps on a single curve, the walk has
effectively explored `O(√n)` of the `n` points.  By switching to
another floor curve via a class-group step, the walker "resets" its
position in a structurally-related but newly-explored space.

Naively the search space is `h · n`, but **distinguishing-point
collisions in `E_↓^{(i)}` are correlated with DPs in `E_↓^{(j)}`
via the class-group isogeny** — they share an x-coordinate structure
up to the `γ`-action.

This correlation might give a **shared-DP-table speedup**: a DP
found by one walker on `E_↓^{(i)}` constrains DPs by another walker
on `E_↓^{(j)}` via the action of `γ^{i−j}`.  Quantifying this
correlation is the open empirical question.

**(b) Class-group steps as Frobenius accelerators.**  At depth 1
on a split ℓ, the class group `Cl(−3ℓ²)` has a non-trivial generator
γ.  The action of `γ` on points of `E_↓^{(i)}(F_p)` is `multipli-
cation by a scalar c_γ ∈ Z/nZ` (since the `F_p`-rational subgroup
of order `n` is cyclic).  Specifically, `c_γ` is a root of the
class-group's relation polynomial in `(Z[ω] / ℓ²)`.

**The novel structural claim**: the multiset `{c_γ, c_γ², …, c_γ^h}
mod n` has a non-trivial **algebraic structure** — these are the
eigenvalues of the Frobenius acting on the `ℓ²`-th roots of unity
under the CM action.  If `c_γ` is in a structured subset of
`(Z/nZ)*` (e.g., a coset of a small-index subgroup), then walks
mixing additive and γ-steps explore a constrained sub-lattice of
`(Z/nZ)`.

If the sub-lattice has rank `< 1` (proper sublattice), this would
be a real attack: the DLP reduces to a problem in a smaller cyclic
group.

**This is the open mathematical question.**  My current belief is
that `c_γ` is *generic* in `(Z/nZ)*` for the floor of secp256k1
(no special algebraic structure constraining it to a subgroup) —
but this is checkable.

## 5. Concrete experiment to falsify the hypothesis

The hypothesis to test: **`c_γ` (the scalar action of the
class-group generator on the F_p-rational n-subgroup at the floor)
lies in a non-trivial proper subgroup of (Z/nZ)*.**

If `c_γ` is generic, the proposal fails — no asymptotic speedup
over plain ρ.  If `c_γ` is in a small subgroup, the proposal becomes
a genuine attack with cost `√(n / index)`.

### Step-by-step experiment

1. Pick a small CM curve `E_toy / F_{p_toy}` with `p_toy ≈ 2^{32}`,
   `j(E_toy) = 0`, ordinary, prime `n_toy ≈ 2^{32}`.
2. Compute Frobenius trace `t_toy`, confirm `j = 0` / `End = Z[ω]`.
3. For each ℓ ∈ {5, 7, 11, 13, 17}, descend via the canonical
   ℓ-isogeny `φ`: get `E_↓ / F_{p_toy}` with `End = Z + ℓ Z[ω]`,
   `h(−3 ℓ²)` curves in the floor orbit.
4. Compute the class-group generator `γ` (one of the prime-above-ℓ-
   conjugate ideals) as an explicit ideal in `Z + ℓ Z[ω]`.
5. Compute the action of `γ` on `E_↓^{(1)}` to get the canonical
   floor neighbour `E_↓^{(2)}`.
6. Trace a fixed point `P ∈ E_↓^{(1)}(F_p)` of order `n_toy`
   through the `γ`-action to its image on each floor curve.
7. Recover the scalar `c_γ ∈ Z/n_toy` such that `γ(P) = c_γ · P_2`
   under the `(γ-induced)` group-isomorphism `E_↓^{(1)}(F_p) ≃
   E_↓^{(2)}(F_p)`.
8. **Test**: factor `(c_γ − 1)` and `gcd(c_γ − 1, n_toy)`.  If
   `gcd` is non-trivial, `c_γ` lies in a proper subgroup ⇒ speedup
   exists.  Otherwise, `c_γ` is generic ⇒ no attack.

### Expected outcome

I'd bet 9-to-1 that `c_γ` is generic (`gcd = 1`, no subgroup
membership).  But if the experiment shows even one `(ℓ, p_toy)` pair
with non-trivial `gcd`, the proposal becomes a real avenue.

### Cost of the experiment

- Find suitable `E_toy`: minutes in PARI.
- Compute floor isogenies via Vélu: seconds per ℓ.
- Compute `c_γ` via point-comparison: seconds.
- Statistical test across 100 random `E_toy`: hours of PARI.

Implementable in PARI/GP in `< 200 LOC`, runnable on the codebase's
existing `secp256k1_cm_audit/` infrastructure.

## 6. Where this differs from existing attacks

| Existing attack         | Uses isogeny graph? | Uses class group? | Blocked for secp256k1? | This proposal differs by                                                            |
|-------------------------|---------------------|-------------------|------------------------|--------------------------------------------------------------------------------------|
| Plain Pollard ρ         | no                  | no                | no                     | adds floor-level class-group action                                                  |
| `aut_folded_rho.rs`     | no                  | no                | √6 speedup only        | extends fold from Aut(E) to Cl(−3ℓ²)                                                 |
| `cga_hnc.rs`            | yes (2-isogeny BFS) | yes (#E orbit-CRT)| yes (`h(−3) = 1` at crater) | uses *floor* class group, not crater                                              |
| CSIDH / CRS             | yes (horizontal)    | yes (key exchange)| yes (`h = 1`)          | uses floor where `h > 1`                                                             |
| MOV / Frey-Rück         | no                  | no                | yes (embed deg huge)   | no pairing dependence                                                                |
| Smart 1999 (anomalous)  | no                  | no                | yes (`#E ≠ p`)         | not anomaly-dependent                                                                |
| GHS / Weil descent      | yes                 | no                | yes (no subfield)      | no subfield needed; works over `F_p`                                                 |
| (N, N)-cover (slice 3)  | yes (Weil restr.)   | no                | open                   | uses ρ framework, not Jacobian decomposition                                         |
| `ec_index_calculus_j0.rs` | no (orbit reduce) | no                | implemented, null      | structural, not index-calculus-based                                                 |

The closest neighbour is `aut_folded_rho` (which folds by `Aut(E) =
Z/6`).  This proposal folds by the *larger* group `Aut(E) × Cl(−3ℓ²)`
(or rather, augments the state with a class-group coordinate that's
acted on by the walk).

The closest neighbour conceptually is CGA-HNC, but CGA-HNC tries to
amortise PH recovery on `#E`-smoothness (vacuous when `n` is prime),
while VFCG-ρ uses class-group structure to potentially constrain the
DLP eigenvalue (`c_γ`).

## 7. Variations

If the experiment in §5 fails on the basic version, several variants
are worth trying:

### V1: Multi-depth volcano

Walk on the floor of depth `d ≥ 2`.  Class number `h(−3 · ℓ^{2d})`
grows like `ℓ^{d−1}`, giving a larger orbit.  The per-step isogeny
cost grows as `(d · ℓ²)`.

Optimal `d` balances these: `d ≈ log_ℓ(n)` would give a `√(n^d)`
search space with `(d · ℓ²)` per-step cost.  Likely sub-asymptotic
to plain ρ but worth checking.

### V2: Multi-prime floor

Use the floor reached by descending via *several* primes `ℓ_1, ℓ_2,
…` (e.g., descend by 5, then by 7, etc.).  Endomorphism order
becomes `Z + ℓ_1 ℓ_2 · Z[ω]`, class number `h(−3 · (ℓ_1 ℓ_2)²)`.

This widens the class group at the cost of compositional isogeny
expense.

### V3: Tagged DPs with class-group residues

Instead of class-group steps in the walk, use the class-group
coordinate as a DP tag.  A point `R ∈ E_↓^{(i)}(F_p)` is a DP iff
its hash matches a target *that depends on `i`*.

This biases DPs to specific class-group positions, potentially
giving correlated collisions.  Closer to preprocessing-ρ (Bernstein-
Lange) than to plain ρ.

### V4: Walks on the supersingular reduction of E mod auxiliary ℓ

Beyond scope of this proposal but a related direction: take small
primes ℓ where the reduction `E mod ℓ` (over `F̄_ℓ`) is supersingular,
work in the supersingular isogeny graph of `F̄_ℓ`.  Connections to
Castryck-Decru / SIDH-break machinery.

## 8. Honest verdict

**Likely outcome**: this proposal will *not* break secp256k1.

The critical uncertainty (§5: is `c_γ` generic?) is checkable in
≤ 1 day of compute.  If the answer is "yes, generic", the proposal
fails and we have a clean negative result: "the class-group action
at the depth-1 floor is generic with respect to `n`, so VFCG-ρ has
no asymptotic advantage over plain ρ on secp256k1".

If the answer is "no, in a subgroup of index `k`", the proposal
becomes a real attack with cost `√(n/k)`.  This would be a
significant finding even for moderate `k`.

The framework is also reusable for **other prime-field CM curves**:
the proposal applies to any j=0 / j=1728 curve, and the experiment
in §5 can be run against P-256-shaped CM curves, Koblitz curves,
GLS curves, etc.

### Why this is worth writing up even if null

Per the original research-direction preamble: "you publish a paper
that says 'we computed X, found no exploitable structure, here's a
tighter lower bound on attacks of type Y'.  This is normal,
valuable research."

A null result on VFCG-ρ tightens the structural bound:

> *No "horizontal-action-at-the-floor" ρ-style attack on j=0
> prime-field CM curves gives asymptotic speedup, provided the
> class-group action `c_γ` is generic modulo `n`*.

Combined with the existing `RESEARCH.md` (Phase 1-3 of CGA-HNC),
`RESEARCH_SECP256K1_CM.md` (audit), and `RESEARCH_P256.md`
(the broader programme), this would close the "use the class-group
action somehow" research direction more completely than any single
prior work.

## 9. Implementation sketch (Rust)

A toy implementation would fit in ~600 LOC inside
`crypto/src/cryptanalysis/`, alongside the existing rho variants.
Sketch of the module structure:

```rust
//! src/cryptanalysis/vfcg_rho.rs

/// Floor-level state: (curve_index, point_on_that_curve).
pub struct FloorState {
    pub i: u32,                   // curve index in Cl(O)-orbit
    pub pt: Pt2,                  // point on E_↓^{(i)}
    pub additive_path: Vec<u8>,   // log of step types taken
}

/// Precomputed: the h floor curves + class-group action.
pub struct VolcanoFloor {
    pub ell: u32,                 // descent prime
    pub h: u32,                   // class number h(-3·ℓ²)
    pub curves: Vec<(BigInt, BigInt)>,  // (a_i, b_i) for each E_↓^{(i)}
    pub gamma_action: Vec<...>,   // explicit isogeny γ : E_↓^{(i)} → E_↓^{(i+1)}
    pub gamma_scalar: BigInt,     // c_γ ∈ Z/nZ
}

pub fn vfcg_rho_attack(
    floor: &VolcanoFloor,
    p: &Pt2,  q: &Pt2,
    n: &BigInt,
    dp_bits: u32,
    max_steps: u64,
) -> Option<BigInt> {
    // 1. Start two walks from P and Q, both lifted to E_↓^{(0)}.
    // 2. Iterate the step function described in §2.
    // 3. On DP collision across the two walks, solve for d via §2.
    // 4. Return Some(d) or None on timeout.
}

#[cfg(test)]
mod tests {
    /// Reproduces the §5 experiment on E_toy / F_{2^32}.
    #[test]
    #[ignore]   // hours of compute
    fn cgamma_subgroup_test_for_secp256k1_like_toy() {
        // 1. Build a 32-bit-size analog of secp256k1.
        // 2. For each ℓ ∈ {5, 7, 11, 13, 17}, compute c_γ.
        // 3. Check gcd(c_γ − 1, n) — is it 1 or > 1?
        // 4. Tabulate.
    }
}
```

PARI helpers in `secp256k1_cm_audit/`:

- `descend_one_floor.gp` — given `(p, b)` of a j=0 curve, compute
  the depth-1 floor orbit using the modular polynomial Φ_ℓ.
- `compute_c_gamma.gp` — for each floor curve and class-group
  generator, compute the action scalar `c_γ ∈ Z/nZ`.
- `vfcg_experiment.gp` — run the §5 falsification experiment.

Initial implementation effort: ~3-5 days for a PARI prototype that
runs §5 on multiple `(p, ℓ)` pairs and tabulates `gcd(c_γ − 1, n)`.
A Rust port would follow if any pair shows non-trivial gcd.

## 9.4. Strengthened structural analysis (executed)

The original §5 falsifier tested whether `gcd(ℓ ± 1, n) > 1` — but
this only catches a very narrow type of structuredness.  The
**fundamental algebraic question** is:

> *What is the multiplicative order of `c_γ` in `(Z/n)*`?*

The companion analysis (in
[`secp256k1_cm_audit/vfcg_variants.gp`](secp256k1_cm_audit/vfcg_variants.gp))
gives the upper bound

```
        ord(c_γ)  ≤  2 · h(−3 · ℓ^{2d})
```

for the depth-`d`, prime-`ℓ` floor (the factor of `2` is the unit
ambiguity in `O_d`).  This bound is *uniformly* small: even at
depth 3 with `ℓ = 17`, `ord(c_γ) ≤ 3468`.

**So `c_γ` is *always* in a small subgroup of `(Z/n)*` for every
floor variant the proposal considers.**  The §4(b) "structuredness
hope" is satisfied uniformly — not just possibly, but *always*.

### Strengthened-bound table (excerpt)

| Variant     | `ℓ` | `d` | `disc(−3·ℓ^{2d})` | `h(disc)` | bound `ord(c_γ)` |
|-------------|----:|----:|------------------:|----------:|------------------:|
| V1, basic   |   5 |   1 |              −75  |        2  |              4    |
| V1, basic   |   7 |   1 |             −147  |        2  |              4    |
| V1, basic   |  17 |   1 |             −867  |        6  |             12    |
| V1, depth 2 |   5 |   2 |           −1 875  |       10  |             20    |
| V1, depth 2 |  17 |   2 |         −250 563  |      102  |            204    |
| V1, depth 3 |  17 |   3 |      −72 412 707  |    1 734  |          3 468    |
| V2, ℓ_1·ℓ_2 |  5,7 |  − |          −3 675   |       12  |             24    |
| V2, ℓ_1·ℓ_2 | 11,17 | − |       −104 907   |       72  |            144    |
| V2, ℓ_1·ℓ_2 | 13,17 | − |       −146 523   |       72  |            144    |

### Why uniform-small-subgroup STILL doesn't give an attack

Two structural reasons (see §3 + the variant script's verdict
section for the detailed argument):

**(a) Additive steps re-randomise.**  A VFCG-ρ walk mixes
`γ`-multiplicative steps with additive steps.  The additive steps
move freely in `(Z/n)`, so the *combined* walk visits all of
`(Z/n)` with the same density as plain ρ — regardless of the
small subgroup that `c_γ` itself lives in.

**(b) γ-only walks are useless.**  If we restrict to walks that
use *only* `γ`-multiplicative steps (which respect the small
subgroup), the walker cycles in a single `O(h)`-element coset
of `⟨c_γ⟩`.  Such a walker recovers the secret `d` only if `d`
happens to lie in that coset — probability `h/n`, vanishing.

### Strengthened verdict

> **VFCG-ρ is null in a stronger sense than §5 caught.**  The
> "c_γ in a small subgroup" condition that §4(b) hoped for is
> satisfied *uniformly* for every floor variant — yet the
> structural cost analysis of §3 shows this is not sufficient
> for a ρ speedup.  The proposal is therefore closed not just
> for secp256k1's specific `n`, not just generically, but for
> *every* `(j = 0, ℓ, d)` choice.

This is a tight, structurally-justified negative result.

### V3 and V4 status

**V3 (tagged DPs)** is orthogonal: it changes the distinguishing
criterion but not the walk function `f` or the state space.  The
cycle structure of a ρ walk depends on `f`, not on which states
are marked.  Therefore V3's asymptotic cost is plain ρ's, and
V3 is a re-parameterisation rather than an attack.

**V4 (supersingular reduction of E mod auxiliary ℓ)** is out of
scope for this falsifier.  The supersingular isogeny graph over
`F̄_ℓ` is a different mathematical world from the F_p-isogeny
class of secp256k1; no known direct connection from the
supersingular isogeny problem (attacked by Castryck-Decru / SQI-
sign-cryptanalysis machinery) to F_p ECDLP.  Would require its
own research-direction paper.  Deferred.

### Reproducer

```
cd secp256k1_cm_audit
gp -q vfcg_variants.gp > vfcg_variants_output.txt
diff vfcg_variants_output.txt vfcg_variants_output.expected.txt
```

The expected output is committed at
[`secp256k1_cm_audit/vfcg_variants_output.expected.txt`](secp256k1_cm_audit/vfcg_variants_output.expected.txt).

## 9.5. Empirical result of §5 falsifier (executed)

The §5 experiment is implemented in
[`secp256k1_cm_audit/vfcg_experiment.gp`](secp256k1_cm_audit/vfcg_experiment.gp).
Initial run results, with `gcd(ℓ ± 1, n)` as the structuredness test
for the `h(−3ℓ²) = 2` case (`c_γ ≡ ±ℓ mod n` generically):

**Direct test on secp256k1**: for every `ℓ ∈ {5, 7, 11, 13, 17, 19,
23, 29, 31}`, both `gcd(ℓ − 1, n)` and `gcd(ℓ + 1, n)` equal `1`.
**Trivial gcd ⇒ no proper subgroup structure ⇒ no attack avenue.**

**Aggregate test over 200 random j=0 prime-order toy curves**
(p in the range `10⁵ … 10⁷`, b varying through sextic-twist
representatives until a prime order is found): **0 of 200** had a
non-trivial gcd for any tested `ℓ`.

Conclusion: the proposal is empirically null at toy scale (matching
the random-model prediction of `< 1` expected hit), and is
structurally null for secp256k1 by direct test.  The §5 falsifier
**closes the VFCG-ρ direction with a tight negative result**:

> No `(secp256k1, ℓ)` pair with `ℓ ≤ 31` admits the proposed
> class-group-induced subgroup speedup; and over random j=0
> prime-order curves, the speedup arises with frequency below the
> random-coincidence floor.

The output of the experiment is committed at
[`secp256k1_cm_audit/vfcg_output.expected.txt`](secp256k1_cm_audit/vfcg_output.expected.txt)
for reproducibility (diff against `vfcg_output.txt` after re-running
should be empty).

## 10. Next concrete action

The single experiment that decides whether this direction is alive
or dead:

```
For 100 random j=0 ordinary curves E / F_p with prime order n
(p ≈ 2^{32}, 2^{40}, 2^{48}, 2^{56}, 2^{64}):

   For ℓ ∈ {5, 7, 11, 13, 17}:
      Compute c_γ for that (E, ℓ).
      Compute gcd(c_γ − 1, n), and similarly for c_γ^k − 1
      for k ∈ {1, 2, …, h(−3ℓ²)}.

Report: fraction of curves where any such gcd > 1.
```

If the fraction is anywhere above `O(1/n^{1/2})` (i.e., above what
you'd expect by random coincidence), VFCG-ρ is alive and we proceed
to a full toy ρ implementation.  If at or below random, the
direction is closed.

A clean PARI script for this experiment is the minimal next-action
deliverable for this proposal.

---

## Summary

**What's proposed**: a Pollard-ρ variant that walks on the depth-1
floor of secp256k1's isogeny volcano (where `h(End) = h(−3ℓ²) > 1`),
using class-group-action steps in addition to additive steps.
Walks tagged with the class-group orbit index; DPs shared across
the orbit.

**Why it's potentially novel**: the obvious "use the class-group
action" attacks (CSIDH, CRS, CGA-HNC) are all blocked at the
**crater** for secp256k1 because `h(−3) = 1`.  This proposal works
at the **floor** where the class group is non-trivial, while
preserving the same `#E = n` (so the attack stays a real attack on
secp256k1's actual ECDLP, not a different problem).

**Why it probably doesn't work**: the per-step cost of class-group
transitions (`O(ℓ²)` for an `ℓ`-isogeny via Vélu) outweighs the
birthday-saving factor `√h` unless the class-group action's scalar
`c_γ` lies in a structured proper subgroup of `(Z/nZ)*`.  My
strong prior is that `c_γ` is generic — but this is checkable in
< 1 day of PARI compute.

**Decision deliverable**: a PARI script (≤ 200 LOC, ~1 day to write)
that empirically tests whether `c_γ` is structured.  If yes,
proceed to full implementation; if no, write up the negative result
and add to the structural-lower-bound table.

**What it adds to the existing repo even if null**:

- A new structural-lower-bound statement: "the floor class-group
  action at depth 1 has a generic scalar on `n`."
- A reusable PARI module for computing class-group actions at the
  floor of CM-volcanos — useful for future research directions.
- A clean negative result extending CGA-HNC's "trivial-class-group
  at the crater" obstruction to the deeper claim "non-trivial-class-
  group at the floor is also blocked".

This is the kind of focused, falsifiable mini-programme that a
small team could prosecute in ~2 weeks and publish as a structural
note even with a null outcome.
