# Open IC leads after 2026-10-11: what each costs, and what was settled

This is a search round over the five open leads of
[`RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md`](../../RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)
§9 and the round-1 tree-split follow-ups. Three were run or calculated in this round, and
two were scoped:

| lead | status after this round | evidence |
|:--|:--|:--|
| 1. `m = 4` solving degree at `ℓ = 4…6` | **scoped, not run.** Out of reach on any host the repo has used. The decisive cells and their prices are below | §1 (derived) |
| 2. Equivariant F4 extended to full F4 | **scoped, not run.** Engineering of 3–5 days, and it needs a Frobenius-invariant `m = 4` formulation that does not exist | §2 (from the repo's own scope notes) |
| 3A. Preprocessing rho at equal advice | **calculated.** The IC online class loses online by about 48× in operations at n = 73 and about 3,400× at n = 83. The note's n = 73 figure used the rank-stage probe count | [`../ic_preprocessing_rho_equal_advice_20261011/REPORT.md`](../ic_preprocessing_rho_equal_advice_20261011/REPORT.md) (derived) |
| 3B. Free-relation census | **measured.** Distinct relations, as the rank mod `r`: 49 at n = 61, 1 at n = 73 and 6 at n = 83, against 599 needed. The rank stage is not free, and the note's count is about `n×` too high | [`../ic_free_relation_census_20261011/RESULTS.md`](../ic_free_relation_census_20261011/RESULTS.md) (measured) |
| 4. D₄-symmetrised `S₅` lattice at `m = 4` | **settled by a counting argument; not run.** Symmetrisation cannot move the exponent of `δ`, and 1/80 was never a measurement | §3 (derived) |
| 5. Tree split, round 2 | **pre-registered and run** | [`../ic_tree_split_structured_v2_20261011/RESULTS.md`](../ic_tree_split_structured_v2_20261011/RESULTS.md) |

## 1. Lead 1: which solving-degree cell decides, and what it costs

**What is measured.** The `m = 4` chained-`S₃` ladder
(`research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md`, Result 10,
`research/dreg_m4_chain_20261006/`) has:

| `ℓ` | cells | measured `D` |
|--:|:--|:--|
| 1 | (4–6, 1) | 5 |
| 2 | (5–8, 2) | 6 |
| 3 | (9, 3), (10, 3) | ≥ 7. Degree 7 was out of reach: about 2.1M × 2.8M |

Its fit is `D = ℓ + 4`. The alternative the leads note keeps open is a degree that levels
off at 7.

**The point the note leaves implicit.** At `ℓ ≤ 3` the two hypotheses predict the same
values: 5, 6 and 7. The first `m = 4` cell that separates them is at `ℓ = 4`, where
`ℓ + 4` gives 8 and levelling gives 7. So finishing `(9, 3)` at degree 7 would confirm 7,
but it cannot decide between them.

**Prices** (Macaulay columns `Σ_{k≤D} C(N, k)`; the agent's arithmetic, spot-checked here:
`C(39,8) = 61,523,748`). `N` is the unknown count: `(m−2)n + mℓ` for the chained system, or
`4(ℓ−1) + 1 + 2n` for the symmetrised-chained system
(`koblitz_symmetrised::build_chained_symmetrised_system`).

| cell | what it decides | `N` | columns at the deciding degree | host |
|:--|:--|--:|--:|:--|
| `m = 3` (15, 5) at `D = 7` | levelled at 7 against the semi-regular 8, at `m = 3` | 30 | band after 7 about 110–120 GB (DREG note) | about 128 GB, the cheapest decisive cell in the programme |
| `m = 4` (9, 3) at `D = 7` | confirms 7 at `ℓ = 3`; **does not separate** the hypotheses | 30 | 2.8M | large host |
| `m = 4` (13, 4), symmetrised, at `D = 7` | does the system close at 7? If it does, `ℓ + 4` is refuted | 39 | about 1.9·10⁷ | sparse F4/F5 on TB-class memory; dense is out |
| `m = 4` (13, 4), symmetrised, at `D = 8` | `ℓ + 4` | 39 | about 8.1·10⁷ | out of reach |
| `m = 4` (17–20, 5) and (21–24, 6) at `ℓ + 4` | the registered ask | 51–69 | 4·10⁹ to 6·10¹¹ | out of reach |

**Recommendation.** Run `m = 3` (15, 5) at degree 7 on a 128–256 GB host, pre-registered,
with the existing `dreg_ladder` harness. It is the cheapest cell that tests "levels off at 7"
at all.

- If `m = 3` does not level, the `m = 4` levelling hypothesis loses its only analogue.
- If it does level, the next step is the `m = 4` (13, 4) degree-7 closure test. That needs
  a sparse engine and Frobenius blocks (lead 2), not more memory.

Neither fits this container (4 vCPU, 15 GB).

## 2. Lead 2: equivariant F4 to full F4

**What exists.**
- `examples/g2_equivariant_orbit_probe.rs` (n ≤ 11) and `examples/g2_wide_probe.rs`
  (n ≤ 23); outputs are in `research/g2_equivariant_f4_20261008/`.
- They build one degree-`D` Macaulay matrix for the `m = 2` barrel `S₃` system, and check
  that the trivial block plus the character blocks (a DFT over `F_{2^d}`) reproduce its rank.
  That rank identity is the measured 31–128×.
- There is no S-pair loop, no inverse transform of reduced rows, and no `src/` module.

**What "full F4" needs** (`RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md`, steps 4–5, estimated at
3–5 days):
- the inverse DFT back to Boolean polynomials;
- a `reduce_system` engine in `koblitz_groebner`, checked against the current engine;
- a carry-less `F_{2^d}` kernel for `d` up to 130;
- extension beyond `m = 2`.

**The blocker for lead 1.**
- The character blocks need a Frobenius-invariant system: the barrel or normal-basis
  formulation with selectors.
- The solving-degree ladders deliberately use random subspaces that are not
  Frobenius-stable (`koblitz_symmetrised.rs`, `random_subspace_containing_one_in`).
- So the blocks cannot shrink lead 1's matrices until an `m = 4` barrel formulation exists.
  That is a third piece of work, before the engine.
- And the speedup is a constant (≤ `n`), not an exponent.

**Recommendation.** Build it only as the enabler of the `m = 4` (13, 4) degree-7 closure
test, and only if `m = 3` (15, 5) shows levelling.

## 3. Lead 4: D₄-symmetrised `S₅` lattice, settled without a run

**What 1/80 is.** It is not a measurement. It is the asymptotic Jochemsz–May bound
`δ_max = 1/(m(m+1)·2^{m−2})` for a polynomial filling `[0, 2^{m−1}]^m`
(`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md` §2.2). For `S₅` that gives 1/80, against
the 1/6 the factor base needs. The only lattices ever built are `S₃` (`m = 2`, Probes B
and C, 40- and 56-bit curves):

| quantity | value |
|:--|:--|
| rigorous reach | 0.10 |
| heuristic reach | 0.16–0.18, against the predicted 1/6 |
| symmetric `(e₁, e₂)` lattice threshold | the same, only 3–10× cheaper per reduction |

**What D₄ is.** It is `(Z/2)³ ⋊ S₄`: translating an even number of summands by a rational
2-torsion point, plus permutations. The "630×" figure in the leads note is a Bézout-count
estimate for Gröbner solving over `F_{p⁴}`
(`RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`), not a lattice reach.

**The counting argument.**
1. **The elementary-symmetric rewrite.** It maps the small root `(x₁…x_m)` with
   `|x_i| ≤ B` to `e_k ≤ C(m,k)·B^k`. Each `e`-monomial has the same size bound as the
   `x`-monomials it replaces, `B` to the weighted degree. So the Howgrave-Graham
   determinant condition is unchanged, and only the lattice dimension shrinks, by about
   `m!`. Probe C measured exactly this at `m = 2`: the same threshold, a cheaper reduction.
2. **The `(Z/2)³` part.** In Edwards `y` it acts by `y ↦ −y`. Its invariants `y²` halve the
   degree but square the bound, so `X^d` is unchanged, and so is `δ`. In Weierstrass `x`,
   translation by `T` is a Möbius map that sends small `x` to large values. So the
   invariants are not small at the roots, and there is no small-root problem to solve.
3. **It needs special curves.** Rational 2- or 4-torsion excludes prime-order curves, the
   deployed case, outright.
4. **The gap is too big.** It is a factor of 13 in the exponent of `δ` (1/80 against 1/6).
   A constant-factor change in dimension cannot close it.

**What would reopen it.** A small-root method whose reach does not follow the
Newton-polytope determinant count. That is a new algebraic idea, not this probe.

**Optional sanity check** (about an hour, not done here): add an Edwards-symmetric `S₃`
lattice to `examples/prime_ecdlp_breakthrough_probe.rs` `probe_c`, and confirm its threshold
matches the existing symmetric lattice, which hits zero at about 0.18.

## What this round does not establish

- No exponent, speedup or deployed-curve claim.
- Leads 1 and 2 are priced and scoped, not run.
- Lead 3A is a calculation. Its ledger consequence needs a measured BL run on a curve.
