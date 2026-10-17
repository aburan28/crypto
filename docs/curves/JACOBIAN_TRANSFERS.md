# Jacobian transfer discovery: every route, every registered curve, one unit

[`registry.json`](registry.json) names the curves; [`COVERS.md`](COVERS.md)
proves that same-field covers exist; [`GHS_SCREEN.md`](../GHS_SCREEN.md)
screens one binary curve at a time.  None of them says **whether moving a
curve's discrete logarithm into a Jacobian could ever pay**.
[`jacobian-transfers.json`](jacobian-transfers.json) does: for every curve in
the registry it enumerates the known transfer routes, computes what is exact
about each (magic numbers, genera, map degrees, torsion signals), prices the
destination with one frozen formula-level cost model, and reports the ratio to
Pollard rho on the original subgroup.  `jacobian_transfers` rebuilds and
checks the file; CI fails when the registry gains a curve the file has not
screened.

```bash
cargo build --release --bin jacobian_transfers
./target/release/jacobian_transfers build                 # registry.json → jacobian-transfers.json
./target/release/jacobian_transfers check                 # CI: the committed file is current
./target/release/jacobian_transfers show "Oakley Group 3" # every route of one curve
./target/release/jacobian_transfers rank --top 20         # curves by best ratio
./target/release/jacobian_transfers rank --class formula_advantage
./target/release/jacobian_transfers table                 # the Markdown table below
./target/release/jacobian_transfers chart --svg chart.svg --tikz chart.tex
```

Implementation: [`src/cryptanalysis/jacobian_transfer/`](../../src/cryptanalysis/jacobian_transfer/mod.rs)
(`cost.rs` the model, `binary.rs` the GHS and isogeny-extension screens,
`report.rs` the table and charts), [`src/bin/jacobian_transfers.rs`](../../src/bin/jacobian_transfers.rs),
workflow [`jacobian-transfers.yml`](../../.github/workflows/jacobian-transfers.yml).

## What a row is, and is not

A row is a **formula**.  Its `ratio_log2` is `log2(rho on E) − log2(route)`
under the cost model in the file's `cost_model` block; its
`exponent_ratio_log2` is the same difference with every constant stripped
(`0.5·log2 r` against the solver's leading term), so an exponent gain can be
told from a constant.  The classes are:

| class | meaning |
|:--|:--|
| `formula_advantage` | the literature's complexity formulas, with the stated constants, put the route at least `advantage_margin_bits` (8) below rho and the screen's necessary size/degree conditions do not rule out subgroup transport. This is only a conditional theoretical lead: it is not a verified catalog transfer, speedup, measurement, or security claim. |
| `formula_marginal` | the formulas place the route below rho by less than the margin: the exponents favour it and the constants decide, so the row is as strong as the model's constants.  Read `exponent_ratio_log2` beside `ratio_log2`. |
| `boundary` | the formulas themselves place the route at or above rho.  No constant-factor engineering on the route can cross rho, so the row is a boundary in the sense of `AGENTS.md` §1. |
| `not_carried` | the destination cannot contain the subgroup (`too_small`), is `E` itself over a subfield (`degenerate`), or does not exist over the subfield (`descent_obstructed`, MMT Lemma 6 below): the route solves a different problem. |

Each route also carries an `evidence` level: `exact_invariant` (the geometry
was computed on this curve; only the solver cost is a formula),
`construction` (a known construction family supplies the destination; not
exercised here), `proposed` (existence is open; the row prices the most
favourable genus the literature allows, so it is a lower bound on the route),
`heuristic` (counting and class-membership heuristics enter the price).

Every obligation in [`../endomorphism-rules.md`](../endomorphism-rules.md) stays
open on a `formula_advantage` row: the cover is not constructed, the map's
action on the subgroup is not verified, no relation is collected, nothing is
timed.  Turning a row into a result means a separately registered,
end-to-end experiment under the `AGENTS.md` boundary/table/ratio rules.

## Routes

| field | route (`kind`) | exact part | estimated part |
|:--|:--|:--|:--|
| `GF(2^N)` | `ghs_weil_descent`, one per factorisation `N = n·l` including `l = 1` | magic number `m = dim Span{(1, σ^i(√b))}`, type, genus `2^(m−1)` or `2^(m−1) − 1`, cover degree `2^m` | solver on `Jac(C)/GF(2^l)` |
| `GF(2^N)` | `ghs_isogeny_extension` (Galbraith–Hess–Smart) | the number of `b′` whose `√b′` has Frobenius annihilator of degree `d` for the tower: divisors of `x^n − 1` of degree `d` times `2^(l·d)` | that the isogeny class is a random set of `√|t² − 4q| / π` curves; the walk cost until one is hit |
| `GF(2^N)` | `same_field_cover` | genus 3, degree 3 ([`COVERS.md`](COVERS.md)) | solver on `Jac(H)/GF(2^N)` |
| `GF(p)` | `same_field_cover` | genus 2, degree 2 ([`COVERS.md`](COVERS.md)) | solver on `Jac(H)/GF(p)` |
| `GF(p)` | `base_change_descent`, `k = 2, 3` | nothing: existence open | solver at the most favourable genus `g = k` |
| `GF(p^k)` | `extension_index_calculus` (Gaudry–Diem, factor base from `GF(p)`) | — | `k! · p^(2 − 2/k)` decomposition attempts, each at the measured `C_k` (`k = 3, 4`) or at least `D³`, `D = 2^(k(k−1))` the Bézout bound of the Semaev system |
| `GF(p^k)` | `cover_descent` to a cover over `GF(p)` | whether `#E` is even (a rational 2-torsion point, the usual precondition at `k = 3`) | solver at the most favourable genus: `g = k`, or the GHS genus floor where that is higher (`genus_basis`) |

**The magic number and `a`.**  `E_(a,b) ≅ E_(a′,b)` whenever `a + a′` is of
the form `α² + α`, so `a` enters only through its absolute trace, and it
never changes the magic number: `m = dim Span{(1, σ^i(√b))}` depends on `b`
alone (Maurer–Menezes–Teske, LMS J. Comput. Math. 5 (2002), eq. (5)).  What
`a` decides is whether the descent exists.  MMT Lemma 6, condition (8): the
cover descends to `F_(2^l)` when `n` is odd, or `Tr(a) = 0`, or
`Tr_(F_(2^N)/F_(2^l))(b) ≠ 0`.  When all three fail, the sum of the `n`
Artin–Schreier conjugates of `z² + z = x + a + b/x²` is the constant
`Tr_(K/k)(a)`, whose root lies in `F_(2^(2l))` but not in `F_(2^l)`.  Every
lift of the `2^l`-Frobenius then has order `2n`, and there is no descent.
Each GHS row records `trace_a`, `relative_trace_b_zero` and
`mmt_condition_8`.  A row where (8) fails has transport
`descent_obstructed`.  Because `Tr(a)` is an isogeny invariant (it fixes
`#E mod 4`), the isogeny-extension rows for such a tower count only targets
`b′` that satisfy (8): their annihilator contains the full power of `x − 1`
in `x^n − 1`, and their magic number is then exactly the annihilator degree.

Before 2026-10-13 the screen read `Tr(a) = 1` with `n` even as "`m_b` or
`m_b + 1`", priced both as `heuristic`, and had an `unresolved` class for
pairs that disagreed.  Both are gone (schema `jacobian-transfers/2`).
c2pnb176w1, c2pnb272w1 and c2pnb304w1 have `Tr(a) = 1` but satisfy (8) on
every tower, so their `m = 8` rows over `GF(4)` are exact.  The only
registry curve where (8) fails is an 18-bit test curve with `b ∈ F_2`.
`ec_trapdoor::mmt_condition_8` is the same check, and `audit_curve` records
it in `DescentRow::descends`.  Each row records the magic number
`audit_curve` reports as `audit_curve_magic_number`; the row is not priced
with it.

**2026-10-15 correction: the genus a `cover_descent` row is priced at.**  Until this date every `cover_descent` row was priced at genus `k`, the least genus a cover of the Weil restriction can have in general.  The odd-characteristic study ([`research/odd_char_descent_genus3_20261014/`](../../research/odd_char_descent_genus3_20261014/REPORT.md)) shows that no GHS cover of `E/GF(p^k)` with `F_p(j) = GF(p^k)` has genus below 9, 49 and 41 at `k` = 6, 10 and 14 (`ghs_genus_floor`); Joux–Vitse (Eurocrypt 2012, App. C) shapes 12/13 give 9 as the minimum at `k` = 6.  The rows are now priced at the floor, the label says so ("descent to a genus-9 cover …"), `genus_basis` records why, and each keeps what it said at genus `k` as `class_at_genus_k` and `ratio_log2_at_genus_k`.  Of the 846 `cover_descent` rows that were `formula_advantage`, 599 still are (576 at `k` = 6, 23 at `k` = 14), 41 are now `formula_marginal` and 206 are `boundary` (all 126 at `k` = 10, 80 at `k` = 14); of the 46 that were `formula_marginal`, 22 still are and 24 are `boundary`.  Curve verdicts move from 848 / 50 / 1201 (`formula_advantage` / `formula_marginal` / `boundary`) to 601 / 67 / 1431.  The floor is a lower bound for GHS covers only: a curve whose `j` lies in a proper subfield may have a smaller cover, a non-GHS cover is not covered by it, and whether a particular curve attains the floor is verified at no `k`.  *Correction (2026-10-16):* this paragraph first said it was verified at `k` = 6; Joux–Vitse's shapes 12/13 give the minimum genus, not the existence of a genus-9 cover for a given curve, so every `k` = 6 row now carries `conditional_on` = "a genus-9 cover of this curve being exhibited".  Rows with `genus_basis` = `ghs_genus_floor` are still `proposed`: no cover is exhibited.

**2026-10-16 correction: the measured Jacobian operation, the extension index-calculus price, and `k` = 6 rows made conditional.**  A claim audit found three stale or unvalidated inputs.  (1) The Jacobian operation was still `10·g²` multiplications although [`research/jacobian_op_count_20261014/`](../../research/jacobian_op_count_20261014/README.md) (#1846) counted textbook Cantor 2.0 bits (genus 2) to 3.3 bits (genus 32) higher; the model now uses those counts (`jacobian_group_op_mults_measured`).  (2) `extension_index_calculus` priced each decomposition at `D = 2^(k(k−1))` elliptic operations with no trial factor; it now charges `k!` attempts per relation and the measured `C_3 = 0.876·10⁶`, `C_4 = 1.227·10¹²` `F_p` multiplications, or at least `D³` where no `C_k` is measured.  (3) The `k` = 6 `cover_descent` rows carry `conditional_on` = "a genus-9 cover of this curve being exhibited": Joux–Vitse's shapes 12/13 give 9 as the minimum genus, not existence.  Effect on the 815 `formula_advantage` routes: **598 survive** — 574 `k` = 6 `cover_descent` (every one conditional on an unexhibited genus-9 cover), 22 `k` = 14 `cover_descent` at the genus-41 GHS floor (also `proposed`, no cover exhibited) and 2 binary `ghs_weil_descent`.  All 214 `extension_index_calculus` advantages and all 173 marginals are now `boundary` (best ratio −15.0 bits); 3 `cover_descent` advantages fell to `formula_marginal`; 12 `cover_descent` and 1 `ghs_weil_descent` marginal are `boundary`.  Curve verdicts move from 601 / 67 / 1444 to 598 / 57 / 1457.  Only the 2 binary GHS rows (`exact_invariant`) do not rest on a cover nobody has built; the other 596 are `proposed`.

Routes the scan does **not** enumerate, and why: Hess's generalised
Artin–Schreier covers beyond the classical GHS family (needs a per-curve
function-field search); pairing transfers to `GF(q^k)^*` (a different
destination, covered by the `embedding` trait in [`TRAITS.md`](TRAITS.md));
G2 groups of pairing curves (not in the registry).  The isogeny-extension
rows never exhibit an isogenous curve: they say how many are expected.

## The cost model

All prices are `log2` of estimated bit operations, so source and destination
share a unit.  The parameters are written into the file and are estimates.

| quantity | formula |
|:--|:--|
| field multiplication, `w` bits | `w^1.585` (Karatsuba) |
| elliptic group operation | 10 multiplications |
| genus-`g` Jacobian operation | the counted textbook-Cantor multiplications of [`research/jacobian_op_count_20261014/`](../../research/jacobian_op_count_20261014/README.md) (#1846; addition over `F_1009`, inversion = 8 multiplications): 164.7 at `g = 2` … 102,282.5 at `g = 32`, interpolated in `log2 g`, extrapolated as `g^2.347` past 32 (a lower bound there: the local exponent is about 2.75); `g = 1` is the elliptic constant |
| reference: rho on `E` | `√(πr/4)` elliptic operations over `GF(q)` |
| `rho_on_jacobian` | `√(πr/4)` Jacobian operations over the destination field |
| `hyperelliptic_index_calculus` | `g! · q + g² · q^(2 − 2/g)` Jacobian operations: Gaudry's relation collection (a divisor splits over the degree-one factor base with probability about `1/g!`) plus double-large-prime linear algebra (Gaudry–Thomé–Thériault–Diem); `g ≥ 2`.  The `g!` term confines the method to genus small beside `q`. |
| `enge_gaudry_smoothness` | Enge–Gaudry priced by an explicit count, for genus up to 1024: factor base of the places of degree at most `b` (half the irreducibles of each degree), relations at the exact fraction of `b`-smooth monic polynomials of degree `g`, sparse linear algebra `g·F²` operations modulo `r`; minimised over `b` |
| `enge_gaudry` | `L_N(1/2, √2)` Jacobian operations, `N = q^g`; only above genus 1024 and when `g ≥ log2 q` |
| `extension_index_calculus` | `k! · p^(2 − 2/k)` decompositions, each `C_k` multiplications over `GF(p)`: measured `C_3 = 0.876·10⁶`, `C_4 = 1.227·10¹²` (`RESEARCH_RESIDUAL_WALKS.md` §11.8, §11.17), else `D³` with `D = 2^(k(k−1))` (a floor; both measured values lie above it) |
| transfer | `g` elliptic operations over the source field (one conorm–norm evaluation) |
| search (isogeny extension) | `2^(N+1) / #candidates` isogeny steps, each one multiplication over `GF(2^N)` times `N` |

The destination price is the minimum over the solvers that apply; the route's
total is the log-sum of search, transfer and destination.  The subgroup is
`size_admits` when `g · log2(q_dest) + 1 ≥ log2 r`, a necessary condition
only. For an actual finite same-field cover with verified
`norm ∘ pullback = [d]` and exact subgroup order `r`, coprimality of `r` and
the map degree `d` would make pullback injective on that subgroup. The scan's
“carried” label records this conditional theorem input; it is not an executed
map or subgroup certificate. A route on `E` itself (`extension_index_calculus`,
`base_change_descent`) has `same_group` transport.

Changing a parameter changes every row; `check` rebuilds with the committed
file's parameters, so a parameter change is a reviewed diff of the whole
file.

## Reading the result

`jacobian_transfers rank` orders curves by the best route classified as
conditionally carried by this screen. Under
the default model the registry splits as the file's `summary` block says;
the curves with a `formula_advantage` verdict are listed there by slug with
the route that gives it, and the chart in
[`../../research/jacobian_transfer_scan_20261010/`](../../research/jacobian_transfer_scan_20261010/README.md)
plots every curve's best ratio against its subgroup size.

Two readings are safe and the rest are not:

- A **`boundary` verdict on every route is a closed question for those
  routes**: the formulas, not the engineering, put every enumerated transfer
  at or above rho.  That is the state of every deployed prime-field curve
  and every prime-degree binary curve in the registry, agreeing with
  `research/cover_route_prime_binary_20261008/`.
- A **`formula_advantage` or `formula_marginal` verdict is a lead**, with the
  route's `evidence` level saying how much of it is established.  See
  [`research/jacobian_transfer_scan_20261010/`](../../research/jacobian_transfer_scan_20261010/README.md)
  for the curves in each class and what each still needs.

Not safe: reading `ratio_log2` as a speedup (it is a formula ratio, §2 of
`AGENTS.md` admits only end-to-end measurements as speed); comparing rows
across different cost-model parameters; treating a heuristic
`ghs_isogeny_extension` row as the existence of an isogenous weak curve;
treating a `proposed` row's genus as attained.

## Summary table

Generated by `jacobian_transfers table` from the committed file: every curve
with a standard name and every `formula_advantage` curve, by best ratio.

<!-- jacobian_transfers table -->
| curve | field | r bits | rho log2 | best conditionally carried route | genus | total log2 | exponent ratio | ratio log2 | class |
|:--|:--|--:|--:|:--|--:|--:|--:|--:|:--|
| c2pnb272w1 (`icv1-f2m272-t114182707820475866940422329614334265958971-985955bd`) | GF(2^272) | 256 | 144.0 | GHS descent to Jac(C)/GF(2^2), N = 136 * 2 | 128 | 74.9 | +69.9 | +69.1 | formula_advantage |
| `icv1-fp42k6-t45058067596985537030509149035791573678-385b2709` | GF(p^6), p 42 bits | 250 | 140.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.8 | +39.8 | +43.0 | formula_advantage |
| `icv1-fp42k6-tm127882077426822819848891979507313072634-ce17ec21` | GF(p^6), p 42 bits | 250 | 140.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.8 | +39.8 | +43.0 | formula_advantage |
| `icv1-fp40k6-tm292901180069185163850972135730463082-86d9ad4b` | GF(p^6), p 40 bits | 238 | 134.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +37.8 | +41.0 | formula_advantage |
| `icv1-fp42k6-t8930769277479864512670044108354100538-1ae4bc32` | GF(p^6), p 42 bits | 244 | 137.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +37.0 | +40.2 | formula_advantage |
| `icv1-fp42k6-tm44110276869417112049281439078164335142-3a830fff` | GF(p^6), p 42 bits | 244 | 137.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +36.6 | +39.9 | formula_advantage |
| `icv1-fp42k6-tm44110276869417112049281439078164335142-5b7c78bb` | GF(p^6), p 42 bits | 244 | 137.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +36.6 | +39.9 | formula_advantage |
| `icv1-fp40k6-tm13709315706324085134265846307139378-8f5afad0` | GF(p^6), p 40 bits | 236 | 133.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +36.7 | +39.9 | formula_advantage |
| `icv1-fp40k6-t598456382823960921215959826004744218-c2be60b4` | GF(p^6), p 40 bits | 235 | 133.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +36.5 | +39.7 | formula_advantage |
| `icv1-fp40k6-t1990693073699655257162695508951311630-2e0a5586` | GF(p^6), p 40 bits | 235 | 133.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +36.4 | +39.6 | formula_advantage |
| `icv1-fp40k6-t1112861023040603223163073207108369490-de187046` | GF(p^6), p 40 bits | 235 | 133.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +36.3 | +39.5 | formula_advantage |
| `icv1-fp40k6-tm1079691786106281078201479722463943470-5ed88ea1` | GF(p^6), p 40 bits | 235 | 133.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +36.3 | +39.5 | formula_advantage |
| `icv1-fp40k6-tm1079691786106281078201479722463943470-df6a6e29` | GF(p^6), p 40 bits | 235 | 133.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +36.3 | +39.5 | formula_advantage |
| `icv1-fp40k6-tm290230737312338524932089875922384242-794dd1c9` | GF(p^6), p 40 bits | 235 | 133.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +36.2 | +39.4 | formula_advantage |
| `icv1-fp40k6-t1022985513504314449990276264241491654-0e5eea4d` | GF(p^6), p 40 bits | 235 | 133.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.7 | +36.1 | +39.3 | formula_advantage |
| `icv1-fp42k6-tm35931195335088482800968663222430596806-21d8f293` | GF(p^6), p 42 bits | 242 | 137.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +36.1 | +39.3 | formula_advantage |
| `icv1-fp42k6-t42108441828852030315173363066246162586-8161411a` | GF(p^6), p 42 bits | 242 | 136.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +35.9 | +39.1 | formula_advantage |
| `icv1-fp42k6-tm20881010826321473227414487444218254566-25483bc6` | GF(p^6), p 42 bits | 242 | 136.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +35.6 | +38.9 | formula_advantage |
| `icv1-fp42k6-tm20881010826321473227414487444218254566-4772a700` | GF(p^6), p 42 bits | 242 | 136.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +35.6 | +38.9 | formula_advantage |
| `icv1-fp40k6-tm189981468464338651169319038939670950-c99d5cd2` | GF(p^6), p 40 bits | 234 | 132.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +35.6 | +38.8 | formula_advantage |
| `icv1-fp42k6-tm139397368184640235945825692291306341094-56623779` | GF(p^6), p 42 bits | 241 | 136.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +35.4 | +38.7 | formula_advantage |
| `icv1-fp42k6-tm139397368184640235945825692291306341094-f7af9aa6` | GF(p^6), p 42 bits | 241 | 136.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +35.4 | +38.7 | formula_advantage |
| `icv1-fp38k6-tm28555588863730399043033619153456926-332a0998` | GF(p^6), p 38 bits | 225 | 128.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.6 | +35.3 | +38.5 | formula_advantage |
| `icv1-fp40k6-tm1138014948113436022385175216268138942-11cac4ac` | GF(p^6), p 40 bits | 233 | 132.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +35.2 | +38.5 | formula_advantage |
| `icv1-fp40k6-t665405924970026162701491545965324990-f8473e80` | GF(p^6), p 40 bits | 233 | 132.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +35.2 | +38.4 | formula_advantage |
| `icv1-fp40k6-t1863913783112918683701033322387770322-5b3a9c38` | GF(p^6), p 40 bits | 233 | 132.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +35.2 | +38.4 | formula_advantage |
| `icv1-fp42k6-t1449329879519965966669504876259048666-13b27537` | GF(p^6), p 42 bits | 240 | 136.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +35.0 | +38.3 | formula_advantage |
| `icv1-fp42k6-t1449329879519965966669504876259048666-daed92e0` | GF(p^6), p 42 bits | 240 | 136.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +35.0 | +38.3 | formula_advantage |
| `icv1-fp40k6-t494811548722146648595995475931754494-796b5505` | GF(p^6), p 40 bits | 232 | 131.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +35.0 | +38.2 | formula_advantage |
| `icv1-fp38k6-t11958786489917378873825658058094890-9fb5c761` | GF(p^6), p 38 bits | 224 | 127.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +34.8 | +38.0 | formula_advantage |
| `icv1-fp38k6-t39647719120351838730753721191063210-e0475c1a` | GF(p^6), p 38 bits | 224 | 127.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +34.8 | +38.0 | formula_advantage |
| `icv1-fp38k6-t39647719120351838730753721191063210-e5042833` | GF(p^6), p 38 bits | 224 | 127.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +34.8 | +38.0 | formula_advantage |
| `icv1-fp40k6-tm1048669530144502480127388528082528010-cd064fac` | GF(p^6), p 40 bits | 232 | 131.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +34.7 | +38.0 | formula_advantage |
| `icv1-fp38k6-tm155314119424755961078517515381458-93c1dbde` | GF(p^6), p 38 bits | 224 | 127.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +34.7 | +37.9 | formula_advantage |
| `icv1-fp42k6-tm57664958956439933191375515704719752298-714af95b` | GF(p^6), p 42 bits | 239 | 135.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +34.2 | +37.5 | formula_advantage |
| `icv1-fp42k6-t135044599113790442697853745259521825354-23ab444a` | GF(p^6), p 42 bits | 239 | 135.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +34.2 | +37.4 | formula_advantage |
| `icv1-fp42k6-t135044599113790442697853745259521825354-3ddeb182` | GF(p^6), p 42 bits | 239 | 135.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +34.2 | +37.4 | formula_advantage |
| `icv1-fp38k6-tm36212044362955483013997068193391274-ec247554` | GF(p^6), p 38 bits | 222 | 126.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +34.0 | +37.2 | formula_advantage |
| `icv1-fp40k6-tm932743203579365771526612219541923910-606ce65b` | GF(p^6), p 40 bits | 230 | 130.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +33.9 | +37.1 | formula_advantage |
| `icv1-fp36k6-t262812689477160887897309508847366-1543f6ae` | GF(p^6), p 36 bits | 214 | 122.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +33.8 | +37.0 | formula_advantage |
| `icv1-fp42k6-t150717022691811933176635677505680003042-306bfd23` | GF(p^6), p 42 bits | 237 | 134.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +33.6 | +36.9 | formula_advantage |
| `icv1-fp42k6-t96647665839480929326068241848685971530-1d8c267e` | GF(p^6), p 42 bits | 236 | 133.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +32.8 | +36.1 | formula_advantage |
| `icv1-fp42k6-t96647665839480929326068241848685971530-275b58bf` | GF(p^6), p 42 bits | 236 | 133.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +32.8 | +36.1 | formula_advantage |
| `icv1-fp40k6-t630557234967504570346093466502926814-da304588` | GF(p^6), p 40 bits | 228 | 129.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +32.8 | +36.1 | formula_advantage |
| `icv1-fp42k6-tm103424882089013968735935169172205724210-6f9e733b` | GF(p^6), p 42 bits | 236 | 133.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.7 | +32.7 | +36.0 | formula_advantage |
| `icv1-fp36k6-tm171229324292869415575893259583430-2d50770c` | GF(p^6), p 36 bits | 212 | 121.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +32.8 | +36.0 | formula_advantage |
| `icv1-fp36k6-tm171229324292869415575893259583430-566e3c75` | GF(p^6), p 36 bits | 212 | 121.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +32.8 | +36.0 | formula_advantage |
| `icv1-fp40k6-tm602411274482119859657366840582215118-0c890e3b` | GF(p^6), p 40 bits | 227 | 129.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +32.4 | +35.7 | formula_advantage |
| `icv1-fp40k6-tm602411274482119859657366840582215118-636f0036` | GF(p^6), p 40 bits | 227 | 129.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +32.4 | +35.7 | formula_advantage |
| `icv1-fp36k6-tm363622537811854311743487219063234-f3e78a13` | GF(p^6), p 36 bits | 211 | 121.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +32.4 | +35.6 | formula_advantage |
| `icv1-fp38k6-t2717978647929456783301116716814934-036a98dd` | GF(p^6), p 38 bits | 219 | 124.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +32.2 | +35.4 | formula_advantage |
| `icv1-fp36k6-t205253235936326958388310903933766-876f4b2b` | GF(p^6), p 36 bits | 211 | 120.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +32.1 | +35.3 | formula_advantage |
| `icv1-fp42k6-tm34100444908108954056602732293678277022-8c2c828f` | GF(p^6), p 42 bits | 234 | 132.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +31.9 | +35.2 | formula_advantage |
| `icv1-fp34k6-t589547474161273649282298788622-303031ca` | GF(p^6), p 34 bits | 202 | 116.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +31.8 | +35.0 | formula_advantage |
| `icv1-fp42k6-t70524409183775436746206207966618142326-971f838a` | GF(p^6), p 42 bits | 233 | 132.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +31.5 | +34.9 | formula_advantage |
| `icv1-fp38k6-tm4817483968765705745667162592660870-1d5d9451` | GF(p^6), p 38 bits | 217 | 124.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +31.5 | +34.8 | formula_advantage |
| `icv1-fp38k6-t23985362100570983542664784392022202-cad6e077` | GF(p^6), p 38 bits | 217 | 124.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +31.5 | +34.8 | formula_advantage |
| `icv1-fp38k6-t23985362100570983542664784392022202-d0d3391b` | GF(p^6), p 38 bits | 217 | 124.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +31.5 | +34.8 | formula_advantage |
| `icv1-fp40k6-tm308018148972052635901808025126390510-4bf8882f` | GF(p^6), p 40 bits | 225 | 128.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +31.3 | +34.6 | formula_advantage |
| `icv1-fp40k6-tm308018148972052635901808025126390510-b803bd28` | GF(p^6), p 40 bits | 225 | 128.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.6 | +31.3 | +34.6 | formula_advantage |
| `icv1-fp42k6-t20881010826321473227414487444218254566-8fc9a8b8` | GF(p^6), p 42 bits | 233 | 132.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +31.2 | +34.6 | formula_advantage |
| `icv1-fp38k6-t9784145319066463981943187610684410-d843344b` | GF(p^6), p 38 bits | 217 | 124.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +31.3 | +34.5 | formula_advantage |
| `icv1-fp38k6-t9784145319066463981943187610684410-f0319e81` | GF(p^6), p 38 bits | 217 | 124.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +31.3 | +34.5 | formula_advantage |
| `icv1-fp36k6-t478755799246825535986385953954326-cbe72568` | GF(p^6), p 36 bits | 209 | 119.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +31.2 | +34.5 | formula_advantage |
| `icv1-fp36k6-tm373757992317786448733652447629906-2e771f13` | GF(p^6), p 36 bits | 209 | 119.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +31.2 | +34.5 | formula_advantage |
| `icv1-fp40k6-tm2005756000445522048033269270111173070-7cb310e3` | GF(p^6), p 40 bits | 224 | 127.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +31.0 | +34.3 | formula_advantage |
| `icv1-fp42k6-tm4824240785322900275944834622558049178-cddeec0b` | GF(p^6), p 42 bits | 232 | 131.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +31.0 | +34.3 | formula_advantage |
| `icv1-fp38k6-t3051338254343228924269559952880010-a7da3d9f` | GF(p^6), p 38 bits | 216 | 123.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +31.0 | +34.3 | formula_advantage |
| `icv1-fp38k6-t3051338254343228924269559952880010-f893a3d7` | GF(p^6), p 38 bits | 216 | 123.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +31.0 | +34.3 | formula_advantage |
| `icv1-fp36k6-t236490246835594525795424392331978-547d5d69` | GF(p^6), p 36 bits | 208 | 119.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +31.0 | +34.3 | formula_advantage |
| `icv1-fp36k6-t236490246835594525795424392331978-56b2eecd` | GF(p^6), p 36 bits | 208 | 119.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +31.0 | +34.3 | formula_advantage |
| `icv1-fp36k6-t236490246835594525795424392331978-e8a4a42d` | GF(p^6), p 36 bits | 208 | 119.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +31.0 | +34.3 | formula_advantage |
| `icv1-fp40k6-tm696410845406644362359444856252571462-844b3ee1` | GF(p^6), p 40 bits | 224 | 127.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +31.0 | +34.3 | formula_advantage |
| `icv1-fp42k6-tm64576643174671653517239917461229492414-95e3d230` | GF(p^6), p 42 bits | 232 | 131.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +30.8 | +34.1 | formula_advantage |
| `icv1-fp40k6-tm20621316128534799417148382679580734-17fdbad9` | GF(p^6), p 40 bits | 224 | 127.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +30.7 | +34.1 | formula_advantage |
| `icv1-fp40k6-tm20621316128534799417148382679580734-e29ed42c` | GF(p^6), p 40 bits | 224 | 127.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +30.7 | +34.1 | formula_advantage |
| `icv1-fp34k6-t5399082817581527030251892957314-f9557687` | GF(p^6), p 34 bits | 200 | 115.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.8 | +34.0 | formula_advantage |
| `icv1-fp34k6-tm2894955473854795684072614329918-61a6fa7b` | GF(p^6), p 34 bits | 200 | 115.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.8 | +34.0 | formula_advantage |
| `icv1-fp34k6-tm2894955473854795684072614329918-f7898ef7` | GF(p^6), p 34 bits | 200 | 115.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.8 | +34.0 | formula_advantage |
| `icv1-fp34k6-tm698614099488561323391098718782-8369d03e` | GF(p^6), p 34 bits | 200 | 115.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.8 | +34.0 | formula_advantage |
| `icv1-fp36k6-t237789834949214153779567448808626-75bcb815` | GF(p^6), p 36 bits | 208 | 119.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +30.7 | +33.9 | formula_advantage |
| `icv1-fp40k6-t292033252515325002994064967366715738-9f30c710` | GF(p^6), p 40 bits | 223 | 127.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +30.5 | +33.9 | formula_advantage |
| `icv1-fp42k6-t8951448592536445311481466008335881530-408f25b2` | GF(p^6), p 42 bits | 231 | 131.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +30.3 | +33.6 | formula_advantage |
| `icv1-fp42k6-t8951448592536445311481466008335881530-f4d06a80` | GF(p^6), p 42 bits | 231 | 131.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +30.3 | +33.6 | formula_advantage |
| `icv1-fp38k6-t3632153308929655318308623828590042-5ee6c3c9` | GF(p^6), p 38 bits | 215 | 123.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +30.3 | +33.6 | formula_advantage |
| `icv1-fp38k6-t3632153308929655318308623828590042-d56fc627` | GF(p^6), p 38 bits | 215 | 123.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.5 | +30.3 | +33.6 | formula_advantage |
| `icv1-fp42k6-t21541954469755788973632133873550188182-025b4ab3` | GF(p^6), p 42 bits | 231 | 131.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +30.2 | +33.6 | formula_advantage |
| `icv1-fp34k6-t10034611511325963943441764652114-654f53d0` | GF(p^6), p 34 bits | 199 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.3 | +33.5 | formula_advantage |
| `icv1-fp34k6-t10034611511325963943441764652114-88cf6c79` | GF(p^6), p 34 bits | 199 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.3 | +33.5 | formula_advantage |
| `icv1-fp34k6-t3478139482010906749904310941586-2ca76047` | GF(p^6), p 34 bits | 199 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.3 | +33.5 | formula_advantage |
| `icv1-fp34k6-t3478139482010906749904310941586-da352b1d` | GF(p^6), p 34 bits | 199 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.3 | +33.5 | formula_advantage |
| `icv1-fp34k6-t6827801092898994424990489661394-05016338` | GF(p^6), p 34 bits | 199 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.3 | +33.5 | formula_advantage |
| `icv1-fp34k6-t6827801092898994424990489661394-a5f5309f` | GF(p^6), p 34 bits | 199 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +30.3 | +33.5 | formula_advantage |
| `icv1-fp38k6-tm3649156825249135653754775396351866-83d34b0a` | GF(p^6), p 38 bits | 215 | 123.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +30.2 | +33.5 | formula_advantage |
| `icv1-fp36k6-t466054673267443778344436740368950-5adf6593` | GF(p^6), p 36 bits | 207 | 118.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +30.1 | +33.4 | formula_advantage |
| `icv1-fp42k6-t90673943121076988551967377154538423126-5574e587` | GF(p^6), p 42 bits | 230 | 130.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +30.0 | +33.3 | formula_advantage |
| `icv1-fp40k6-t945757229327528023559427763154585934-fca3ec8e` | GF(p^6), p 40 bits | 222 | 126.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +29.9 | +33.3 | formula_advantage |
| `icv1-fp40k6-t501599052820971152177975224127972002-01e670a9` | GF(p^6), p 40 bits | 222 | 126.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +29.9 | +33.2 | formula_advantage |
| `icv1-fp40k6-t501599052820971152177975224127972002-98306107` | GF(p^6), p 40 bits | 222 | 126.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +29.9 | +33.2 | formula_advantage |
| `icv1-fp36k6-tm447302069326943222206923221257670-0baa51ac` | GF(p^6), p 36 bits | 206 | 118.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.4 | +29.9 | +33.1 | formula_advantage |
| `icv1-fp40k6-t879961038577238578773521469545766294-76e62059` | GF(p^6), p 40 bits | 222 | 126.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +29.8 | +33.1 | formula_advantage |
| `icv1-fp34k6-tm3042422134874163011242603665370-60c1cde4` | GF(p^6), p 34 bits | 198 | 114.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +29.9 | +33.1 | formula_advantage |
| `icv1-fp34k6-t460781918388327973060639259786-879f5ac0` | GF(p^6), p 34 bits | 198 | 114.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +29.7 | +33.0 | formula_advantage |
| `icv1-fp34k6-tm7853522135711894671428769205310-3d102641` | GF(p^6), p 34 bits | 198 | 114.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +29.7 | +32.9 | formula_advantage |
| `icv1-fp34k6-tm7853522135711894671428769205310-6b1e40ed` | GF(p^6), p 34 bits | 198 | 114.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +29.7 | +32.9 | formula_advantage |
| `icv1-fp38k6-tm23777983176963980758757184102360830-f50928fe` | GF(p^6), p 38 bits | 213 | 122.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +29.6 | +32.9 | formula_advantage |
| `icv1-fp40k6-tm1008121143751749127358300311755128398-adc75721` | GF(p^6), p 40 bits | 221 | 126.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +29.5 | +32.8 | formula_advantage |
| `icv1-fp40k6-tm1008121143751749127358300311755128398-bb447121` | GF(p^6), p 40 bits | 221 | 126.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +29.5 | +32.8 | formula_advantage |
| `icv1-fp40k6-tm1008121143751749127358300311755128398-e448bb1d` | GF(p^6), p 40 bits | 221 | 126.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +29.5 | +32.8 | formula_advantage |
| `icv1-fp38k6-tm8507701258781288579819180548798550-11d5fb25` | GF(p^6), p 38 bits | 213 | 122.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +29.5 | +32.8 | formula_advantage |
| `icv1-fp38k6-tm8507701258781288579819180548798550-e423712c` | GF(p^6), p 38 bits | 213 | 122.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +29.5 | +32.8 | formula_advantage |
| `icv1-fp38k6-t8507701258781288579819180548798550-16a3ee39` | GF(p^6), p 38 bits | 213 | 122.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +29.4 | +32.7 | formula_advantage |
| `icv1-fp34k6-tm3415120380487737561606789344386-35dcbc87` | GF(p^6), p 34 bits | 197 | 113.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.3 | +29.3 | +32.5 | formula_advantage |
| `icv1-fp36k6-t395211335361671592541987723781962-0f9e0267` | GF(p^6), p 36 bits | 205 | 117.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +29.3 | +32.5 | formula_advantage |
| `icv1-fp36k6-t395211335361671592541987723781962-8441ccca` | GF(p^6), p 36 bits | 205 | 117.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +29.3 | +32.5 | formula_advantage |
| `icv1-fp38k6-t9075672861665887541181946242606758-a20f520b` | GF(p^6), p 38 bits | 213 | 122.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +29.2 | +32.5 | formula_advantage |
| `icv1-fp34k6-tm18143208094393754906718265710-3b62a2e9` | GF(p^6), p 34 bits | 197 | 113.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +29.2 | +32.4 | formula_advantage |
| `icv1-fp34k6-t4910622335240915623989319497818-d780a6ab` | GF(p^6), p 34 bits | 197 | 113.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +29.1 | +32.4 | formula_advantage |
| `icv1-fp42k6-tm106214774323094524311483593733601712194-12b7d1a0` | GF(p^6), p 42 bits | 228 | 129.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +29.0 | +32.3 | formula_advantage |
| `icv1-fp42k6-t4601270379285945057251505653555594498-79f15d5b` | GF(p^6), p 42 bits | 228 | 129.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +28.9 | +32.3 | formula_advantage |
| `icv1-fp42k6-tm108109203991845710648600038206604674838-57193233` | GF(p^6), p 42 bits | 228 | 129.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +28.9 | +32.2 | formula_advantage |
| `icv1-fp36k6-t468168962925386983033169963896942-855bdcee` | GF(p^6), p 36 bits | 204 | 117.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +29.0 | +32.2 | formula_advantage |
| `icv1-fp32k6-tm63078555702759151532015341450-e49df7e1` | GF(p^6), p 32 bits | 188 | 109.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.2 | +29.0 | +32.2 | formula_advantage |
| `icv1-fp38k6-tm3518509357595230630226586171738438-71ec83a7` | GF(p^6), p 38 bits | 212 | 121.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +28.9 | +32.2 | formula_advantage |
| `icv1-fp38k6-tm3518509357595230630226586171738438-7a6b2941` | GF(p^6), p 38 bits | 212 | 121.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +28.9 | +32.2 | formula_advantage |
| `icv1-fp40k6-t2098435873855354521056169846574056274-8c1efd1f` | GF(p^6), p 40 bits | 220 | 125.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +28.8 | +32.1 | formula_advantage |
| `icv1-fp40k6-t2098435873855354521056169846574056274-c9471c32` | GF(p^6), p 40 bits | 220 | 125.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +28.8 | +32.1 | formula_advantage |
| `icv1-fp42k6-tm47625218402626040157351434367347524950-621637c4` | GF(p^6), p 42 bits | 228 | 129.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +28.7 | +32.1 | formula_advantage |
| `icv1-fp42k6-tm47625218402626040157351434367347524950-abe16089` | GF(p^6), p 42 bits | 228 | 129.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +28.7 | +32.1 | formula_advantage |
| `icv1-fp40k6-t30734519825362929505328835000769346-c47d0491` | GF(p^6), p 40 bits | 220 | 125.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +28.7 | +32.0 | formula_advantage |
| `icv1-fp40k6-t30734519825362929505328835000769346-e3e3af8a` | GF(p^6), p 40 bits | 220 | 125.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +28.7 | +32.0 | formula_advantage |
| `icv1-fp40k6-tm1871955401415393083467310545608280238-d3bbad43` | GF(p^6), p 40 bits | 220 | 125.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +28.6 | +32.0 | formula_advantage |
| `icv1-fp32k6-tm66356026174625047198089592338-309e821e` | GF(p^6), p 32 bits | 188 | 109.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +28.7 | +31.9 | formula_advantage |
| `icv1-fp38k6-tm39322260952504820046357871321031762-9b21cb30` | GF(p^6), p 38 bits | 211 | 121.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +28.5 | +31.9 | formula_advantage |
| `icv1-fp38k6-t18121727529971619060583448516174618-77ca168a` | GF(p^6), p 38 bits | 211 | 121.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +28.4 | +31.8 | formula_advantage |
| `icv1-fp38k6-tm8549936347722772953821962323339526-af332fd4` | GF(p^6), p 38 bits | 211 | 121.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +28.3 | +31.7 | formula_advantage |
| `icv1-fp34k6-tm1249817435402031484951825321466-5f7c638c` | GF(p^6), p 34 bits | 195 | 112.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +28.4 | +31.6 | formula_advantage |
| `icv1-fp42k6-tm48704942505053017245356289841019693990-7148ede6` | GF(p^6), p 42 bits | 226 | 129.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +28.1 | +31.4 | formula_advantage |
| `icv1-fp40k6-t2509949356857103406014342743482993022-55520c6c` | GF(p^6), p 40 bits | 218 | 124.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +28.1 | +31.4 | formula_advantage |
| `icv1-fp34k6-t8187850875526795465613921623090-32888faa` | GF(p^6), p 34 bits | 195 | 112.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +28.2 | +31.4 | formula_advantage |
| `icv1-fp34k6-t8187850875526795465613921623090-a6b3c922` | GF(p^6), p 34 bits | 195 | 112.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +28.2 | +31.4 | formula_advantage |
| `icv1-fp38k6-t22677321730373046260519994290263790-fc11da72` | GF(p^6), p 38 bits | 210 | 120.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +28.0 | +31.3 | formula_advantage |
| `icv1-fp36k6-t396134771930714512098166445913650-92a597c8` | GF(p^6), p 36 bits | 202 | 116.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +27.9 | +31.2 | formula_advantage |
| `icv1-fp42k6-tm113956595982069419360252349593637121274-0752a82e` | GF(p^6), p 42 bits | 226 | 128.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +27.8 | +31.2 | formula_advantage |
| `icv1-fp34k6-tm5449485995720001767339637263634-4b832ff0` | GF(p^6), p 34 bits | 194 | 112.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +27.9 | +31.2 | formula_advantage |
| `icv1-fp38k6-tm15263219706848806464205372905794706-c4823d22` | GF(p^6), p 38 bits | 210 | 120.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +27.8 | +31.2 | formula_advantage |
| `icv1-fp40k6-tm510675299706126771882283823057427614-42c8685e` | GF(p^6), p 40 bits | 218 | 124.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +27.8 | +31.1 | formula_advantage |
| `icv1-fp40k6-tm510675299706126771882283823057427614-b63f60a4` | GF(p^6), p 40 bits | 218 | 124.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +27.8 | +31.1 | formula_advantage |
| `icv1-fp38k6-t18125997103966405216085283473714470-55eb3283` | GF(p^6), p 38 bits | 210 | 120.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +27.7 | +31.0 | formula_advantage |
| `icv1-fp34k6-t5449485995720001767339637263634-babb41ef` | GF(p^6), p 34 bits | 194 | 112.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +27.6 | +30.9 | formula_advantage |
| `icv1-fp34k6-t5449485995720001767339637263634-eb68d886` | GF(p^6), p 34 bits | 194 | 112.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +27.6 | +30.9 | formula_advantage |
| `icv1-fp32k6-t78347245741357633901675975990-52f4b347` | GF(p^6), p 32 bits | 185 | 107.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +27.5 | +30.7 | formula_advantage |
| `icv1-fp32k6-t131394180832725690985430230374-b963d0e6` | GF(p^6), p 32 bits | 185 | 107.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +27.4 | +30.6 | formula_advantage |
| `icv1-fp36k6-tm77545060890920344875815742589942-5c95d140` | GF(p^6), p 36 bits | 201 | 115.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +27.3 | +30.6 | formula_advantage |
| `icv1-fp36k6-tm77545060890920344875815742589942-96213031` | GF(p^6), p 36 bits | 201 | 115.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +27.3 | +30.6 | formula_advantage |
| `icv1-fp42k6-tm69027207743865244667769899216784826010-4c253581` | GF(p^6), p 42 bits | 225 | 128.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.6 | +27.1 | +30.6 | formula_advantage |
| `icv1-fp34k6-tm997252301807352666642538383922-6f7e4fd0` | GF(p^6), p 34 bits | 193 | 111.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +27.2 | +30.5 | formula_advantage |
| `icv1-fp32k6-tm84464324299017033324873114198-5abe3ab6` | GF(p^6), p 32 bits | 185 | 107.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +27.2 | +30.4 | formula_advantage |
| `icv1-fp40k6-t875739874724797115605150218161230798-93e8dfe7` | GF(p^6), p 40 bits | 216 | 123.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +26.9 | +30.3 | formula_advantage |
| `icv1-fp42k6-tm94337795294083411234614731411795056966-f6a5b1ff` | GF(p^6), p 42 bits | 224 | 127.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +26.9 | +30.3 | formula_advantage |
| `icv1-fp40k6-tm712456293763396208466715541382314542-8913ae1e` | GF(p^6), p 40 bits | 216 | 123.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +26.8 | +30.2 | formula_advantage |
| `icv1-fp40k6-tm712456293763396208466715541382314542-bb413729` | GF(p^6), p 40 bits | 216 | 123.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.5 | +26.8 | +30.2 | formula_advantage |
| `icv1-fp32k6-tm14484914975528427710430752602-d2d56953` | GF(p^6), p 32 bits | 184 | 107.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +26.9 | +30.2 | formula_advantage |
| `icv1-fp32k6-t74332334055026805220953711114-c63f58a9` | GF(p^6), p 32 bits | 184 | 107.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +26.8 | +30.1 | formula_advantage |
| `icv1-fp42k6-tm90673943121076988551967377154538423126-7e441e3a` | GF(p^6), p 42 bits | 223 | 127.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +26.6 | +30.0 | formula_advantage |
| `icv1-fp42k6-tm70524409183775436746206207966618142326-12080389` | GF(p^6), p 42 bits | 223 | 127.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +26.5 | +29.9 | formula_advantage |
| `icv1-fp34k6-tm1023235065186341755971623767326-f8a9f9fe` | GF(p^6), p 34 bits | 192 | 111.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +26.6 | +29.9 | formula_advantage |
| `icv1-fp38k6-tm4333180390218985965181562993066582-241501bd` | GF(p^6), p 38 bits | 207 | 119.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +26.4 | +29.8 | formula_advantage |
| `icv1-fp38k6-tm4333180390218985965181562993066582-c1d643ad` | GF(p^6), p 38 bits | 207 | 119.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +26.4 | +29.8 | formula_advantage |
| `icv1-fp32k6-tm59621373577963920613733060330-8fd16810` | GF(p^6), p 32 bits | 183 | 106.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +26.5 | +29.8 | formula_advantage |
| `icv1-fp38k6-tm15135809390506173800996957253016490-b092ee07` | GF(p^6), p 38 bits | 207 | 119.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +26.4 | +29.8 | formula_advantage |
| `icv1-fp36k6-tm245692035398173207078702327659770-b2f1a843` | GF(p^6), p 36 bits | 199 | 114.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +26.3 | +29.7 | formula_advantage |
| `icv1-fp38k6-t11577268662598471418975454407051866-4f485d16` | GF(p^6), p 38 bits | 207 | 118.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +26.2 | +29.5 | formula_advantage |
| `icv1-fp38k6-t11577268662598471418975454407051866-cea6dea2` | GF(p^6), p 38 bits | 207 | 118.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +26.2 | +29.5 | formula_advantage |
| `icv1-fp34k6-tm2361525654024792672386354888502-03cd5615` | GF(p^6), p 34 bits | 191 | 110.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +26.2 | +29.5 | formula_advantage |
| `icv1-fp42k6-tm64034012845409394729220743976057956986-66ee9158` | GF(p^6), p 42 bits | 222 | 127.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +26.0 | +29.4 | formula_advantage |
| `icv1-fp40k6-tm1112861023040603223163073207108369490-7ee06db4` | GF(p^6), p 40 bits | 214 | 122.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +26.0 | +29.4 | formula_advantage |
| `icv1-fp42k6-tm81933216907660852315187676140185091926-37a1e378` | GF(p^6), p 42 bits | 222 | 126.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +25.9 | +29.4 | formula_advantage |
| `icv1-fp42k6-tm81933216907660852315187676140185091926-e6b4016d` | GF(p^6), p 42 bits | 222 | 126.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +25.9 | +29.4 | formula_advantage |
| `icv1-fp34k6-t1626025026650907714798762396866-0ab10c8c` | GF(p^6), p 34 bits | 190 | 110.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +26.1 | +29.4 | formula_advantage |
| `icv1-fp34k6-t1626025026650907714798762396866-d0fd8380` | GF(p^6), p 34 bits | 190 | 110.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +26.1 | +29.4 | formula_advantage |
| `icv1-fp36k6-tm153289125619373908359139010306838-2bc0f14a` | GF(p^6), p 36 bits | 198 | 114.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +26.0 | +29.3 | formula_advantage |
| `icv1-fp36k6-tm153289125619373908359139010306838-c636617a` | GF(p^6), p 36 bits | 198 | 114.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +26.0 | +29.3 | formula_advantage |
| `icv1-fp38k6-t13340058275007123218476566335469658-24757bcc` | GF(p^6), p 38 bits | 206 | 118.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +26.0 | +29.3 | formula_advantage |
| `icv1-fp38k6-t13340058275007123218476566335469658-9326aa73` | GF(p^6), p 38 bits | 206 | 118.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.4 | +26.0 | +29.3 | formula_advantage |
| `icv1-fp32k6-tm78347245741357633901675975990-52cc2d69` | GF(p^6), p 32 bits | 182 | 106.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +26.0 | +29.2 | formula_advantage |
| `icv1-fp32k6-tm78347245741357633901675975990-692a2cf8` | GF(p^6), p 32 bits | 182 | 106.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +26.0 | +29.2 | formula_advantage |
| `icv1-fp34k6-tm5395699078413387959179524439862-c906094b` | GF(p^6), p 34 bits | 190 | 110.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +25.9 | +29.2 | formula_advantage |
| `icv1-fp36k6-t109858165670618724298017225484562-4b937d50` | GF(p^6), p 36 bits | 198 | 114.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +25.8 | +29.2 | formula_advantage |
| `icv1-fp34k6-tm2052099295465299955628705112306-af3980c3` | GF(p^6), p 34 bits | 190 | 110.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +25.9 | +29.2 | formula_advantage |
| `icv1-fp36k6-t430636157326068713228412158502122-7a94c462` | GF(p^6), p 36 bits | 198 | 114.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +25.7 | +29.1 | formula_advantage |
| `icv1-fp36k6-t430636157326068713228412158502122-831cac03` | GF(p^6), p 36 bits | 198 | 114.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +25.7 | +29.1 | formula_advantage |
| `icv1-fp32k6-t15813541453371577054675899158-6f884100` | GF(p^6), p 32 bits | 182 | 106.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +25.8 | +29.1 | formula_advantage |
| `icv1-fp36k6-t447302069326943222206923221257670-67c91a1f` | GF(p^6), p 36 bits | 198 | 114.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.3 | +25.7 | +29.0 | formula_advantage |
| `icv1-fp42k6-tm2491251533173327737941972795111971366-526ba8ed` | GF(p^6), p 42 bits | 221 | 126.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +25.5 | +28.9 | formula_advantage |
| `icv1-fp42k6-tm2491251533173327737941972795111971366-67ed6d43` | GF(p^6), p 42 bits | 221 | 126.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +25.5 | +28.9 | formula_advantage |
| `icv1-fp34k6-tm3715044537733610511139125809350-1c24da01` | GF(p^6), p 34 bits | 190 | 110.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.2 | +25.6 | +28.9 | formula_advantage |
| `icv1-fp42k6-tm147495772641788131041904934022788125210-273c9b1c` | GF(p^6), p 42 bits | 221 | 126.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +25.4 | +28.8 | formula_advantage |
| `icv1-fp40k6-tm945757229327528023559427763154585934-20743ef0` | GF(p^6), p 40 bits | 213 | 122.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +25.4 | +28.8 | formula_advantage |
| `icv1-fp40k6-tm945757229327528023559427763154585934-b80f6867` | GF(p^6), p 40 bits | 213 | 122.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +25.4 | +28.8 | formula_advantage |
| `icv1-fp42k6-t28601062462356134170141359486915257214-a293aeb2` | GF(p^6), p 42 bits | 221 | 126.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +25.3 | +28.7 | formula_advantage |
| `icv1-fp40k6-tm1195900778100877146770415087895037282-14fa62e6` | GF(p^6), p 40 bits | 213 | 122.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +25.2 | +28.6 | formula_advantage |
| `icv1-fp38k6-t18784062855248225110074820977041502-f018f3b1` | GF(p^6), p 38 bits | 205 | 118.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +25.2 | +28.6 | formula_advantage |
| `icv1-fp36k6-t245692035398173207078702327659770-0b1c9513` | GF(p^6), p 36 bits | 197 | 113.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +25.2 | +28.6 | formula_advantage |
| `icv1-fp36k6-t245692035398173207078702327659770-82d068bd` | GF(p^6), p 36 bits | 197 | 113.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +25.2 | +28.6 | formula_advantage |
| `icv1-fp32k6-t32816149465264987946732403194-bea7e680` | GF(p^6), p 32 bits | 181 | 105.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +25.2 | +28.5 | formula_advantage |
| `icv1-fp32k6-t32816149465264987946732403194-df1a8832` | GF(p^6), p 32 bits | 181 | 105.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +25.2 | +28.5 | formula_advantage |
| `icv1-fp32k6-tm93288269984182716629364799690-7e9768db` | GF(p^6), p 32 bits | 181 | 105.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +25.1 | +28.4 | formula_advantage |
| `icv1-fp42k6-tm96457303662662337409351348340650293994-e4919e97` | GF(p^6), p 42 bits | 220 | 125.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +24.9 | +28.4 | formula_advantage |
| `icv1-fp32k6-t87449312845250641537531312086-813b1b3f` | GF(p^6), p 32 bits | 180 | 105.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.1 | +25.1 | +28.4 | formula_advantage |
| `icv1-fp42k6-t4824240785322900275944834622558049178-f2acb567` | GF(p^6), p 42 bits | 220 | 125.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +24.7 | +28.2 | formula_advantage |
| `icv1-fp42k6-t4824240785322900275944834622558049178-f648919c` | GF(p^6), p 42 bits | 220 | 125.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +24.7 | +28.2 | formula_advantage |
| `icv1-fp40k6-t290230737312338524932089875922384242-1cbbae2d` | GF(p^6), p 40 bits | 212 | 121.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +24.7 | +28.1 | formula_advantage |
| `icv1-fp34k6-tm8243999478733545754824242480078-a04a7eab` | GF(p^6), p 34 bits | 188 | 109.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +24.7 | +28.1 | formula_advantage |
| `icv1-fp34k6-tm8243999478733545754824242480078-f2c500a2` | GF(p^6), p 34 bits | 188 | 109.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +24.7 | +28.1 | formula_advantage |
| `icv1-fp32k6-tm16589822727447947892730390010-982ef622` | GF(p^6), p 32 bits | 180 | 105.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +24.8 | +28.0 | formula_advantage |
| `icv1-fp34k6-tm5173745615859079190160365337022-365bc820` | GF(p^6), p 34 bits | 188 | 109.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +24.6 | +28.0 | formula_advantage |
| `icv1-fp34k6-tm5173745615859079190160365337022-91eb85c6` | GF(p^6), p 34 bits | 188 | 109.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +24.6 | +28.0 | formula_advantage |
| `icv1-fp38k6-t3249817234800780839498103527633866-68a07786` | GF(p^6), p 38 bits | 203 | 117.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +24.6 | +28.0 | formula_advantage |
| `icv1-fp38k6-t37257395515011200599426178614034846-ece9b366` | GF(p^6), p 38 bits | 203 | 117.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +24.4 | +27.8 | formula_advantage |
| `icv1-fp32k6-t125938448560524821128630910618-aee6c766` | GF(p^6), p 32 bits | 179 | 104.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +24.5 | +27.8 | formula_advantage |
| `icv1-fp32k6-t125938448560524821128630910618-d8da3956` | GF(p^6), p 32 bits | 179 | 104.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +24.5 | +27.8 | formula_advantage |
| `icv1-fp40k6-t1238296160463272838724820471522147074-9ee73ff2` | GF(p^6), p 40 bits | 211 | 121.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +24.3 | +27.7 | formula_advantage |
| `icv1-fp40k6-t1238296160463272838724820471522147074-d5b6d994` | GF(p^6), p 40 bits | 211 | 121.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +24.3 | +27.7 | formula_advantage |
| `icv1-fp42k6-tm94425043260925418064943765136461443902-2954d8b9` | GF(p^6), p 42 bits | 219 | 125.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +24.2 | +27.7 | formula_advantage |
| `icv1-fp40k6-tm1248579780941707831635671970573221790-b59e9999` | GF(p^6), p 40 bits | 211 | 121.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +24.2 | +27.7 | formula_advantage |
| `icv1-fp40k6-tm1248579780941707831635671970573221790-f439b6f4` | GF(p^6), p 40 bits | 211 | 121.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +24.2 | +27.7 | formula_advantage |
| `icv1-fp40k6-t448236941744719717984179615361417866-d8e165e9` | GF(p^6), p 40 bits | 211 | 121.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +24.2 | +27.6 | formula_advantage |
| `icv1-fp38k6-t10369598231134931509632480969249882-6ba214f3` | GF(p^6), p 38 bits | 203 | 116.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +24.2 | +27.6 | formula_advantage |
| `icv1-fp38k6-t10369598231134931509632480969249882-daf54f11` | GF(p^6), p 38 bits | 203 | 116.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +24.2 | +27.6 | formula_advantage |
| `icv1-fp40k6-tm1927832409848495291239213980680019890-aa92dff1` | GF(p^6), p 40 bits | 211 | 121.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +24.1 | +27.6 | formula_advantage |
| `icv1-fp38k6-tm21899198471095324845165564550270322-1bfc805e` | GF(p^6), p 38 bits | 203 | 116.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +24.1 | +27.5 | formula_advantage |
| `icv1-fp42k6-tm164816127623881732294837773186127712426-dabbf73d` | GF(p^6), p 42 bits | 218 | 125.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +24.0 | +27.5 | formula_advantage |
| `icv1-fp34k6-tm1768157185091978969274182572622-d750d67c` | GF(p^6), p 34 bits | 187 | 108.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +24.1 | +27.4 | formula_advantage |
| `icv1-fp34k6-t2778058223101701185177335516194-2d20f60b` | GF(p^6), p 34 bits | 186 | 108.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +24.0 | +27.3 | formula_advantage |
| `icv1-fp34k6-t2778058223101701185177335516194-f5dd9176` | GF(p^6), p 34 bits | 186 | 108.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +24.0 | +27.3 | formula_advantage |
| `icv1-fp32k6-t135366469999176815019213964042-38f5cf73` | GF(p^6), p 32 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +24.0 | +27.3 | formula_advantage |
| `icv1-fp32k6-tm25593149993636054719270006726-58cc5d39` | GF(p^6), p 32 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +24.0 | +27.3 | formula_advantage |
| `icv1-fp36k6-t416059758384670072236771416265090-8443330c` | GF(p^6), p 36 bits | 194 | 112.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +23.8 | +27.2 | formula_advantage |
| `icv1-fp34k6-tm8343469462759645670504383165214-3fe77ea8` | GF(p^6), p 34 bits | 186 | 108.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +23.8 | +27.2 | formula_advantage |
| `icv1-fp34k6-tm8343469462759645670504383165214-88ecdbd5` | GF(p^6), p 34 bits | 186 | 108.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +23.8 | +27.2 | formula_advantage |
| `icv1-fp32k6-tm37810418333563058468914351474-68b7b0bb` | GF(p^6), p 32 bits | 178 | 104.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +23.8 | +27.1 | formula_advantage |
| `icv1-fp36k6-tm423060820223822621613427538943606-b60b16a2` | GF(p^6), p 36 bits | 194 | 112.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +23.7 | +27.1 | formula_advantage |
| `icv1-fp36k6-tm423060820223822621613427538943606-ce453b1d` | GF(p^6), p 36 bits | 194 | 112.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +23.7 | +27.1 | formula_advantage |
| `icv1-fp34k6-t7882711358008644359153620611790-a4e9beb2` | GF(p^6), p 34 bits | 186 | 108.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +23.7 | +27.0 | formula_advantage |
| `icv1-fp38k6-t33653313885195831439390070846015878-6d238eb6` | GF(p^6), p 38 bits | 201 | 116.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +23.4 | +26.8 | formula_advantage |
| `icv1-fp34k6-t1852025694741637258036688563202-4246dd27` | GF(p^6), p 34 bits | 185 | 107.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +23.5 | +26.8 | formula_advantage |
| `icv1-fp34k6-t1852025694741637258036688563202-50976e59` | GF(p^6), p 34 bits | 185 | 107.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +23.5 | +26.8 | formula_advantage |
| `icv1-fp36k6-t56961645353625724903712231054462-28a5b421` | GF(p^6), p 36 bits | 193 | 112.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +23.4 | +26.8 | formula_advantage |
| `icv1-fp36k6-t400576880020473356175501275633114-a7137fc6` | GF(p^6), p 36 bits | 193 | 111.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +23.3 | +26.7 | formula_advantage |
| `icv1-fp36k6-t349583752751787938470935204879638-84988901` | GF(p^6), p 36 bits | 193 | 111.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +23.3 | +26.7 | formula_advantage |
| `icv1-fp33k6-tm24932403969225860243224553758-8c249042` | GF(p^6), p 32 bits | 177 | 103.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +23.4 | +26.7 | formula_advantage |
| `icv1-fp33k6-tm24932403969225860243224553758-92482ac0` | GF(p^6), p 32 bits | 177 | 103.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +23.4 | +26.7 | formula_advantage |
| `icv1-fp36k6-t170593691427670977961111241816650-0f07b579` | GF(p^6), p 36 bits | 193 | 111.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +23.2 | +26.7 | formula_advantage |
| `icv1-fp38k6-t11470009159992702495635638686870810-686bc72e` | GF(p^6), p 38 bits | 201 | 115.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +23.2 | +26.6 | formula_advantage |
| `icv1-fp38k6-t35129226029225194675644221938173586-21cb7109` | GF(p^6), p 38 bits | 201 | 115.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +23.1 | +26.5 | formula_advantage |
| `icv1-fp40k6-tm1495053359913849976434210433711990670-7796c4f5` | GF(p^6), p 40 bits | 208 | 119.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +23.0 | +26.5 | formula_advantage |
| `icv1-fp36k6-tm204722971118342323583719762246230-c3892924` | GF(p^6), p 36 bits | 193 | 111.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +23.1 | +26.5 | formula_advantage |
| `icv1-fp42k6-t127882077426822819848891979507313072634-2c367235` | GF(p^6), p 42 bits | 216 | 123.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +23.0 | +26.5 | formula_advantage |
| `icv1-fp42k6-t127882077426822819848891979507313072634-e22938e7` | GF(p^6), p 42 bits | 216 | 123.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +23.0 | +26.5 | formula_advantage |
| `icv1-fp34k6-t7165476592591050354659853469990-bd90687d` | GF(p^6), p 34 bits | 185 | 107.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +23.1 | +26.5 | formula_advantage |
| `icv1-fp40k6-tm201274051287991827256524896036165882-1073078b` | GF(p^6), p 40 bits | 208 | 119.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +22.9 | +26.4 | formula_advantage |
| `icv1-fp42k6-t19461047747900396722106726129360388450-cf86fac4` | GF(p^6), p 42 bits | 216 | 123.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +22.7 | +26.2 | formula_advantage |
| `icv1-fp40k6-tm1875261934287626114672193986136586826-d432179f` | GF(p^6), p 40 bits | 208 | 119.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +22.6 | +26.1 | formula_advantage |
| `icv1-fp32k6-t18289037840015464640675956522-62cce425` | GF(p^6), p 32 bits | 176 | 103.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +22.7 | +26.1 | formula_advantage |
| `icv1-fp32k6-t18289037840015464640675956522-a27ec798` | GF(p^6), p 32 bits | 176 | 103.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +22.7 | +26.1 | formula_advantage |
| `icv1-fp38k6-t36212044362955483013997068193391274-616aaceb` | GF(p^6), p 38 bits | 199 | 115.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +22.6 | +26.0 | formula_advantage |
| `icv1-fp38k6-t15135809390506173800996957253016490-49160f02` | GF(p^6), p 38 bits | 199 | 115.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +22.4 | +25.9 | formula_advantage |
| `icv1-fp38k6-t15135809390506173800996957253016490-8cc67011` | GF(p^6), p 38 bits | 199 | 115.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +22.4 | +25.9 | formula_advantage |
| `icv1-fp34k6-tm8640604617662406951608052462898-c4095f6a` | GF(p^6), p 34 bits | 183 | 107.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +22.5 | +25.9 | formula_advantage |
| `icv1-fp42k6-t4811780584012495637473434679770450170-28ed7233` | GF(p^6), p 42 bits | 215 | 123.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +22.4 | +25.9 | formula_advantage |
| `icv1-fp42k6-t4811780584012495637473434679770450170-523f8289` | GF(p^6), p 42 bits | 215 | 123.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +22.4 | +25.9 | formula_advantage |
| `icv1-fp36k6-tm116719443578082754024866611148310-11a3cb72` | GF(p^6), p 36 bits | 191 | 111.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +22.4 | +25.9 | formula_advantage |
| `icv1-fp36k6-tm116719443578082754024866611148310-19795633` | GF(p^6), p 36 bits | 191 | 111.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +22.4 | +25.9 | formula_advantage |
| `icv1-fp42k6-t31828309486290994963515878471386978850-3e66b1c7` | GF(p^6), p 42 bits | 215 | 123.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.5 | +22.3 | +25.8 | formula_advantage |
| `icv1-fp42k6-tm26651597370808610731249890859521747770-42d67a40` | GF(p^6), p 42 bits | 215 | 123.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +22.2 | +25.8 | formula_advantage |
| `icv1-fp36k6-t559600199465457298174643799810698-5eaee8b2` | GF(p^6), p 36 bits | 191 | 110.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +22.2 | +25.6 | formula_advantage |
| `icv1-fp38k6-tm21778159553746265407774699488272586-b3e2c8a0` | GF(p^6), p 38 bits | 198 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +22.0 | +25.5 | formula_advantage |
| `icv1-fp36k6-tm333499376444505757633055028910694-ba9c94e3` | GF(p^6), p 36 bits | 190 | 110.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +22.0 | +25.5 | formula_advantage |
| `icv1-fp36k6-tm333499376444505757633055028910694-edf4a289` | GF(p^6), p 36 bits | 190 | 110.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +22.0 | +25.5 | formula_advantage |
| `icv1-fp38k6-t31116608428445849544095091174863146-5c8c30c6` | GF(p^6), p 38 bits | 198 | 114.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +22.0 | +25.4 | formula_advantage |
| `icv1-fp38k6-t31116608428445849544095091174863146-e17377c9` | GF(p^6), p 38 bits | 198 | 114.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +22.0 | +25.4 | formula_advantage |
| `icv1-fp36k6-tm61402725279656295982520444550682-c413d15b` | GF(p^6), p 36 bits | 190 | 110.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +22.0 | +25.4 | formula_advantage |
| `icv1-fp42k6-t6129872371660289382057012330257869986-e88bc63e` | GF(p^6), p 42 bits | 214 | 122.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +21.9 | +25.4 | formula_advantage |
| `icv1-fp40k6-tm665405924970026162701491545965324990-0a76b04f` | GF(p^6), p 40 bits | 206 | 118.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.4 | +21.8 | +25.3 | formula_advantage |
| `icv1-fp32k6-t56939309847052024180839379882-a7505d96` | GF(p^6), p 32 bits | 174 | 102.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +21.9 | +25.3 | formula_advantage |
| `icv1-fp34k6-t2052099295465299955628705112306-fda7d398` | GF(p^6), p 34 bits | 182 | 106.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +21.9 | +25.3 | formula_advantage |
| `icv1-fp38k6-tm15309988203816246634582542258998526-9a49339d` | GF(p^6), p 38 bits | 198 | 114.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.3 | +21.8 | +25.2 | formula_advantage |
| `icv1-fp32k6-tm36397365622551860598626024294-07705586` | GF(p^6), p 32 bits | 174 | 102.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +21.9 | +25.2 | formula_advantage |
| `icv1-fp32k6-tm36397365622551860598626024294-b6453a99` | GF(p^6), p 32 bits | 174 | 102.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +21.9 | +25.2 | formula_advantage |
| `icv1-fp40k6-tm494309492534452699850592604986578094-51cd3a28` | GF(p^6), p 40 bits | 206 | 118.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +21.7 | +25.2 | formula_advantage |
| `icv1-fp40k6-t1927832409848495291239213980680019890-0f520c8d` | GF(p^6), p 40 bits | 206 | 118.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +21.6 | +25.1 | formula_advantage |
| `icv1-fp40k6-tm630557234967504570346093466502926814-47e35d62` | GF(p^6), p 40 bits | 205 | 118.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +21.5 | +25.0 | formula_advantage |
| `icv1-fp40k6-tm630557234967504570346093466502926814-9067537f` | GF(p^6), p 40 bits | 205 | 118.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +21.5 | +25.0 | formula_advantage |
| `icv1-fp36k6-t61402725279656295982520444550682-136947b4` | GF(p^6), p 36 bits | 189 | 110.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +21.6 | +25.0 | formula_advantage |
| `icv1-fp36k6-t61402725279656295982520444550682-b08cc903` | GF(p^6), p 36 bits | 189 | 110.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +21.6 | +25.0 | formula_advantage |
| `icv1-fp42k6-t18845429001969576848180501317741958110-db075bbc` | GF(p^6), p 42 bits | 213 | 122.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +21.5 | +25.0 | formula_advantage |
| `icv1-fp36k6-t151269242012529176447176043075558-615c30d2` | GF(p^6), p 36 bits | 189 | 110.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.2 | +21.6 | +25.0 | formula_advantage |
| `icv1-fp40k6-tm2509949356857103406014342743482993022-7abfae23` | GF(p^6), p 40 bits | 205 | 118.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +21.4 | +24.9 | formula_advantage |
| `icv1-fp40k6-tm2509949356857103406014342743482993022-f51c41c1` | GF(p^6), p 40 bits | 205 | 118.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +21.4 | +24.9 | formula_advantage |
| `icv1-fp42k6-t147495772641788131041904934022788125210-a86c41a7` | GF(p^6), p 42 bits | 213 | 122.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +21.3 | +24.8 | formula_advantage |
| `icv1-fp42k6-t147495772641788131041904934022788125210-cf9ae874` | GF(p^6), p 42 bits | 213 | 122.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +21.3 | +24.8 | formula_advantage |
| `icv1-fp38k6-t12628575110007103759117203599212102-8f67fc26` | GF(p^6), p 38 bits | 197 | 114.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +21.3 | +24.8 | formula_advantage |
| `icv1-fp42k6-t75417726566101703640201649601320175710-f964f935` | GF(p^6), p 42 bits | 213 | 122.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +21.1 | +24.7 | formula_advantage |
| `icv1-fp38k6-t25043404913375920520837933282540958-d351a320` | GF(p^6), p 38 bits | 196 | 113.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +21.1 | +24.6 | formula_advantage |
| `icv1-fp34k6-t4867682541869492820951426833146-cbdb7f1e` | GF(p^6), p 34 bits | 181 | 105.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.1 | +21.1 | +24.6 | formula_advantage |
| `icv1-fp42k6-t16061877539138736143564201885824856426-d5bb8b43` | GF(p^6), p 42 bits | 212 | 122.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +21.0 | +24.5 | formula_advantage |
| `icv1-fp32k6-tm83035224272727947564063586618-95c4c0fc` | GF(p^6), p 32 bits | 172 | 101.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 77.0 | +21.0 | +24.4 | formula_advantage |
| `icv1-fp42k6-t140923858286471124179460900847362548778-5d4e3592` | GF(p^6), p 42 bits | 212 | 121.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +20.7 | +24.3 | formula_advantage |
| `icv1-fp32k6-tm141867637087080037808722483910-83273f42` | GF(p^6), p 32 bits | 172 | 101.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +20.8 | +24.2 | formula_advantage |
| `icv1-fp32k6-tm141867637087080037808722483910-aa2cf7f1` | GF(p^6), p 32 bits | 172 | 101.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +20.8 | +24.2 | formula_advantage |
| `icv1-fp36k6-t434632891845129734969098332602618-42ac0f33` | GF(p^6), p 36 bits | 188 | 109.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.1 | +20.7 | +24.2 | formula_advantage |
| `icv1-fp40k6-t1178780120850042221354987923466218202-a3cd03c0` | GF(p^6), p 40 bits | 203 | 117.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +20.6 | +24.1 | formula_advantage |
| `icv1-fp36k6-t266347467794506259203439006117222-c57ab59d` | GF(p^6), p 36 bits | 188 | 109.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.1 | +20.6 | +24.1 | formula_advantage |
| `icv1-fp32k6-t33054116898269288753305348542-8b209c9e` | GF(p^6), p 32 bits | 172 | 101.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +20.6 | +24.0 | formula_advantage |
| `icv1-fp34k6-t5173745615859079190160365337022-7eb64057` | GF(p^6), p 34 bits | 179 | 105.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +20.5 | +24.0 | formula_advantage |
| `icv1-fp38k6-tm18125997103966405216085283473714470-4f66d150` | GF(p^6), p 38 bits | 195 | 113.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +20.4 | +23.9 | formula_advantage |
| `icv1-fp38k6-tm18125997103966405216085283473714470-843c4598` | GF(p^6), p 38 bits | 195 | 113.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +20.4 | +23.9 | formula_advantage |
| `icv1-fp38k6-tm22151669180698970899803872554108238-3bc646c6` | GF(p^6), p 38 bits | 195 | 113.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +20.4 | +23.9 | formula_advantage |
| `icv1-fp36k6-t204722971118342323583719762246230-0d537e88` | GF(p^6), p 36 bits | 187 | 108.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.1 | +20.3 | +23.8 | formula_advantage |
| `icv1-fp40k6-tm1941187569253126322303441098237856698-e88a9ceb` | GF(p^6), p 40 bits | 203 | 117.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +20.2 | +23.7 | formula_advantage |
| `icv1-fp42k6-t96457303662662337409351348340650293994-23e39bc1` | GF(p^6), p 42 bits | 211 | 121.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +20.1 | +23.7 | formula_advantage |
| `icv1-fp42k6-t96457303662662337409351348340650293994-d9af91f0` | GF(p^6), p 42 bits | 211 | 121.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +20.1 | +23.7 | formula_advantage |
| `icv1-fp40k6-tm455440728487123980231477854918140510-f746edd1` | GF(p^6), p 40 bits | 203 | 117.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +20.1 | +23.6 | formula_advantage |
| `icv1-fp32k6-tm49818530323077121909575684326-437a7d02` | GF(p^6), p 32 bits | 171 | 100.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +20.2 | +23.6 | formula_advantage |
| `icv1-fp34k6-tm1626025026650907714798762396866-f4aafaf9` | GF(p^6), p 34 bits | 178 | 104.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +20.0 | +23.4 | formula_advantage |
| `icv1-fp32k6-t30810004660830079984930045046-ac852a38` | GF(p^6), p 32 bits | 170 | 100.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +20.0 | +23.4 | formula_advantage |
| `icv1-fp34k6-t1174294105908176066788580091062-cade2114` | GF(p^6), p 34 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +19.8 | +23.3 | formula_advantage |
| `icv1-fp34k6-t222467787949046006723433574990-ad47291b` | GF(p^6), p 34 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +19.8 | +23.3 | formula_advantage |
| `icv1-fp34k6-t274801603712586561732682558770-3175e34f` | GF(p^6), p 34 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +19.8 | +23.3 | formula_advantage |
| `icv1-fp34k6-tm7882711358008644359153620611790-2745225a` | GF(p^6), p 34 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +19.8 | +23.2 | formula_advantage |
| `icv1-fp34k6-tm7882711358008644359153620611790-d0c57ce5` | GF(p^6), p 34 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +19.8 | +23.2 | formula_advantage |
| `icv1-fp42k6-t123044564208745112288712406413410987534-24c80bbf` | GF(p^6), p 42 bits | 210 | 120.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +19.7 | +23.2 | formula_advantage |
| `icv1-fp40k6-t272141120984694834845539129520103906-8bc053d2` | GF(p^6), p 40 bits | 202 | 116.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +19.7 | +23.2 | formula_advantage |
| `icv1-fp40k6-t272141120984694834845539129520103906-e3d362a2` | GF(p^6), p 40 bits | 202 | 116.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +19.7 | +23.2 | formula_advantage |
| `icv1-fp38k6-t20384121038916434198934221667625930-01417586` | GF(p^6), p 38 bits | 194 | 112.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +19.6 | +23.1 | formula_advantage |
| `icv1-fp38k6-t20384121038916434198934221667625930-0dfd2572` | GF(p^6), p 38 bits | 194 | 112.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +19.6 | +23.1 | formula_advantage |
| `icv1-fp32k6-t139745304813034977667967079314-295dd669` | GF(p^6), p 32 bits | 169 | 99.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +19.5 | +23.0 | formula_advantage |
| `icv1-fp40k6-tm1883083775382747833993667748916468750-5cc05c11` | GF(p^6), p 40 bits | 201 | 116.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +19.3 | +22.8 | formula_advantage |
| `icv1-fp40k6-tm1883083775382747833993667748916468750-64246701` | GF(p^6), p 40 bits | 201 | 116.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +19.3 | +22.8 | formula_advantage |
| `icv1-fp34k6-tm5017136059121885029961175382086-c61b1bd3` | GF(p^6), p 34 bits | 177 | 103.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +19.3 | +22.8 | formula_advantage |
| `icv1-fp36k6-t490498038359306366180694999614966-da7dc0e5` | GF(p^6), p 36 bits | 185 | 107.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.1 | +19.2 | +22.7 | formula_advantage |
| `icv1-fp38k6-tm25703862405918691689370241743508566-f13c7f7d` | GF(p^6), p 38 bits | 193 | 111.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +19.1 | +22.7 | formula_advantage |
| `icv1-fp42k6-t85965374489496676757941558323912879818-c784c69f` | GF(p^6), p 42 bits | 208 | 120.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +19.0 | +22.6 | formula_advantage |
| `icv1-fp42k6-t85965374489496676757941558323912879818-c8528336` | GF(p^6), p 42 bits | 208 | 120.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +19.0 | +22.6 | formula_advantage |
| `icv1-fp32k6-t29242045006735922381596573882-8a1fca43` | GF(p^6), p 32 bits | 168 | 99.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +19.1 | +22.5 | formula_advantage |
| `icv1-fp42k6-tm117264918480432704958797206118676183454-fed7e2bc` | GF(p^6), p 42 bits | 208 | 119.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +18.9 | +22.5 | formula_advantage |
| `icv1-fp40k6-t1079691786106281078201479722463943470-0942650d` | GF(p^6), p 40 bits | 200 | 115.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +18.9 | +22.4 | formula_advantage |
| `icv1-fp38k6-tm11577268662598471418975454407051866-33edb44c` | GF(p^6), p 38 bits | 192 | 111.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +18.9 | +22.4 | formula_advantage |
| `icv1-fp40k6-tm1254970134597607642355604557049934798-12c7c3ab` | GF(p^6), p 40 bits | 200 | 115.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +18.8 | +22.4 | formula_advantage |
| `icv1-fp40k6-tm1254970134597607642355604557049934798-298c3db1` | GF(p^6), p 40 bits | 200 | 115.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +18.8 | +22.4 | formula_advantage |
| `icv1-fp36k6-tm255886181995749826340143739575682-53d2fda7` | GF(p^6), p 36 bits | 184 | 107.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.1 | +18.8 | +22.3 | formula_advantage |
| `icv1-fp32k6-t63078555702759151532015341450-80b9eed1` | GF(p^6), p 32 bits | 168 | 99.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +18.6 | +22.1 | formula_advantage |
| `icv1-fp34k6-tm9510370307821810800508195287014-c12d9890` | GF(p^6), p 34 bits | 175 | 103.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +18.6 | +22.1 | formula_advantage |
| `icv1-fp32k6-tm18289037840015464640675956522-9c95fc0a` | GF(p^6), p 32 bits | 167 | 98.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +18.5 | +22.0 | formula_advantage |
| `icv1-fp38k6-t6767891398895265965022010568115226-8e55eef0` | GF(p^6), p 38 bits | 191 | 111.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +18.4 | +22.0 | formula_advantage |
| `icv1-fp32k6-t8210865583189936379968334058-7d679ba4` | GF(p^6), p 32 bits | 167 | 98.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +18.5 | +22.0 | formula_advantage |
| `icv1-fp40k6-tm980612957738896330886104157610531838-54f5b89c` | GF(p^6), p 40 bits | 199 | 115.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +18.3 | +21.9 | formula_advantage |
| `icv1-fp40k6-tm980612957738896330886104157610531838-874a5c0f` | GF(p^6), p 40 bits | 199 | 115.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +18.3 | +21.9 | formula_advantage |
| `icv1-fp32k6-tm83045631604073078514661537106-60bbd785` | GF(p^6), p 32 bits | 167 | 98.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +18.4 | +21.9 | formula_advantage |
| `icv1-fp34k6-t8640604617662406951608052462898-22662afe` | GF(p^6), p 34 bits | 175 | 102.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +18.4 | +21.9 | formula_advantage |
| `icv1-fp42k6-t90676903788957252138290740157732070074-4a026d10` | GF(p^6), p 42 bits | 207 | 119.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +18.2 | +21.8 | formula_advantage |
| `icv1-fp32k6-t26895410567024103251795854646-5e67e877` | GF(p^6), p 32 bits | 167 | 98.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +18.3 | +21.8 | formula_advantage |
| `icv1-fp34k6-t4889968707445081258109231685378-63750ab3` | GF(p^6), p 34 bits | 175 | 102.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +18.2 | +21.7 | formula_advantage |
| `icv1-fp42k6-t26651597370808610731249890859521747770-60b8120b` | GF(p^6), p 42 bits | 206 | 119.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.4 | +18.0 | +21.6 | formula_advantage |
| `icv1-fp34k6-tm6232445055316855885358383888106-c63c85b1` | GF(p^6), p 34 bits | 175 | 102.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +18.1 | +21.6 | formula_advantage |
| `icv1-fp38k6-tm12045386832711253957245317395179814-623e3c26` | GF(p^6), p 38 bits | 190 | 110.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +18.0 | +21.6 | formula_advantage |
| `icv1-fp38k6-tm12045386832711253957245317395179814-99d69e5a` | GF(p^6), p 38 bits | 190 | 110.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +18.0 | +21.6 | formula_advantage |
| `icv1-fp38k6-tm12045386832711253957245317395179814-d7a2cd3e` | GF(p^6), p 38 bits | 190 | 110.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +18.0 | +21.6 | formula_advantage |
| `icv1-fp32k6-t127812454930946669422317632054-eb59bb0c` | GF(p^6), p 32 bits | 167 | 98.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +18.1 | +21.6 | formula_advantage |
| `icv1-fp40k6-tm78757769114834478014016497882486494-7aebc9da` | GF(p^6), p 40 bits | 198 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +17.9 | +21.5 | formula_advantage |
| `icv1-fp40k6-tm78757769114834478014016497882486494-9c1df671` | GF(p^6), p 40 bits | 198 | 114.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +17.9 | +21.5 | formula_advantage |
| `icv1-fp32k6-t8516682081799136568269031802-2d61cd0c` | GF(p^6), p 32 bits | 166 | 98.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +18.0 | +21.5 | formula_advantage |
| `icv1-fp32k6-t8516682081799136568269031802-d00cf5ff` | GF(p^6), p 32 bits | 166 | 98.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +18.0 | +21.5 | formula_advantage |
| `icv1-fp34k6-tm1060509355543416520796049580354-af99a71e` | GF(p^6), p 34 bits | 174 | 102.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +17.9 | +21.4 | formula_advantage |
| `icv1-fp32k6-tm8210865583189936379968334058-440ae268` | GF(p^6), p 32 bits | 166 | 98.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +17.9 | +21.4 | formula_advantage |
| `icv1-fp34k6-tm4224784695443668058230807120510-65b8ff98` | GF(p^6), p 34 bits | 174 | 102.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +17.9 | +21.4 | formula_advantage |
| `icv1-fp34k6-tm6827801092898994424990489661394-02ee041c` | GF(p^6), p 34 bits | 174 | 102.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +17.8 | +21.3 | formula_advantage |
| `icv1-fp40k6-t346133854454224542686020282566086050-bfe3c820` | GF(p^6), p 40 bits | 198 | 114.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.3 | +17.7 | +21.3 | formula_advantage |
| `icv1-fp36k6-tm108019297965626695735775020187222-6845c54b` | GF(p^6), p 36 bits | 182 | 106.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.1 | +17.8 | +21.3 | formula_advantage |
| `icv1-fp38k6-tm12728477339114286346496772952805918-d946d97c` | GF(p^6), p 38 bits | 190 | 110.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +17.7 | +21.2 | formula_advantage |
| `icv1-fp38k6-t9652483448788076169284986507797854-b32bf5bd` | GF(p^6), p 38 bits | 189 | 110.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.2 | +17.6 | +21.1 | formula_advantage |
| `icv1-fp36k6-t350990977962642348083811637516458-5d8adf3f` | GF(p^6), p 36 bits | 181 | 106.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.1 | +17.5 | +21.0 | formula_advantage |
| `icv1-fp36k6-t350990977962642348083811637516458-bc72a33f` | GF(p^6), p 36 bits | 181 | 106.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.1 | +17.5 | +21.0 | formula_advantage |
| `icv1-fp34k6-tm5054121035902274520965406734810-545d6854` | GF(p^6), p 34 bits | 173 | 101.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 81.0 | +17.5 | +21.0 | formula_advantage |
| `icv1-fp32k6-t46110035577588264290417783898-aec1528a` | GF(p^6), p 32 bits | 165 | 97.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +17.5 | +21.0 | formula_advantage |
| `icv1-fp32k6-tm127993243901090866891695182754-e2ada1e5` | GF(p^6), p 32 bits | 165 | 97.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.9 | +17.5 | +20.9 | formula_advantage |
| `icv1-fp40k6-t1746682707660615888470782891457259390-69feba67` | GF(p^6), p 40 bits | 196 | 113.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +17.0 | +20.6 | formula_advantage |
| `icv1-fp36k6-tm115876558483567622403462425082082-3fca665b` | GF(p^6), p 36 bits | 180 | 105.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +17.0 | +20.6 | formula_advantage |
| `icv1-fp36k6-t457754664664390214216967658012570-0a4faca1` | GF(p^6), p 36 bits | 180 | 105.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +16.8 | +20.4 | formula_advantage |
| `icv1-fp36k6-t457754664664390214216967658012570-2dad9111` | GF(p^6), p 36 bits | 180 | 105.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +16.8 | +20.4 | formula_advantage |
| `icv1-fp32k6-t28950477264453945225407686866-db0a584c` | GF(p^6), p 32 bits | 164 | 97.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +16.9 | +20.4 | formula_advantage |
| `icv1-fp36k6-tm321040139606615739387616393541050-3d7c8ad5` | GF(p^6), p 36 bits | 180 | 105.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +16.8 | +20.4 | formula_advantage |
| `icv1-fp36k6-tm434810334729911370096206102964806-35f04903` | GF(p^6), p 36 bits | 180 | 105.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +16.8 | +20.3 | formula_advantage |
| `icv1-fp36k6-tm434810334729911370096206102964806-5a5600d8` | GF(p^6), p 36 bits | 180 | 105.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +16.8 | +20.3 | formula_advantage |
| `icv1-fp40k6-tm159434457918589210674937944722440790-3eed8b93` | GF(p^6), p 40 bits | 196 | 113.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +16.7 | +20.3 | formula_advantage |
| c2pnb176w1 (`icv1-f2m176-t38161931374797611715595947-5639a6d4`) | GF(2^176) | 160 | 95.0 | GHS descent to Jac(C)/GF(2^2), N = 88 * 2 | 128 | 74.8 | +21.9 | +20.1 | formula_advantage |
| `icv1-fp38k6-tm22590552414788011593839745152427434-d4fae634` | GF(p^6), p 38 bits | 187 | 109.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +16.5 | +20.1 | formula_advantage |
| `icv1-fp34k6-t9155156567012472540378130718678-f978945a` | GF(p^6), p 34 bits | 171 | 100.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.9 | +16.2 | +19.7 | formula_advantage |
| `icv1-fp40k6-tm1863913783112918683701033322387770322-50685d0a` | GF(p^6), p 40 bits | 194 | 112.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +16.1 | +19.7 | formula_advantage |
| `icv1-fp42k6-tm149883592739821990689174352259577762446-37928db7` | GF(p^6), p 42 bits | 202 | 116.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.3 | +15.9 | +19.6 | formula_advantage |
| `icv1-fp32k6-t83035224272727947564063586618-0af66eb6` | GF(p^6), p 32 bits | 162 | 96.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +16.0 | +19.5 | formula_advantage |
| `icv1-fp38k6-tm17150065526389465328380381471097614-6ec15154` | GF(p^6), p 38 bits | 186 | 108.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +15.9 | +19.5 | formula_advantage |
| `icv1-fp36k6-tm518602515804779808627202906021202-7fb07e79` | GF(p^6), p 36 bits | 178 | 104.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +15.9 | +19.5 | formula_advantage |
| `icv1-fp32k6-tm26895410567024103251795854646-66d9cbf6` | GF(p^6), p 32 bits | 162 | 96.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +15.9 | +19.4 | formula_advantage |
| `icv1-fp32k6-tm26895410567024103251795854646-b23bda52` | GF(p^6), p 32 bits | 162 | 96.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +15.9 | +19.4 | formula_advantage |
| `icv1-fp42k6-t113956595982069419360252349593637121274-a07293a9` | GF(p^6), p 42 bits | 202 | 116.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.3 | +15.7 | +19.4 | formula_advantage |
| `icv1-fp42k6-t113956595982069419360252349593637121274-a119c665` | GF(p^6), p 42 bits | 202 | 116.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.3 | +15.7 | +19.4 | formula_advantage |
| `icv1-fp40k6-t783481071997334004792008993315150290-d5fe4ea4` | GF(p^6), p 40 bits | 194 | 112.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +15.7 | +19.3 | formula_advantage |
| `icv1-fp32k6-tm56939309847052024180839379882-cde6c446` | GF(p^6), p 32 bits | 162 | 96.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +15.7 | +19.2 | formula_advantage |
| `icv1-fp36k6-tm551406749111006501290465557387398-9cc21384` | GF(p^6), p 36 bits | 177 | 104.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +15.5 | +19.1 | formula_advantage |
| `icv1-fp36k6-tm551406749111006501290465557387398-ae668a76` | GF(p^6), p 36 bits | 177 | 104.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +15.5 | +19.1 | formula_advantage |
| `icv1-fp38k6-tm15617623173157607185862112772323686-e5340038` | GF(p^6), p 38 bits | 185 | 108.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +15.3 | +18.9 | formula_advantage |
| `icv1-fp40k6-tm2184879796728907784779547342147671030-debbb703` | GF(p^6), p 40 bits | 193 | 111.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +15.1 | +18.7 | formula_advantage |
| `icv1-fp40k6-tm1397487368265228690211894845558797786-738eb37b` | GF(p^6), p 40 bits | 192 | 111.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +15.1 | +18.7 | formula_advantage |
| `icv1-fp40k6-tm409168471856020180342889123576194214-554a2975` | GF(p^6), p 40 bits | 192 | 111.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +14.9 | +18.6 | formula_advantage |
| `icv1-fp38k6-tm5004407166360428023714402207672450-cdf9b534` | GF(p^6), p 38 bits | 184 | 107.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +15.0 | +18.6 | formula_advantage |
| `icv1-fp32k6-t93288269984182716629364799690-3cad6dc2` | GF(p^6), p 32 bits | 160 | 95.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +14.9 | +18.5 | formula_advantage |
| `icv1-fp32k6-t93288269984182716629364799690-abbc4349` | GF(p^6), p 32 bits | 160 | 95.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +14.9 | +18.5 | formula_advantage |
| `icv1-fp40k6-t481560125841730908343892522527130474-63c78358` | GF(p^6), p 40 bits | 192 | 111.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +14.8 | +18.4 | formula_advantage |
| `icv1-fp42k6-tm158547333510555401346139301564490547882-b8428891` | GF(p^6), p 42 bits | 200 | 115.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.3 | +14.8 | +18.4 | formula_advantage |
| `icv1-fp42k6-t55486372619441781840536659695190091330-f5de4252` | GF(p^6), p 42 bits | 200 | 115.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.3 | +14.6 | +18.3 | formula_advantage |
| `icv1-fp38k6-tm39647719120351838730753721191063210-d0bba3ff` | GF(p^6), p 38 bits | 184 | 107.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +14.6 | +18.2 | formula_advantage |
| `icv1-fp32k6-tm79862672573156897455635774390-97e2d874` | GF(p^6), p 32 bits | 160 | 95.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +14.7 | +18.2 | formula_advantage |
| `icv1-fp34k6-tm5399082817581527030251892957314-6f528201` | GF(p^6), p 34 bits | 168 | 99.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.9 | +14.6 | +18.2 | formula_advantage |
| `icv1-fp38k6-tm35083093762090154370752461436721094-9e2f77d4` | GF(p^6), p 38 bits | 184 | 107.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +14.6 | +18.2 | formula_advantage |
| `icv1-fp38k6-tm35083093762090154370752461436721094-d0ae921d` | GF(p^6), p 38 bits | 184 | 107.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +14.6 | +18.2 | formula_advantage |
| `icv1-fp34k6-t6697515153415246422853117310550-14c3c823` | GF(p^6), p 34 bits | 167 | 98.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.9 | +14.4 | +18.0 | formula_advantage |
| `icv1-fp38k6-t9792496143475352021826444200954378-0eee15f6` | GF(p^6), p 38 bits | 183 | 107.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +14.3 | +17.9 | formula_advantage |
| `icv1-fp38k6-t9792496143475352021826444200954378-67114235` | GF(p^6), p 38 bits | 183 | 107.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +14.3 | +17.9 | formula_advantage |
| `icv1-fp36k6-tm248387321097047254046941931850846-cde3ed22` | GF(p^6), p 36 bits | 175 | 102.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +14.3 | +17.9 | formula_advantage |
| `icv1-fp36k6-tm490498038359306366180694999614966-ea313251` | GF(p^6), p 36 bits | 175 | 102.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +14.3 | +17.9 | formula_advantage |
| `icv1-fp36k6-tm490498038359306366180694999614966-ea655c18` | GF(p^6), p 36 bits | 175 | 102.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +14.3 | +17.9 | formula_advantage |
| `icv1-fp38k6-t17707299750593443726436498443849798-f3e50591` | GF(p^6), p 38 bits | 183 | 106.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +14.1 | +17.8 | formula_advantage |
| `icv1-fp34k6-tm3832143756343999437227421675406-4d140473` | GF(p^6), p 34 bits | 167 | 98.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.9 | +14.2 | +17.8 | formula_advantage |
| `icv1-fp34k6-tm3832143756343999437227421675406-62690bce` | GF(p^6), p 34 bits | 167 | 98.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.9 | +14.2 | +17.8 | formula_advantage |
| `icv1-fp34k6-tm3832143756343999437227421675406-9e31e5f6` | GF(p^6), p 34 bits | 167 | 98.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.9 | +14.2 | +17.8 | formula_advantage |
| `icv1-fp32k6-t1263056764387787201622557442-b76b5f16` | GF(p^6), p 32 bits | 159 | 94.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +14.2 | +17.8 | formula_advantage |
| `icv1-fp34k6-t18143208094393754906718265710-235aa288` | GF(p^6), p 34 bits | 167 | 98.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.9 | +14.1 | +17.7 | formula_advantage |
| `icv1-fp36k6-tm404978269892502675061882653125086-5069e170` | GF(p^6), p 36 bits | 174 | 102.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +14.1 | +17.7 | formula_advantage |
| `icv1-fp32k6-tm30810004660830079984930045046-1141eaaf` | GF(p^6), p 32 bits | 158 | 94.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +14.0 | +17.6 | formula_advantage |
| `icv1-fp32k6-tm30810004660830079984930045046-ddbdcce3` | GF(p^6), p 32 bits | 158 | 94.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +14.0 | +17.6 | formula_advantage |
| `icv1-fp38k6-t22590552414788011593839745152427434-2212b80b` | GF(p^6), p 38 bits | 182 | 106.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +13.9 | +17.6 | formula_advantage |
| `icv1-fp38k6-t22590552414788011593839745152427434-39fbf9c2` | GF(p^6), p 38 bits | 182 | 106.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.1 | +13.9 | +17.6 | formula_advantage |
| `icv1-fp42k6-t3757559998291327151324708095133477130-26568e34` | GF(p^6), p 42 bits | 198 | 114.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.3 | +13.7 | +17.4 | formula_advantage |
| `icv1-fp42k6-t3757559998291327151324708095133477130-b97d288c` | GF(p^6), p 42 bits | 198 | 114.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.3 | +13.7 | +17.4 | formula_advantage |
| `icv1-fp18k14-tm68395056172132912118537176034894976142-90e1a450` | GF(p^14), p 18 bits | 247 | 139.5 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +15.5 | +17.4 | formula_advantage |
| `icv1-fp32k6-tm32816149465264987946732403194-825e0267` | GF(p^6), p 32 bits | 158 | 94.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.8 | +13.7 | +17.3 | formula_advantage |
| `icv1-fp40k6-tm1990693073699655257162695508951311630-71d70d5b` | GF(p^6), p 40 bits | 189 | 110.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +13.6 | +17.2 | formula_advantage |
| `icv1-fp40k6-tm1990693073699655257162695508951311630-9e700f84` | GF(p^6), p 40 bits | 189 | 110.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +13.6 | +17.2 | formula_advantage |
| `icv1-fp42k6-tm85965374489496676757941558323912879818-b0b7b0e2` | GF(p^6), p 42 bits | 197 | 114.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.3 | +13.4 | +17.1 | formula_advantage |
| `icv1-fp36k6-tm260583852068831826778382984658130-7402b3bc` | GF(p^6), p 36 bits | 173 | 102.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 85.0 | +13.4 | +17.0 | formula_advantage |
| `icv1-fp40k6-t1327766184040927902463371280455357170-6ed928c6` | GF(p^6), p 40 bits | 189 | 110.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.2 | +13.3 | +17.0 | formula_advantage |
| `icv1-fp34k6-t1060509355543416520796049580354-266d8fdd` | GF(p^6), p 34 bits | 165 | 97.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +13.4 | +17.0 | formula_advantage |
| `icv1-fp34k6-t1060509355543416520796049580354-ae125061` | GF(p^6), p 34 bits | 165 | 97.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +13.4 | +17.0 | formula_advantage |
| `icv1-fp32k6-tm39645023259520707621979274842-d38d2667` | GF(p^6), p 32 bits | 157 | 93.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +13.4 | +17.0 | formula_advantage |
| `icv1-fp32k6-t10780078424225041273621456910-c93319a0` | GF(p^6), p 32 bits | 157 | 93.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +13.4 | +17.0 | formula_advantage |
| `icv1-fp34k6-t8343469462759645670504383165214-af0df25f` | GF(p^6), p 34 bits | 165 | 97.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +13.1 | +16.8 | formula_advantage |
| `icv1-fp18k14-tm139738101568416365238979260318246451498-eda2c725` | GF(p^14), p 18 bits | 246 | 138.8 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +14.9 | +16.7 | formula_advantage |
| `icv1-fp43k6-t80369948951667839190628429049567884210-1a345239` | GF(p^6), p 42 bits | 196 | 113.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +12.8 | +16.5 | formula_advantage |
| `icv1-fp43k6-t80369948951667839190628429049567884210-8a669356` | GF(p^6), p 42 bits | 196 | 113.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +12.8 | +16.5 | formula_advantage |
| `icv1-fp32k6-t10656495889732328400009566070-a889c117` | GF(p^6), p 32 bits | 156 | 93.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +12.7 | +16.3 | formula_advantage |
| `icv1-fp42k6-t57664958956439933191375515704719752298-1c1c2637` | GF(p^6), p 42 bits | 195 | 113.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +12.4 | +16.1 | formula_advantage |
| `icv1-fp32k6-tm97949610295501091330605680258-e4a45239` | GF(p^6), p 32 bits | 155 | 92.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +12.4 | +16.1 | formula_advantage |
| `icv1-fp36k6-tm348163876081752902381128240571926-72201f11` | GF(p^6), p 36 bits | 171 | 101.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +12.4 | +16.0 | formula_advantage |
| `icv1-fp38k6-tm2717978647929456783301116716814934-0c867b8a` | GF(p^6), p 38 bits | 179 | 104.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +12.2 | +15.9 | formula_advantage |
| `icv1-fp38k6-tm2717978647929456783301116716814934-86c2c3ee` | GF(p^6), p 38 bits | 179 | 104.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +12.2 | +15.9 | formula_advantage |
| `icv1-fp38k6-tm24644369446264111877218863115945618-3973f0bb` | GF(p^6), p 38 bits | 178 | 104.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +11.9 | +15.7 | formula_advantage |
| `icv1-fp18k14-t23685949668102613212425493219205765758-1744fed4` | GF(p^14), p 18 bits | 244 | 137.8 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +13.8 | +15.6 | formula_advantage |
| `icv1-fp42k6-t164816127623881732294837773186127712426-2846d75c` | GF(p^6), p 42 bits | 194 | 112.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +11.9 | +15.6 | formula_advantage |
| `icv1-fp42k6-t164816127623881732294837773186127712426-8a07de52` | GF(p^6), p 42 bits | 194 | 112.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +11.9 | +15.6 | formula_advantage |
| `icv1-fp38k6-t32207907774093674990821893566327454-d2a3b684` | GF(p^6), p 38 bits | 178 | 104.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +11.8 | +15.5 | formula_advantage |
| `icv1-fp42k6-tm42685072077033061909020167905304827814-ed81cabb` | GF(p^6), p 42 bits | 194 | 112.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +11.8 | +15.5 | formula_advantage |
| `icv1-fp32k6-tm48482309647343732500014249610-25ae5ee5` | GF(p^6), p 32 bits | 154 | 92.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +11.9 | +15.5 | formula_advantage |
| `icv1-fp34k6-t974880199944949524815491736754-060daf26` | GF(p^6), p 34 bits | 162 | 96.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +11.8 | +15.5 | formula_advantage |
| `icv1-fp34k6-t974880199944949524815491736754-610c01ef` | GF(p^6), p 34 bits | 162 | 96.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +11.8 | +15.5 | formula_advantage |
| `icv1-fp34k6-t698614099488561323391098718782-7b2e2e1c` | GF(p^6), p 34 bits | 162 | 96.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +11.8 | +15.5 | formula_advantage |
| `icv1-fp34k6-t3847521778517058181815528842402-744ce002` | GF(p^6), p 34 bits | 162 | 96.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +11.8 | +15.4 | formula_advantage |
| `icv1-fp34k6-t3847521778517058181815528842402-b3283bbd` | GF(p^6), p 34 bits | 162 | 96.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +11.8 | +15.4 | formula_advantage |
| `icv1-fp40k6-tm494811548722146648595995475931754494-4912fe8e` | GF(p^6), p 40 bits | 186 | 108.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.1 | +11.6 | +15.4 | formula_advantage |
| `icv1-fp32k6-tm127812454930946669422317632054-0087f3f7` | GF(p^6), p 32 bits | 154 | 92.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +11.7 | +15.3 | formula_advantage |
| `icv1-fp32k6-tm127812454930946669422317632054-fdfff12c` | GF(p^6), p 32 bits | 154 | 92.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +11.7 | +15.3 | formula_advantage |
| `icv1-fp38k6-tm16701024292395226668393629918546758-80458c25` | GF(p^6), p 38 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +11.6 | +15.3 | formula_advantage |
| `icv1-fp38k6-tm16701024292395226668393629918546758-9d8f8dd3` | GF(p^6), p 38 bits | 178 | 104.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +11.6 | +15.3 | formula_advantage |
| `icv1-fp38k6-tm8927393401872064762098407276967838-f9f4315f` | GF(p^6), p 38 bits | 177 | 104.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +11.5 | +15.2 | formula_advantage |
| `icv1-fp36k6-tm131474894390697220914257827196726-3321f215` | GF(p^6), p 36 bits | 169 | 100.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +11.5 | +15.2 | formula_advantage |
| `icv1-fp36k6-tm131474894390697220914257827196726-6b43aeaa` | GF(p^6), p 36 bits | 169 | 100.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +11.5 | +15.2 | formula_advantage |
| `icv1-fp36k6-t574682551485435693406210031183306-5d2fc0f6` | GF(p^6), p 36 bits | 169 | 100.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +11.5 | +15.2 | formula_advantage |
| `icv1-fp18k14-t63438034467637535885309421389268171802-20b71da1` | GF(p^14), p 18 bits | 243 | 137.2 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +13.2 | +15.1 | formula_advantage |
| `icv1-fp18k14-t63438034467637535885309421389268171802-8544736e` | GF(p^14), p 18 bits | 243 | 137.2 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +13.2 | +15.1 | formula_advantage |
| `icv1-fp32k6-tm29242045006735922381596573882-c2f26c97` | GF(p^6), p 32 bits | 153 | 91.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +11.4 | +15.0 | formula_advantage |
| `icv1-fp42k6-tm18224734467133666853515094586548756450-be9d285c` | GF(p^6), p 42 bits | 193 | 112.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +11.2 | +14.9 | formula_advantage |
| `icv1-fp32k6-tm15813541453371577054675899158-3832dc05` | GF(p^6), p 32 bits | 153 | 91.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +11.2 | +14.9 | formula_advantage |
| `icv1-fp32k6-tm15813541453371577054675899158-96a479b8` | GF(p^6), p 32 bits | 153 | 91.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +11.2 | +14.9 | formula_advantage |
| `icv1-fp18k14-t33740855631107397075996375079172915058-03057729` | GF(p^14), p 18 bits | 242 | 136.8 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +12.9 | +14.7 | formula_advantage |
| `icv1-fp34k6-t2941176509346770792185269096662-c56a249b` | GF(p^6), p 34 bits | 160 | 95.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +10.8 | +14.5 | formula_advantage |
| `icv1-fp32k6-tm144849842750884344644479759526-38668902` | GF(p^6), p 32 bits | 152 | 91.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +10.8 | +14.5 | formula_advantage |
| `icv1-fp36k6-tm210075224010930853915881415035702-42b5d3c6` | GF(p^6), p 36 bits | 168 | 99.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +10.8 | +14.5 | formula_advantage |
| `icv1-fp36k6-tm210075224010930853915881415035702-5eeae033` | GF(p^6), p 36 bits | 168 | 99.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +10.8 | +14.5 | formula_advantage |
| `icv1-fp18k14-tm18873402740655092695313128842762072054-663f050a` | GF(p^14), p 18 bits | 241 | 136.5 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +12.5 | +14.4 | formula_advantage |
| `icv1-fp32k6-t16589822727447947892730390010-9a3fe949` | GF(p^6), p 32 bits | 152 | 91.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +10.7 | +14.4 | formula_advantage |
| `icv1-fp42k6-tm1449329879519965966669504876259048666-466f3e6e` | GF(p^6), p 42 bits | 191 | 111.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +10.5 | +14.3 | formula_advantage |
| `icv1-fp36k6-t97034101501832544519362719691618-02e835a5` | GF(p^6), p 36 bits | 167 | 99.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +10.5 | +14.3 | formula_advantage |
| `icv1-fp34k6-t289226834666250011557687440226-8487a9fe` | GF(p^6), p 34 bits | 159 | 95.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +10.5 | +14.2 | formula_advantage |
| `icv1-fp34k6-t289226834666250011557687440226-c22685a5` | GF(p^6), p 34 bits | 159 | 95.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +10.5 | +14.2 | formula_advantage |
| `icv1-fp34k6-tm589547474161273649282298788622-30d84d46` | GF(p^6), p 34 bits | 159 | 95.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +10.5 | +14.2 | formula_advantage |
| `icv1-fp36k6-tm349583752751787938470935204879638-0355401d` | GF(p^6), p 36 bits | 167 | 99.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +10.4 | +14.2 | formula_advantage |
| `icv1-fp36k6-tm349583752751787938470935204879638-3e318d2e` | GF(p^6), p 36 bits | 167 | 99.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +10.4 | +14.2 | formula_advantage |
| `icv1-fp32k6-tm87449312845250641537531312086-1e4f6767` | GF(p^6), p 32 bits | 151 | 90.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +10.4 | +14.1 | formula_advantage |
| `icv1-fp32k6-tm87449312845250641537531312086-7e928d53` | GF(p^6), p 32 bits | 151 | 90.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.7 | +10.4 | +14.1 | formula_advantage |
| `icv1-fp38k6-t21778159553746265407774699488272586-77a259b3` | GF(p^6), p 38 bits | 175 | 103.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +10.3 | +14.0 | formula_advantage |
| `icv1-fp32k6-t12160693092578909686110443866-de626226` | GF(p^6), p 32 bits | 151 | 90.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +10.2 | +13.9 | formula_advantage |
| `icv1-fp32k6-tm10656495889732328400009566070-093e24db` | GF(p^6), p 32 bits | 151 | 90.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +10.2 | +13.9 | formula_advantage |
| `icv1-fp40k6-tm717313670595015088776681899430335030-28aa6341` | GF(p^6), p 40 bits | 183 | 107.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.1 | +10.1 | +13.9 | formula_advantage |
| `icv1-fp42k6-tm29792387245275262226064940024003903270-4a90d2c2` | GF(p^6), p 42 bits | 191 | 111.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +10.1 | +13.9 | formula_advantage |
| `icv1-fp42k6-tm29792387245275262226064940024003903270-5788c993` | GF(p^6), p 42 bits | 191 | 111.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.2 | +10.1 | +13.9 | formula_advantage |
| `icv1-fp32k6-t46788840575897536961008224826-b87427d0` | GF(p^6), p 32 bits | 151 | 90.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +10.2 | +13.9 | formula_advantage |
| `icv1-fp34k6-tm8187850875526795465613921623090-3c79c0ad` | GF(p^6), p 34 bits | 159 | 94.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.8 | +10.1 | +13.8 | formula_advantage |
| `icv1-fp38k6-tm5849615460972165484256291471884378-ac989a5f` | GF(p^6), p 38 bits | 174 | 102.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 89.0 | +9.9 | +13.7 | formula_advantage |
| `icv1-fp36k6-tm466054673267443778344436740368950-ae6f9d4d` | GF(p^6), p 36 bits | 165 | 98.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.9 | +9.6 | +13.3 | formula_advantage |
| `icv1-fp34k6-tm3847521778517058181815528842402-77b526ad` | GF(p^6), p 34 bits | 157 | 94.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +9.5 | +13.3 | formula_advantage |
| `icv1-fp34k6-t1654236123937356050030227467706-305d14f0` | GF(p^6), p 34 bits | 157 | 94.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +9.5 | +13.2 | formula_advantage |
| `icv1-fp18k14-t87408846862324945293719340911362845178-419a5a75` | GF(p^14), p 18 bits | 239 | 135.2 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +11.2 | +13.1 | formula_advantage |
| `icv1-fp18k14-t87408846862324945293719340911362845178-77e9b73e` | GF(p^14), p 18 bits | 239 | 135.2 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +11.2 | +13.1 | formula_advantage |
| `icv1-fp18k14-tm24692850434072448451442441489471858942-c6324482` | GF(p^14), p 18 bits | 239 | 135.2 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +11.2 | +13.1 | formula_advantage |
| `icv1-fp32k6-tm33815349211666425312826761746-7796b3c4` | GF(p^6), p 32 bits | 149 | 89.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +9.2 | +12.9 | formula_advantage |
| `icv1-fp18k14-t40901129643050925884153966080961450506-e2c0acb3` | GF(p^14), p 18 bits | 238 | 135.0 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +11.0 | +12.9 | formula_advantage |
| `icv1-fp18k14-t129136710913319553601414566440572657082-05da67b5` | GF(p^14), p 18 bits | 238 | 134.8 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +10.8 | +12.7 | formula_advantage |
| `icv1-fp18k14-t129136710913319553601414566440572657082-fcfc140c` | GF(p^14), p 18 bits | 238 | 134.8 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +10.8 | +12.7 | formula_advantage |
| `icv1-fp34k6-tm222467787949046006723433574990-1b6feaeb` | GF(p^6), p 34 bits | 156 | 93.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +8.9 | +12.7 | formula_advantage |
| `icv1-fp34k6-t783439549876049572437540640018-5202f3ba` | GF(p^6), p 34 bits | 156 | 93.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +8.8 | +12.6 | formula_advantage |
| `icv1-fp34k6-t783439549876049572437540640018-bbf921c2` | GF(p^6), p 34 bits | 156 | 93.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +8.8 | +12.6 | formula_advantage |
| `icv1-fp42k6-t71034451563579483295596031468022725890-9bfd94a7` | GF(p^6), p 42 bits | 188 | 109.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +8.6 | +12.4 | formula_advantage |
| `icv1-fp36k6-t321040139606615739387616393541050-85445ca6` | GF(p^6), p 36 bits | 163 | 97.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.8 | +8.6 | +12.3 | formula_advantage |
| `icv1-fp36k6-t321040139606615739387616393541050-bb6f3ed2` | GF(p^6), p 36 bits | 163 | 97.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.8 | +8.6 | +12.3 | formula_advantage |
| `icv1-fp42k6-t33767737116018284517469893751560331242-6dd5a751` | GF(p^6), p 42 bits | 187 | 109.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +8.5 | +12.3 | formula_advantage |
| `icv1-fp42k6-t33767737116018284517469893751560331242-70c1a98d` | GF(p^6), p 42 bits | 187 | 109.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +8.5 | +12.3 | formula_advantage |
| `icv1-fp42k6-t33767737116018284517469893751560331242-9aa7cea5` | GF(p^6), p 42 bits | 187 | 109.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +8.5 | +12.3 | formula_advantage |
| `icv1-fp38k6-tm31116608428445849544095091174863146-9758fae0` | GF(p^6), p 38 bits | 171 | 101.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.9 | +8.4 | +12.2 | formula_advantage |
| `icv1-fp36k6-tm417243273485249998640231485044782-d326535c` | GF(p^6), p 36 bits | 163 | 97.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.8 | +8.4 | +12.2 | formula_advantage |
| `icv1-fp18k14-tm154128700048226856270319739242039973614-69549480` | GF(p^14), p 18 bits | 237 | 134.3 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +10.3 | +12.2 | formula_advantage |
| `icv1-fp38k6-t37724595884107578164257832721671138-907d1ac2` | GF(p^6), p 38 bits | 171 | 100.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.9 | +8.2 | +12.0 | formula_advantage |
| `icv1-fp32k6-tm132331583280069005213348294498-ea4339de` | GF(p^6), p 32 bits | 147 | 88.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +8.1 | +11.9 | formula_advantage |
| `icv1-fp18k14-t95683512112052388831930406714508243870-82e54984` | GF(p^14), p 18 bits | 236 | 133.8 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +9.9 | +11.8 | formula_advantage |
| `icv1-fp32k6-t74857419110022024462382264786-d439417e` | GF(p^6), p 32 bits | 146 | 88.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +7.9 | +11.6 | formula_advantage |
| `icv1-fp32k6-t48482309647343732500014249610-0ee38ac6` | GF(p^6), p 32 bits | 146 | 88.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +7.9 | +11.6 | formula_advantage |
| `icv1-fp32k6-tm50709883291525901472845329138-74ceca7e` | GF(p^6), p 32 bits | 146 | 88.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +7.7 | +11.5 | formula_advantage |
| `icv1-fp34k6-t2949540786897853555642082792130-0d1502a4` | GF(p^6), p 34 bits | 154 | 92.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +7.7 | +11.5 | formula_advantage |
| `icv1-fp34k6-t2949540786897853555642082792130-6c7be75e` | GF(p^6), p 34 bits | 154 | 92.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +7.7 | +11.5 | formula_advantage |
| `icv1-fp34k6-t3415120380487737561606789344386-3e4cca0d` | GF(p^6), p 34 bits | 153 | 91.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +7.4 | +11.2 | formula_advantage |
| `icv1-fp34k6-t3415120380487737561606789344386-6cadb70a` | GF(p^6), p 34 bits | 153 | 91.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +7.4 | +11.2 | formula_advantage |
| `icv1-fp34k6-tm1096768777423997627962161807350-536bafb9` | GF(p^6), p 34 bits | 153 | 91.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +7.3 | +11.1 | formula_advantage |
| `icv1-fp38k6-tm12628575110007103759117203599212102-05bbf693` | GF(p^6), p 38 bits | 169 | 99.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.9 | +7.2 | +11.1 | formula_advantage |
| `icv1-fp38k6-tm12628575110007103759117203599212102-162092c1` | GF(p^6), p 38 bits | 169 | 99.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.9 | +7.2 | +11.1 | formula_advantage |
| `icv1-fp34k6-t2371564097477158673272135904890-984422a9` | GF(p^6), p 34 bits | 153 | 91.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +7.2 | +11.0 | formula_advantage |
| `icv1-fp42k6-tm40975292624622716692299527262814166850-95be9ad7` | GF(p^6), p 42 bits | 185 | 108.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +7.1 | +11.0 | formula_advantage |
| `icv1-fp32k6-t59621373577963920613733060330-0b23dffa` | GF(p^6), p 32 bits | 145 | 87.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +7.1 | +10.9 | formula_advantage |
| `icv1-fp32k6-t59621373577963920613733060330-65231570` | GF(p^6), p 32 bits | 145 | 87.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.6 | +7.1 | +10.9 | formula_advantage |
| `icv1-fp38k6-t5849615460972165484256291471884378-09913842` | GF(p^6), p 38 bits | 168 | 99.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.9 | +7.1 | +10.9 | formula_advantage |
| `icv1-fp38k6-t5849615460972165484256291471884378-602e1868` | GF(p^6), p 38 bits | 168 | 99.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.9 | +7.1 | +10.9 | formula_advantage |
| `icv1-fp42k6-t139397368184640235945825692291306341094-55d568f8` | GF(p^6), p 42 bits | 184 | 108.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +7.0 | +10.9 | formula_advantage |
| `icv1-fp40k6-t1138014948113436022385175216268138942-cc40d3c0` | GF(p^6), p 40 bits | 176 | 103.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.0 | +7.0 | +10.8 | formula_advantage |
| `icv1-fp42k6-t69027207743865244667769899216784826010-263c4311` | GF(p^6), p 42 bits | 184 | 107.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +6.8 | +10.7 | formula_advantage |
| `icv1-fp42k6-t69027207743865244667769899216784826010-961830b9` | GF(p^6), p 42 bits | 184 | 107.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +6.8 | +10.7 | formula_advantage |
| `icv1-fp18k14-t139738101568416365238979260318246451498-b39a6c89` | GF(p^14), p 18 bits | 234 | 132.7 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +8.7 | +10.6 | formula_advantage |
| `icv1-fp18k14-t139738101568416365238979260318246451498-c75c943d` | GF(p^14), p 18 bits | 234 | 132.7 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +8.7 | +10.6 | formula_advantage |
| `icv1-fp34k6-t746969533550109445373866422690-f094927f` | GF(p^6), p 34 bits | 152 | 91.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.7 | +6.7 | +10.5 | formula_advantage |
| `icv1-fp40k6-t1469511672384069112375443404508791802-80262b90` | GF(p^6), p 40 bits | 176 | 103.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.0 | +6.6 | +10.5 | formula_advantage |
| `icv1-fp38k6-tm17707299750593443726436498443849798-9061bf8a` | GF(p^6), p 38 bits | 167 | 99.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.9 | +6.4 | +10.3 | formula_advantage |
| `icv1-fp38k6-tm17707299750593443726436498443849798-ad68d576` | GF(p^6), p 38 bits | 167 | 99.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.9 | +6.4 | +10.3 | formula_advantage |
| `icv1-fp36k6-t373548937192775747729002977951618-c3f1cf60` | GF(p^6), p 36 bits | 158 | 94.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.8 | +6.0 | +9.8 | formula_advantage |
| `icv1-fp40k6-t13709315706324085134265846307139378-138a3148` | GF(p^6), p 40 bits | 174 | 102.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.0 | +5.8 | +9.7 | formula_advantage |
| `icv1-fp40k6-t13709315706324085134265846307139378-8d52cf8d` | GF(p^6), p 40 bits | 174 | 102.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.0 | +5.8 | +9.7 | formula_advantage |
| `icv1-fp32k6-t121298768510003002706056856930-0eb56463` | GF(p^6), p 32 bits | 142 | 86.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.5 | +5.8 | +9.6 | formula_advantage |
| `icv1-fp36k6-tm379043372837942653810799561381086-bae1552e` | GF(p^6), p 36 bits | 158 | 94.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.7 | +5.7 | +9.5 | formula_advantage |
| `icv1-fp42k6-tm149433351988728729398325948808771718870-06af03af` | GF(p^6), p 42 bits | 182 | 106.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +5.6 | +9.5 | formula_advantage |
| `icv1-fp42k6-tm149433351988728729398325948808771718870-ef8bdf5a` | GF(p^6), p 42 bits | 182 | 106.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +5.6 | +9.5 | formula_advantage |
| `icv1-fp32k6-t95346967583027824387079081114-2b20be20` | GF(p^6), p 32 bits | 142 | 85.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.5 | +5.6 | +9.4 | formula_advantage |
| `icv1-fp32k6-t95346967583027824387079081114-df0d818a` | GF(p^6), p 32 bits | 142 | 85.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.5 | +5.6 | +9.4 | formula_advantage |
| `icv1-fp40k6-tm1539986654890541381168714181014265694-b3c77886` | GF(p^6), p 40 bits | 173 | 102.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 93.0 | +5.2 | +9.1 | formula_advantage |
| `icv1-fp42k6-tm158444610310601711509643723342538028630-4fe3c2bb` | GF(p^6), p 42 bits | 181 | 106.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +5.1 | +9.1 | formula_advantage |
| `icv1-fp42k6-tm158444610310601711509643723342538028630-f6366d3e` | GF(p^6), p 42 bits | 181 | 106.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 97.1 | +5.1 | +9.1 | formula_advantage |
| `icv1-fp38k6-t15589384784277950280908897000583130-9230dadf` | GF(p^6), p 38 bits | 165 | 97.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.8 | +5.1 | +9.0 | formula_advantage |
| `icv1-fp38k6-t15589384784277950280908897000583130-b9063068` | GF(p^6), p 38 bits | 165 | 97.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.8 | +5.1 | +9.0 | formula_advantage |
| `icv1-fp36k6-tm430636157326068713228412158502122-04f165ac` | GF(p^6), p 36 bits | 156 | 93.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.7 | +5.1 | +9.0 | formula_advantage |
| `icv1-fp32k6-tm109493713542271466535275947182-6175dc54` | GF(p^6), p 32 bits | 140 | 85.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.5 | +5.1 | +8.9 | formula_advantage |
| `icv1-fp38k6-t2706859881483543508812126344841126-c357cc9e` | GF(p^6), p 38 bits | 164 | 97.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.8 | +4.9 | +8.8 | formula_advantage |
| `icv1-fp32k6-tm46788840575897536961008224826-04c6b389` | GF(p^6), p 32 bits | 140 | 85.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.5 | +4.8 | +8.6 | formula_advantage |
| `icv1-fp40k6-t1195900778100877146770415087895037282-7e17ef8f` | GF(p^6), p 40 bits | 172 | 101.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 92.9 | +4.7 | +8.6 | formula_advantage |
| `icv1-fp40k6-t1195900778100877146770415087895037282-e97f74a2` | GF(p^6), p 40 bits | 172 | 101.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 92.9 | +4.7 | +8.6 | formula_advantage |
| `icv1-fp40k6-tm1233373951829317075315162348727047814-58baec14` | GF(p^6), p 40 bits | 171 | 101.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 92.9 | +4.6 | +8.5 | formula_advantage |
| `icv1-fp18k14-t7256866063702034327501088129872235306-883ad044` | GF(p^14), p 18 bits | 229 | 130.4 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +6.4 | +8.3 | formula_advantage |
| `icv1-fp18k14-t7256866063702034327501088129872235306-b32c10b8` | GF(p^14), p 18 bits | 229 | 130.4 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +6.4 | +8.3 | formula_advantage |
| `icv1-fp18k14-t4408213846726974963031874993567538022-d8451fd3` | GF(p^14), p 18 bits | 229 | 130.3 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +6.3 | +8.2 | formula_advantage |
| `icv1-fp32k6-t59292859294895793137471563082-2c62431a` | GF(p^6), p 32 bits | 139 | 84.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.5 | +4.2 | +8.0 | formula_advantage |
| `icv1-fp32k6-t59292859294895793137471563082-d83726fb` | GF(p^6), p 32 bits | 139 | 84.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.5 | +4.2 | +8.0 | formula_advantage |
| `icv1-fp18k14-t78853452396314113567507418714538409866-08e16989` | GF(p^14), p 18 bits | 229 | 130.1 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +6.1 | +8.0 | formula_advantage |
| `icv1-fp18k14-t78853452396314113567507418714538409866-633b0cff` | GF(p^14), p 18 bits | 229 | 130.1 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +6.1 | +8.0 | formula_advantage |
| `icv1-fp34k6-t6409468687277400196484097419506-16550513` | GF(p^6), p 34 bits | 146 | 88.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.6 | +4.0 | +7.9 | formula_marginal |
| `icv1-fp36k6-t502103286078933538421439499726062-e6318d8d` | GF(p^6), p 36 bits | 154 | 92.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.7 | +3.9 | +7.8 | formula_marginal |
| `icv1-fp34k6-t5502116476237580591784371210930-f3230475` | GF(p^6), p 34 bits | 146 | 88.3 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.6 | +3.8 | +7.7 | formula_marginal |
| `icv1-fp18k14-tm41642808330359025650899351199891334938-79609368` | GF(p^14), p 18 bits | 228 | 129.6 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +5.7 | +7.6 | formula_marginal |
| `icv1-fp34k6-t3773347342789427156506271585042-dfe37e37` | GF(p^6), p 34 bits | 146 | 88.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.6 | +3.6 | +7.5 | formula_marginal |
| `icv1-fp38k6-t25703862405918691689370241743508566-7bcea2e9` | GF(p^6), p 38 bits | 161 | 96.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.8 | +3.3 | +7.3 | formula_marginal |
| `icv1-fp32k6-tm79952562419435043747240341346-93f9aeee` | GF(p^6), p 32 bits | 137 | 83.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.4 | +3.3 | +7.2 | formula_marginal |
| `icv1-fp36k6-tm266347467794506259203439006117222-335574d8` | GF(p^6), p 36 bits | 153 | 91.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.7 | +3.2 | +7.1 | formula_marginal |
| `icv1-fp36k6-tm266347467794506259203439006117222-4bdceb53` | GF(p^6), p 36 bits | 153 | 91.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.7 | +3.2 | +7.1 | formula_marginal |
| `icv1-fp38k6-tm15589384784277950280908897000583130-2f4ca2b8` | GF(p^6), p 38 bits | 161 | 95.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.8 | +3.1 | +7.0 | formula_marginal |
| `icv1-fp36k6-t294614919642168700237241526426890-16738747` | GF(p^6), p 36 bits | 152 | 91.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.7 | +3.0 | +6.9 | formula_marginal |
| `icv1-fp36k6-t294614919642168700237241526426890-996fd5c4` | GF(p^6), p 36 bits | 152 | 91.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.7 | +3.0 | +6.9 | formula_marginal |
| `icv1-fp36k6-t171229324292869415575893259583430-72ee8f04` | GF(p^6), p 36 bits | 152 | 91.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.7 | +3.0 | +6.9 | formula_marginal |
| `icv1-fp21k6-tm927079723928371550-0fec2c4a` | GF(p^6), p 20 bits | 98 | 63.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 56.2 | +7.7 | +6.9 | formula_marginal |
| `icv1-fp21k6-tm927079723928371550-a499fc7b` | GF(p^6), p 20 bits | 98 | 63.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 56.2 | +7.7 | +6.9 | formula_marginal |
| `icv1-fp38k6-tm12905434611021849892188445314456482-fa250cd4` | GF(p^6), p 38 bits | 160 | 95.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.8 | +2.8 | +6.8 | formula_marginal |
| `icv1-fp32k6-tm16233462481222275852176981534-d7bfadc5` | GF(p^6), p 32 bits | 136 | 83.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.4 | +2.7 | +6.6 | formula_marginal |
| `icv1-fp18k14-tm4408213846726974963031874993567538022-bab278c1` | GF(p^14), p 18 bits | 225 | 128.5 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +4.5 | +6.4 | formula_marginal |
| `icv1-fp18k14-tm4408213846726974963031874993567538022-f1c11b32` | GF(p^14), p 18 bits | 225 | 128.5 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.1 | +4.5 | +6.4 | formula_marginal |
| `icv1-fp32k6-tm131394180832725690985430230374-05728ea2` | GF(p^6), p 32 bits | 135 | 82.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.4 | +2.3 | +6.3 | formula_marginal |
| `icv1-fp32k6-tm131394180832725690985430230374-de76c9e2` | GF(p^6), p 32 bits | 135 | 82.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.4 | +2.3 | +6.3 | formula_marginal |
| `icv1-fp40k6-tm1408510537005952335530741596990585198-0db40866` | GF(p^6), p 40 bits | 167 | 99.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 92.9 | +2.2 | +6.2 | formula_marginal |
| `icv1-fp40k6-tm272141120984694834845539129520103906-cea63929` | GF(p^6), p 40 bits | 166 | 98.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 92.9 | +2.0 | +6.0 | formula_marginal |
| `icv1-fp32k6-t28081572847806456169299407758-1b4561f1` | GF(p^6), p 32 bits | 133 | 81.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.4 | +1.5 | +5.4 | formula_marginal |
| `icv1-fp34k6-t2894955473854795684072614329918-4f93ca4c` | GF(p^6), p 34 bits | 141 | 85.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.5 | +1.4 | +5.4 | formula_marginal |
| `icv1-fp38k6-tm9075672861665887541181946242606758-29a0f1e6` | GF(p^6), p 38 bits | 157 | 93.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.7 | +1.2 | +5.2 | formula_marginal |
| `icv1-fp38k6-tm9075672861665887541181946242606758-b01197b8` | GF(p^6), p 38 bits | 157 | 93.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.7 | +1.2 | +5.2 | formula_marginal |
| `icv1-fp38k6-tm33653313885195831439390070846015878-6b996fe4` | GF(p^6), p 38 bits | 157 | 93.9 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.7 | +1.2 | +5.2 | formula_marginal |
| `icv1-fp32k6-t124926680297413374833471646370-be81d726` | GF(p^6), p 32 bits | 133 | 81.5 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.4 | +1.2 | +5.1 | formula_marginal |
| `icv1-fp40k6-tm1746682707660615888470782891457259390-28d0b224` | GF(p^6), p 40 bits | 164 | 97.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 92.8 | +0.8 | +4.8 | formula_marginal |
| `icv1-fp40k6-tm1746682707660615888470782891457259390-5c15d384` | GF(p^6), p 40 bits | 164 | 97.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 92.8 | +0.8 | +4.8 | formula_marginal |
| `icv1-fp32k6-t221648548951264495700929522-89d76651` | GF(p^6), p 32 bits | 132 | 81.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.3 | +0.7 | +4.7 | formula_marginal |
| `icv1-fp34k6-t997252301807352666642538383922-24168704` | GF(p^6), p 34 bits | 140 | 85.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.5 | +0.7 | +4.7 | formula_marginal |
| c2pnb368w1 (`icv1-f2m368-tm15466245949934416937841361495935610942842793802409489935-43ee62b3`) | GF(2^368) | 352 | 192.7 | isogenous curve with magic number <= 6 over the same tower, descent to Jac(C)/GF(2^46) | 32 | 188.1 | +13.3 | +4.6 | formula_marginal |
| `icv1-fp32k6-tm70499135745317542528367878850-e6e52034` | GF(p^6), p 32 bits | 131 | 80.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.3 | +0.5 | +4.5 | formula_marginal |
| Oakley Group 4 (`icv1-f2m185-t11157211627747266908830216517-8dbe23c4`) | GF(2^185) | 183 | 106.6 | GHS descent to Jac(C)/GF(2^37), N = 5 * 37 | 16 | 102.4 | +11.2 | +4.2 | formula_marginal |
| `icv1-fp32k6-t14484914975528427710430752602-7aec23dc` | GF(p^6), p 32 bits | 131 | 80.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.3 | +0.1 | +4.1 | formula_marginal |
| `icv1-fp32k6-t46036983227369337013558516990-9111b659` | GF(p^6), p 32 bits | 130 | 80.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.3 | +0.1 | +4.1 | formula_marginal |
| `icv1-fp42k6-t44110276869417112049281439078164335142-93d6ed6a` | GF(p^6), p 42 bits | 170 | 100.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 96.9 | -0.4 | +3.7 | formula_marginal |
| `icv1-fp36k6-t430833983857612770463572771532730-c3743b1d` | GF(p^6), p 36 bits | 145 | 88.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 84.6 | -0.6 | +3.5 | formula_marginal |
| `icv1-fp38k6-tm2706859881483543508812126344841126-59e222cc` | GF(p^6), p 38 bits | 153 | 92.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.7 | -0.7 | +3.4 | formula_marginal |
| `icv1-fp38k6-tm2706859881483543508812126344841126-7444197f` | GF(p^6), p 38 bits | 153 | 92.1 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.7 | -0.7 | +3.4 | formula_marginal |
| `icv1-fp18k14-tm106874586433925198625400569232147924590-698bd52c` | GF(p^14), p 18 bits | 219 | 125.3 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.0 | +1.3 | +3.3 | formula_marginal |
| `icv1-fp18k14-t105347263935672296906898885769304887882-8ebd9b57` | GF(p^14), p 18 bits | 218 | 124.9 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.0 | +1.0 | +2.9 | formula_marginal |
| `icv1-fp18k14-t105347263935672296906898885769304887882-b3440967` | GF(p^14), p 18 bits | 218 | 124.9 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.0 | +1.0 | +2.9 | formula_marginal |
| `icv1-fp32k6-tm34954230878988953321201855330-fbd976c1` | GF(p^6), p 32 bits | 127 | 78.8 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.3 | -1.5 | +2.6 | formula_marginal |
| `icv1-fp34k6-t1023235065186341755971623767326-5a99c63d` | GF(p^6), p 34 bits | 135 | 82.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.4 | -1.9 | +2.2 | formula_marginal |
| `icv1-fp18k14-tm40901129643050925884153966080961450506-2c66095b` | GF(p^14), p 18 bits | 215 | 123.3 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.0 | -0.7 | +1.2 | formula_marginal |
| `icv1-fp32k6-t39645023259520707621979274842-48f3071b` | GF(p^6), p 32 bits | 124 | 77.4 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 76.2 | -3.0 | +1.2 | formula_marginal |
| `icv1-fp34k6-tm800498701409058538299316831430-6ef9d975` | GF(p^6), p 34 bits | 131 | 81.0 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.3 | -3.5 | +0.7 | formula_marginal |
| `icv1-fp38k6-t3649156825249135653754775396351866-70a9f0de` | GF(p^6), p 38 bits | 147 | 89.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.6 | -3.5 | +0.7 | formula_marginal |
| `icv1-fp38k6-t3649156825249135653754775396351866-870dbe39` | GF(p^6), p 38 bits | 147 | 89.2 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 88.6 | -3.5 | +0.7 | formula_marginal |
| `icv1-fp18k14-tm14271219016216206938830791680039958586-6832ca0e` | GF(p^14), p 18 bits | 213 | 122.5 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 122.0 | -1.5 | +0.5 | formula_marginal |
| `icv1-fp34k6-tm4889968707445081258109231685378-bb7e12ef` | GF(p^6), p 34 bits | 131 | 80.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.3 | -3.7 | +0.4 | formula_marginal |
| `icv1-fp34k6-t9111335029992790914857964064330-15ce9d59` | GF(p^6), p 34 bits | 131 | 80.7 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.3 | -3.8 | +0.4 | formula_marginal |
| `icv1-fp34k6-tm783439549876049572437540640018-4d34d6b2` | GF(p^6), p 34 bits | 131 | 80.6 | descent to a genus-9 cover C/GF(p) of E/GF(p^6) | 9 | 80.3 | -3.9 | +0.3 | formula_marginal |
| `icv1-fp14k14-t20646852030805506221878730122-0929919d` | GF(p^14), p 14 bits | 186 | 108.0 | descent to a genus-41 cover C/GF(p) of E/GF(p^14) | 41 | 107.9 | +8.0 | +0.1 | formula_marginal |
| c2pnb304w1 (`icv1-f2m304-tm1949050066876782037155448244017940336576672053-02fbe089`) | GF(2^304) | 288 | 160.2 | GHS descent to Jac(C)/GF(2^76), N = 4 * 76 | 8 | 160.3 | +11.0 | -0.0 | boundary |
| c2pnb208w1 (`icv1-f2m208-t32421385432645202936390897948121-225f5036`) | GF(2^208) | 192 | 111.4 | GHS descent to Jac(C)/GF(2^104), N = 2 * 104 | 2 | 113.8 | +0.0 | -2.5 | boundary |
| secp112r2 (`icv1-fp112-t72213667414400864-3ecbdc5f`) | GF(p), 112 bits | 110 | 68.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 72.9 | +0.0 | -4.0 | boundary |
| wap-wsg-idm-ecid-wtls8 (`icv1-fp112-tm138743339785140737-cbe1d29c`) | GF(p), 112 bits | 112 | 69.9 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 74.0 | +0.0 | -4.0 | boundary |
| secp112r1 (`icv1-fp112-tm4407293269000505-7a11ff60`) | GF(p), 112 bits | 112 | 69.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 73.9 | +0.0 | -4.0 | boundary |
| secp128r2 (`icv1-fp128-t19023044537425668468-618c7c52`) | GF(p), 128 bits | 126 | 77.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 81.3 | +0.0 | -4.0 | boundary |
| secp128r1 (`icv1-fp128-tm8476633335676313877-1c931050`) | GF(p), 128 bits | 128 | 78.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 82.3 | +0.0 | -4.0 | boundary |
| brainpoolP160t1 (`icv1-fp160-t519972310379544251229703-a481858f`) | GF(p), 160 bits | 160 | 94.7 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 98.7 | +0.0 | -4.0 | boundary |
| brainpoolP160r1 (`icv1-fp160-t519972310379544251229703-b8575dce`) | GF(p), 160 bits | 160 | 94.7 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 98.7 | +0.0 | -4.0 | boundary |
| secp160k1 (`icv1-fp160-tm2082454586705745521609279-77aaac0a`) | GF(p), 160 bits | 160 | 94.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 98.8 | +0.0 | -4.0 | boundary |
| wap-wsg-idm-ecid-wtls9 (`icv1-fp160-tm2180728751409538655993507-410c12a9`) | GF(p), 160 bits | 160 | 94.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 98.8 | +0.0 | -4.0 | boundary |
| secp160r1 (`icv1-fp160-tm2364890543760230858433111-6b15890e`) | GF(p), 160 bits | 160 | 94.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 98.8 | +0.0 | -4.0 | boundary |
| secp160r2 (`icv1-fp160-tm250855509121001538581751-eed0581b`) | GF(p), 160 bits | 160 | 94.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 98.8 | +0.0 | -4.0 | boundary |
| prime192v2 (`icv1-fp192-t128840479891808162805939905231-7976b179`) | GF(p), 192 bits | 192 | 111.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 115.2 | +0.0 | -4.0 | boundary |
| secp192k1 (`icv1-fp192-t146402144145231529258894028971-0a61824b`) | GF(p), 192 bits | 192 | 111.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 115.2 | +0.0 | -4.0 | boundary |
| P-192 (`icv1-fp192-t31607402316713927207482677199-52e4af59`) | GF(p), 192 bits | 192 | 111.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 115.2 | +0.0 | -4.0 | boundary |
| prime192v3 (`icv1-fp192-t41351533396743040039645025261-5495245e`) | GF(p), 192 bits | 192 | 111.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 115.2 | +0.0 | -4.0 | boundary |
| brainpoolP192t1 (`icv1-fp192-t75885465139255996133178324439-0e200e9d`) | GF(p), 192 bits | 192 | 111.0 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 115.0 | +0.0 | -4.0 | boundary |
| brainpoolP192r1 (`icv1-fp192-t75885465139255996133178324439-bdd2bc3f`) | GF(p), 192 bits | 192 | 111.0 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 115.0 | +0.0 | -4.0 | boundary |
| brainpoolP224t1 (`icv1-fp224-t4460773185803614235113038911974753-d0798e5e`) | GF(p), 224 bits | 224 | 127.4 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 131.4 | +0.0 | -4.0 | boundary |
| brainpoolP224r1 (`icv1-fp224-t4460773185803614235113038911974753-e1f4ff20`) | GF(p), 224 bits | 224 | 127.4 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 131.4 | +0.0 | -4.0 | boundary |
| P-224 (`icv1-fp224-t4733100108545601916421827343930821-fec01e99`) | GF(p), 224 bits | 224 | 127.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 131.6 | +0.0 | -4.0 | boundary |
| secp224k1 (`icv1-fp224-tm9672873182660579502067895643393161-19ec5a97`) | GF(p), 224 bits | 224 | 127.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 131.6 | +0.0 | -4.0 | boundary |
| prime239v1 (`icv1-fp239-t506926255561332610517105435906892533-3a5a5c8a`) | GF(p), 239 bits | 239 | 135.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 139.2 | +0.0 | -4.0 | boundary |
| prime239v3 (`icv1-fp239-t543286227487037812274483689551665839-dbcebd3e`) | GF(p), 239 bits | 239 | 135.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 139.2 | +0.0 | -4.0 | boundary |
| prime239v2 (`icv1-fp239-tm1078211060294033492355314871841595491-312285e0`) | GF(p), 239 bits | 239 | 135.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 139.2 | +0.0 | -4.0 | boundary |
| W-25519 (`icv1-fp255-tm221938542218978828286815502327069187962-36f15f4b`) | GF(p), 255 bits | 252 | 141.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.9 | +0.0 | -4.0 | boundary |
| Curve25519 (`icv1-fp255-tm221938542218978828286815502327069187962-5eb6d947`) | GF(p), 255 bits | 252 | 141.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.9 | +0.0 | -4.0 | boundary |
| Edwards25519 (`icv1-fp255-tm221938542218978828286815502327069187962-f864fecb`) | GF(p), 255 bits | 252 | 141.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.9 | +0.0 | -4.0 | boundary |
| id-GostR3410-2001-CryptoPro-C-ParamSet (`icv1-fp256-t158633982881099435961745385353335106019-09cbae2d`) | GF(p), 255 bits | 255 | 143.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.5 | +0.0 | -4.0 | boundary |
| SM2 (`icv1-fp256-t188730266966446886577384576996245946077-80e68e0d`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.9 | +0.0 | -4.0 | boundary |
| id-GostR3410-2001-CryptoPro-A-ParamSet (`icv1-fp256-t196221757141320814079428202834060788997-51c05000`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.9 | +0.0 | -4.0 | boundary |
| brainpoolP256r1 (`icv1-fp256-t300418416528525664980082381967979838673-6ad773ce`) | GF(p), 255 bits | 255 | 143.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.6 | +0.0 | -4.0 | boundary |
| brainpoolP256t1 (`icv1-fp256-t300418416528525664980082381967979838673-bd148c7d`) | GF(p), 255 bits | 255 | 143.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.6 | +0.0 | -4.0 | boundary |
| secp256k1 (`icv1-fp256-t432420386565659656852420866390673177327-76dadd18`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.9 | +0.0 | -4.0 | boundary |
| P-256 (`icv1-fp256-t89188191154553853111372247798585809583-f188c491`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.9 | +0.0 | -4.0 | boundary |
| FRP256v1 (`icv1-fp256-tm35197163533674495870879651530057169373-68fd7b2c`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.8 | +0.0 | -4.0 | boundary |
| gost256 (`icv1-fp256-tm447942250905168602652181588181066248577-f512821f`) | GF(p), 255 bits | 255 | 143.3 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.4 | +0.0 | -4.0 | boundary |
| id-GostR3410-2001-CryptoPro-B-ParamSet (`icv1-fp256-tm467140827435544800340052492131167243509-c96d2aed`) | GF(p), 255 bits | 255 | 143.3 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 147.4 | +0.0 | -4.0 | boundary |
| id-tc26-gost-3410-2012-256-paramSetA (`icv1-fp256-tm84256526728449730591029627228991796228-3ff7fff5`) | GF(p), 256 bits | 254 | 142.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 146.9 | +0.0 | -4.0 | boundary |
| brainpoolP320r1 (`icv1-fp320-t1829129012274291271856083956306349739816241699607-60ebe171`) | GF(p), 320 bits | 320 | 176.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 180.2 | +0.0 | -4.0 | boundary |
| brainpoolP320t1 (`icv1-fp320-t1829129012274291271856083956306349739816241699607-fbf419a7`) | GF(p), 320 bits | 320 | 176.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 180.2 | +0.0 | -4.0 | boundary |
| P-384 (`icv1-fp384-t1388124618062372383606759648309780106643088307173319169677-cf66c6c9`) | GF(p), 384 bits | 384 | 208.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 212.8 | +0.0 | -4.0 | boundary |
| brainpoolP384t1 (`icv1-fp384-t5973228999478432667446284273866865475630836211106699249391-05fe695b`) | GF(p), 383 bits | 383 | 208.3 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 212.4 | +0.0 | -4.0 | boundary |
| brainpoolP384r1 (`icv1-fp384-t5973228999478432667446284273866865475630836211106699249391-975ab6a0`) | GF(p), 383 bits | 383 | 208.3 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 212.4 | +0.0 | -4.0 | boundary |
| Curve448 (`icv1-fp448-t28312320572429821613362531907042076847709625476988141958474579766324-3537278c`) | GF(p), 448 bits | 446 | 240.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 244.2 | +0.0 | -4.0 | boundary |
| E448 (`icv1-fp448-t28312320572429821613362531907042076847709625476988141958474579766324-3d2274d2`) | GF(p), 448 bits | 446 | 240.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 244.2 | +0.0 | -4.0 | boundary |
| Edwards448 (`icv1-fp448-t28312320572429821613362531907042076847709625476988141958474579766324-53658426`) | GF(p), 448 bits | 446 | 240.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 244.2 | +0.0 | -4.0 | boundary |
| W-448 (`icv1-fp448-t28312320572429821613362531907042076847709625476988141958474579766324-649c48bd`) | GF(p), 448 bits | 446 | 240.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 244.2 | +0.0 | -4.0 | boundary |
| gost512 (`icv1-fp511-t33317690176989408428354063686744165122149707869071609401230178415964786752661-6e80cdd8`) | GF(p), 510 bits | 510 | 272.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 276.5 | +0.0 | -4.0 | boundary |
| brainpoolP512t1 (`icv1-fp512-t133911538952573548431982907995132016398065354869678696041884562716142492272779-c7f8a111`) | GF(p), 511 bits | 511 | 273.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 277.2 | +0.0 | -4.0 | boundary |
| brainpoolP512r1 (`icv1-fp512-t133911538952573548431982907995132016398065354869678696041884562716142492272779-e69cd968`) | GF(p), 511 bits | 511 | 273.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 277.2 | +0.0 | -4.0 | boundary |
| id-tc26-gost-3410-12-512-paramSetA (`icv1-fp512-t97744483583712349266929640403245629889151353128602905529915952558174263790419-b856c63a`) | GF(p), 512 bits | 512 | 273.4 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 277.5 | +0.0 | -4.0 | boundary |
| id-tc26-gost-3410-2012-512-paramSetC (`icv1-fp512-t98513328673960770018330014639034783119468794885417684309186204259789288992276-0206fe7a`) | GF(p), 512 bits | 510 | 272.4 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 276.5 | +0.0 | -4.0 | boundary |
| id-tc26-gost-3410-12-512-paramSetB (`icv1-fp512-tm149097018916599106931901217895018650162542083271457414743815309640053273077069-aac01b1f`) | GF(p), 511 bits | 511 | 272.9 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 277.0 | +0.0 | -4.0 | boundary |
| P-521 (`icv1-fp521-t657877501894328237357444332315020117536923257219387276263472201219398408051703-f767e0d1`) | GF(p), 521 bits | 521 | 278.0 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 282.0 | +0.0 | -4.0 | boundary |
| wap-wsg-idm-ecid-wtls1 (`icv1-f2m113-t1267584991505179-9318c53d`) | GF(2^113) | 112 | 70.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^113) (ordinary_binary_cubic_pullback) | 3 | 75.0 | +0.0 | -5.1 | boundary |
| sect113r1 (`icv1-f2m113-tm122610772499221213-97df4ac6`) | GF(2^113) | 112 | 70.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^113) (ordinary_binary_cubic_pullback) | 3 | 75.0 | +0.0 | -5.1 | boundary |
| sect113r2 (`icv1-f2m113-tm148884003161464613-d77a0def`) | GF(2^113) | 112 | 70.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^113) (ordinary_binary_cubic_pullback) | 3 | 75.0 | +0.0 | -5.1 | boundary |
| ECC2K-130 (`icv1-f2m131-tm22283658519494248867-115e0dc5`) | GF(2^131) | 129 | 78.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^131) (ordinary_binary_cubic_pullback) | 3 | 83.9 | +0.0 | -5.1 | boundary |
| sect131r2 (`icv1-f2m131-tm52073227371480044317-ecd1b011`) | GF(2^131) | 130 | 79.3 | same-field genus-3 cover H -> E of degree 3 over GF(2^131) (ordinary_binary_cubic_pullback) | 3 | 84.4 | +0.0 | -5.1 | boundary |
| sect131r1 (`icv1-f2m131-tm80868651916585429657-35517338`) | GF(2^131) | 130 | 79.3 | same-field genus-3 cover H -> E of degree 3 over GF(2^131) (ordinary_binary_cubic_pullback) | 3 | 84.4 | +0.0 | -5.1 | boundary |
| Oakley Group 3 (`icv1-f2m155-tm414891960790832521345347-b865721c`) | GF(2^155) | 151 | 90.4 | same-field genus-3 cover H -> E of degree 3 over GF(2^155) (ordinary_binary_cubic_pullback) | 3 | 95.5 | +0.0 | -5.1 | boundary |
| sect163r1 (`icv1-f2m163-t1731532666194638619521227-c861cf8f`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 100.9 | +0.0 | -5.1 | boundary |
| c2pnb163v3 (`icv1-f2m163-t4581356680384720711965167-d8d162ce`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 100.9 | +0.0 | -5.1 | boundary |
| c2pnb163v2 (`icv1-f2m163-t4927277371730230283198131-5500578d`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 100.9 | +0.0 | -5.1 | boundary |
| wap-wsg-idm-ecid-wtls5 (`icv1-f2m163-tm4590722519955900204932993-eb7e4d1d`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 100.9 | +0.0 | -5.1 | boundary |
| sect163k1 (`icv1-f2m163-tm4845466632539410776804317-e15e581e`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 100.9 | +0.0 | -5.1 | boundary |
| B-163 (`icv1-f2m163-tm6224022517221266226059365-31d88d74`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 100.9 | +0.0 | -5.1 | boundary |
| DSTU4145-163 (`icv1-f2m163-tm6637329309425161870107289-653384de`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 100.9 | +0.0 | -5.1 | boundary |
| DSTU4145-167 (`icv1-f2m167-t744409597681349541568451-8d4c1932`) | GF(2^167) | 166 | 97.9 | same-field genus-3 cover H -> E of degree 3 over GF(2^167) (ordinary_binary_cubic_pullback) | 3 | 102.9 | +0.0 | -5.1 | boundary |
| DSTU4145-173 (`icv1-f2m173-tm118990531082834462588575939-e0f56426`) | GF(2^173) | 171 | 100.4 | same-field genus-3 cover H -> E of degree 3 over GF(2^173) (ordinary_binary_cubic_pullback) | 3 | 105.5 | +0.0 | -5.1 | boundary |
| DSTU4145-179 (`icv1-f2m179-t170443561202618430526296611-54da5c7e`) | GF(2^179) | 178 | 104.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^179) (ordinary_binary_cubic_pullback) | 3 | 109.1 | +0.0 | -5.1 | boundary |
| c2tnb191v3 (`icv1-f2m191-tm21751307834701491275992758225-17fe9713`) | GF(2^191) | 188 | 109.4 | same-field genus-3 cover H -> E of degree 3 over GF(2^191) (ordinary_binary_cubic_pullback) | 3 | 114.4 | +0.0 | -5.1 | boundary |
| c2tnb191v1 (`icv1-f2m191-tm2867709611233468671243809609-a71ef20d`) | GF(2^191) | 190 | 110.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^191) (ordinary_binary_cubic_pullback) | 3 | 115.2 | +0.0 | -5.1 | boundary |
| DSTU4145-191 (`icv1-f2m191-tm65396783579344579503573601949-86ea3d4e`) | GF(2^191) | 190 | 110.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^191) (ordinary_binary_cubic_pullback) | 3 | 115.2 | +0.0 | -5.1 | boundary |
| c2tnb191v2 (`icv1-f2m191-tm99424717553083670833495147979-0036f99e`) | GF(2^191) | 189 | 109.7 | same-field genus-3 cover H -> E of degree 3 over GF(2^191) (ordinary_binary_cubic_pullback) | 3 | 114.7 | +0.0 | -5.1 | boundary |
| sect193r1 (`icv1-f2m193-tm123763275178439483174011368593-bd4a68b5`) | GF(2^193) | 192 | 111.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^193) (ordinary_binary_cubic_pullback) | 3 | 116.3 | +0.0 | -5.1 | boundary |
| sect193r2 (`icv1-f2m193-tm214577892669885555631611261865-7d7aef6d`) | GF(2^193) | 192 | 111.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^193) (ordinary_binary_cubic_pullback) | 3 | 116.3 | +0.0 | -5.1 | boundary |
| sect233k1 (`icv1-f2m233-tm137381546011108235394987299651366779-08b072d7`) | GF(2^233) | 231 | 131.1 | same-field genus-3 cover H -> E of degree 3 over GF(2^233) (ordinary_binary_cubic_pullback) | 3 | 136.2 | +0.0 | -5.1 | boundary |
| B-233 (`icv1-f2m233-tm206777407530349254000433718821372333-2bb37fd4`) | GF(2^233) | 232 | 131.6 | same-field genus-3 cover H -> E of degree 3 over GF(2^233) (ordinary_binary_cubic_pullback) | 3 | 136.7 | +0.0 | -5.1 | boundary |
| DSTU4145-233 (`icv1-f2m233-tm206777407530349254000433718821372333-65002b7a`) | GF(2^233) | 232 | 131.6 | same-field genus-3 cover H -> E of degree 3 over GF(2^233) (ordinary_binary_cubic_pullback) | 3 | 136.7 | +0.0 | -5.1 | boundary |
| c2tnb239v3 (`icv1-f2m239-t1688252296967058035833017058395749899-677552d4`) | GF(2^239) | 236 | 133.5 | same-field genus-3 cover H -> E of degree 3 over GF(2^239) (ordinary_binary_cubic_pullback) | 3 | 138.6 | +0.0 | -5.1 | boundary |
| c2tnb239v2 (`icv1-f2m239-t775701775700413337265223005219721971-012b5b00`) | GF(2^239) | 236 | 133.9 | same-field genus-3 cover H -> E of degree 3 over GF(2^239) (ordinary_binary_cubic_pullback) | 3 | 139.0 | +0.0 | -5.1 | boundary |
| sect239k1 (`icv1-f2m239-tm1879124296849029916489117776215138963-fe13c2c6`) | GF(2^239) | 237 | 134.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^239) (ordinary_binary_cubic_pullback) | 3 | 139.3 | +0.0 | -5.1 | boundary |
| c2tnb239v1 (`icv1-f2m239-tm317806026669043444686057064183337243-1c7f8907`) | GF(2^239) | 237 | 134.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^239) (ordinary_binary_cubic_pullback) | 3 | 139.3 | +0.0 | -5.1 | boundary |
| DSTU4145-257 (`icv1-f2m257-tm549493087903001725763798196569983228979-1541e241`) | GF(2^257) | 255 | 143.3 | same-field genus-3 cover H -> E of degree 3 over GF(2^257) (ordinary_binary_cubic_pullback) | 3 | 148.4 | +0.0 | -5.1 | boundary |
| B-283 (`icv1-f2m283-t2863663306391796106224371145726066910599667-6a731228`) | GF(2^283) | 282 | 157.1 | same-field genus-3 cover H -> E of degree 3 over GF(2^283) (ordinary_binary_cubic_pullback) | 3 | 162.1 | +0.0 | -5.1 | boundary |
| sect283k1 (`icv1-f2m283-t7777244870872830999287791970962823977569917-ddd7db8c`) | GF(2^283) | 281 | 156.6 | same-field genus-3 cover H -> E of degree 3 over GF(2^283) (ordinary_binary_cubic_pullback) | 3 | 161.6 | +0.0 | -5.1 | boundary |
| DSTU4145-307 (`icv1-f2m307-t2833281535299466325038960383492079652813050003-0d180efa`) | GF(2^307) | 306 | 169.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^307) (ordinary_binary_cubic_pullback) | 3 | 174.3 | +0.0 | -5.1 | boundary |
| c2tnb359v1 (`icv1-f2m359-t892492397682618518080660149488378363069332596639535741-98bc4060`) | GF(2^359) | 353 | 193.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^359) (ordinary_binary_cubic_pullback) | 3 | 198.1 | +0.0 | -5.1 | boundary |
| DSTU4145-367 (`icv1-f2m367-tm29919614407377426422224252239009944793385357330657729169-5613f811`) | GF(2^367) | 366 | 199.7 | same-field genus-3 cover H -> E of degree 3 over GF(2^367) (ordinary_binary_cubic_pullback) | 3 | 204.7 | +0.0 | -5.1 | boundary |
| sect409k1 (`icv1-f2m409-t10457288737315625927447685387048320737638796957687575791173829-87b3ee45`) | GF(2^409) | 407 | 220.4 | same-field genus-3 cover H -> E of degree 3 over GF(2^409) (ordinary_binary_cubic_pullback) | 3 | 225.5 | +0.0 | -5.1 | boundary |
| B-409 (`icv1-f2m409-tm6059503967182126918765909026644927652236777310526686418445029-86b19751`) | GF(2^409) | 408 | 220.9 | same-field genus-3 cover H -> E of degree 3 over GF(2^409) (ordinary_binary_cubic_pullback) | 3 | 226.0 | +0.0 | -5.1 | boundary |
| c2tnb431r1 (`icv1-f2m431-t1799562222397931009847908291975708857846270858398803555164615329-2fc909a8`) | GF(2^431) | 418 | 225.9 | same-field genus-3 cover H -> E of degree 3 over GF(2^431) (ordinary_binary_cubic_pullback) | 3 | 231.0 | +0.0 | -5.1 | boundary |
| DSTU4145-431 (`icv1-f2m431-t57433707327821439195582946117660639423082714303598267113724573283-5a9bd0cb`) | GF(2^431) | 430 | 232.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^431) (ordinary_binary_cubic_pullback) | 3 | 237.1 | +0.0 | -5.1 | boundary |
| B-571 (`icv1-f2m571-t99534385013609758659469819150465382236412396452349171016760770327496674607579419075443-58357d9f`) | GF(2^571) | 570 | 302.7 | same-field genus-3 cover H -> E of degree 3 over GF(2^571) (ordinary_binary_cubic_pullback) | 3 | 307.7 | +0.0 | -5.1 | boundary |
| sect571k1 (`icv1-f2m571-tm148380926981691413899619140297051490364542574180493936232912339534208516828973111459843-3b150b80`) | GF(2^571) | 569 | 302.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^571) (ordinary_binary_cubic_pullback) | 3 | 307.2 | +0.0 | -5.1 | boundary |
| ECC2K-95 (`icv1-f2m97-t378251973071909-f27192f3`) | GF(2^97) | 95 | 61.1 | same-field genus-3 cover H -> E of degree 3 over GF(2^97) (ordinary_binary_cubic_pullback) | 3 | 66.2 | +0.0 | -5.1 | boundary |
