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
| `unresolved` | the two bounds of an invariant the screen could not pin down disagree: one beats rho, the other does not.  Today this happens only for the GHS magic number when `Tr(a) = 1` and `n` is even (below).  The verdict waits on computing that invariant. |
| `boundary` | the formulas themselves place the route at or above rho.  No constant-factor engineering on the route can cross rho, so the row is a boundary in the sense of `AGENTS.md` §1. |
| `not_carried` | the destination cannot contain the subgroup (`too_small`) or is `E` itself over a subfield (`degenerate`): the route solves a different problem. |

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
| `GF(p^k)` | `extension_index_calculus` (Gaudry–Diem, factor base from `GF(p)`) | — | `p^(2 − 2/k)` relations, each decomposition at the Bézout bound `2^(k(k−1))` of the Semaev system |
| `GF(p^k)` | `cover_descent` to a genus-`k` cover over `GF(p)` | whether `#E` is even (a rational 2-torsion point, the usual precondition at `k = 3`) | solver at `g = k` |

**The magic number and `a`.**  `E_(a,b) ≅ E_(a′,b)` whenever `a + a′` is of
the form `α² + α`, so `a` enters only through its absolute trace.  When the
trace is 0, or `n` is odd, `a` reduces into `F_2` and the magic number is
exactly the `b`-only span above, which is what
`research/cover_route_prime_binary_20261008/ghs_magic.rs` computed.  When the
trace is 1 and `n` is even, the reduced `a` can add one generator, so the
screen emits two rows for that tower, at `m_b` and `m_b + 1`, both
`heuristic`, and marks the pair `unresolved` when they disagree.  The
repository's `ec_trapdoor::audit_curve` adds Frobenius differences of the
unreduced `a` and can report a larger number; each row records it as
`audit_curve_magic_number` and is not priced with it.

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
| genus-`g` Jacobian operation | `10·g²` multiplications (Cantor is quadratic in `g`; `g = 1` recovers the elliptic constant) |
| reference: rho on `E` | `√(πr/4)` elliptic operations over `GF(q)` |
| `rho_on_jacobian` | `√(πr/4)` Jacobian operations over the destination field |
| `hyperelliptic_index_calculus` | `g! · q + g² · q^(2 − 2/g)` Jacobian operations: Gaudry's relation collection (a divisor splits over the degree-one factor base with probability about `1/g!`) plus double-large-prime linear algebra (Gaudry–Thomé–Thériault–Diem); `g ≥ 2`.  The `g!` term confines the method to genus small beside `q`. |
| `enge_gaudry_smoothness` | Enge–Gaudry priced by an explicit count, for genus up to 1024: factor base of the places of degree at most `b` (half the irreducibles of each degree), relations at the exact fraction of `b`-smooth monic polynomials of degree `g`, sparse linear algebra `g·F²` operations modulo `r`; minimised over `b` |
| `enge_gaudry` | `L_N(1/2, √2)` Jacobian operations, `N = q^g`; only above genus 1024 and when `g ≥ log2 q` |
| `extension_index_calculus` | `p^(2 − 2/k) · 2^(k(k−1))` elliptic operations over `GF(p^k)` |
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
  route's `evidence` level saying how much of it is established; an
  **`unresolved`** verdict names the invariant that decides it.  See
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
| c2pnb304w1 (`icv1-f2m304-tm1949050066876782037155448244017940336576672053-02fbe089`) | GF(2^304) | 288 | 160.2 | GHS descent to Jac(C)/GF(2^2), N = 152 * 2 [m = 9, upper bound] | 256 | 102.7 | +58.0 | +57.6 | unresolved |
| c2pnb272w1 (`icv1-f2m272-t114182707820475866940422329614334265958971-985955bd`) | GF(2^272) | 256 | 144.0 | GHS descent to Jac(C)/GF(2^2), N = 136 * 2 [m = 9, upper bound] | 256 | 102.6 | +42.0 | +41.3 | unresolved |
| c2pnb176w1 (`icv1-f2m176-t38161931374797611715595947-5639a6d4`) | GF(2^176) | 160 | 95.0 | GHS descent to Jac(C)/GF(2^2), N = 88 * 2 [m = 8, lower bound] | 127 | 71.1 | +21.9 | +23.9 | unresolved |
| Oakley Group 4 (`icv1-f2m185-t11157211627747266908830216517-8dbe23c4`) | GF(2^185) | 183 | 106.6 | GHS descent to Jac(C)/GF(2^37), N = 5 * 37 | 15 | 95.6 | +14.4 | +10.9 | formula_advantage |
| c2pnb368w1 (`icv1-f2m368-tm15466245949934416937841361495935610942842793802409489935-43ee62b3`) | GF(2^368) | 352 | 192.7 | GHS descent to Jac(C)/GF(2^92), N = 4 * 92 | 7 | 182.6 | +18.3 | +10.1 | formula_advantage |
| Oakley Group 3 (`icv1-f2m155-tm414891960790832521345347-b865721c`) | GF(2^155) | 151 | 90.4 | GHS descent to Jac(C)/GF(2^31), N = 5 * 31 | 15 | 89.2 | +5.4 | +1.1 | formula_marginal |
| `icv1-fp12k5-tm296780010-475f4271` | GF(p^5), p 11 bits | 39 | 32.0 | descent to a genus-5 cover C/GF(p) of E/GF(p^5) | 5 | 31.6 | -2.7 | +0.4 | formula_marginal |
| c2pnb208w1 (`icv1-f2m208-t32421385432645202936390897948121-225f5036`) | GF(2^208) | 192 | 111.4 | GHS descent to Jac(C)/GF(2^52), N = 4 * 52 | 7 | 112.7 | +6.9 | -1.4 | boundary |
| secp112r2 (`icv1-fp112-t72213667414400864-3ecbdc5f`) | GF(p), 112 bits | 110 | 68.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 70.8 | +0.0 | -2.0 | boundary |
| wap-wsg-idm-ecid-wtls8 (`icv1-fp112-tm138743339785140737-cbe1d29c`) | GF(p), 112 bits | 112 | 69.9 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 71.9 | +0.0 | -2.0 | boundary |
| secp112r1 (`icv1-fp112-tm4407293269000505-7a11ff60`) | GF(p), 112 bits | 112 | 69.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 71.8 | +0.0 | -2.0 | boundary |
| secp128r2 (`icv1-fp128-t19023044537425668468-618c7c52`) | GF(p), 128 bits | 126 | 77.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 79.2 | +0.0 | -2.0 | boundary |
| secp128r1 (`icv1-fp128-tm8476633335676313877-1c931050`) | GF(p), 128 bits | 128 | 78.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 80.2 | +0.0 | -2.0 | boundary |
| brainpoolP160t1 (`icv1-fp160-t519972310379544251229703-a481858f`) | GF(p), 160 bits | 160 | 94.7 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 96.7 | +0.0 | -2.0 | boundary |
| brainpoolP160r1 (`icv1-fp160-t519972310379544251229703-b8575dce`) | GF(p), 160 bits | 160 | 94.7 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 96.7 | +0.0 | -2.0 | boundary |
| secp160k1 (`icv1-fp160-tm2082454586705745521609279-77aaac0a`) | GF(p), 160 bits | 160 | 94.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 96.8 | +0.0 | -2.0 | boundary |
| wap-wsg-idm-ecid-wtls9 (`icv1-fp160-tm2180728751409538655993507-410c12a9`) | GF(p), 160 bits | 160 | 94.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 96.8 | +0.0 | -2.0 | boundary |
| secp160r1 (`icv1-fp160-tm2364890543760230858433111-6b15890e`) | GF(p), 160 bits | 160 | 94.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 96.8 | +0.0 | -2.0 | boundary |
| secp160r2 (`icv1-fp160-tm250855509121001538581751-eed0581b`) | GF(p), 160 bits | 160 | 94.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 96.8 | +0.0 | -2.0 | boundary |
| prime192v2 (`icv1-fp192-t128840479891808162805939905231-7976b179`) | GF(p), 192 bits | 192 | 111.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 113.2 | +0.0 | -2.0 | boundary |
| secp192k1 (`icv1-fp192-t146402144145231529258894028971-0a61824b`) | GF(p), 192 bits | 192 | 111.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 113.2 | +0.0 | -2.0 | boundary |
| P-192 (`icv1-fp192-t31607402316713927207482677199-52e4af59`) | GF(p), 192 bits | 192 | 111.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 113.2 | +0.0 | -2.0 | boundary |
| prime192v3 (`icv1-fp192-t41351533396743040039645025261-5495245e`) | GF(p), 192 bits | 192 | 111.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 113.2 | +0.0 | -2.0 | boundary |
| brainpoolP192t1 (`icv1-fp192-t75885465139255996133178324439-0e200e9d`) | GF(p), 192 bits | 192 | 111.0 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 113.0 | +0.0 | -2.0 | boundary |
| brainpoolP192r1 (`icv1-fp192-t75885465139255996133178324439-bdd2bc3f`) | GF(p), 192 bits | 192 | 111.0 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 113.0 | +0.0 | -2.0 | boundary |
| brainpoolP224t1 (`icv1-fp224-t4460773185803614235113038911974753-d0798e5e`) | GF(p), 224 bits | 224 | 127.4 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 129.4 | +0.0 | -2.0 | boundary |
| brainpoolP224r1 (`icv1-fp224-t4460773185803614235113038911974753-e1f4ff20`) | GF(p), 224 bits | 224 | 127.4 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 129.4 | +0.0 | -2.0 | boundary |
| P-224 (`icv1-fp224-t4733100108545601916421827343930821-fec01e99`) | GF(p), 224 bits | 224 | 127.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 129.5 | +0.0 | -2.0 | boundary |
| secp224k1 (`icv1-fp224-tm9672873182660579502067895643393161-19ec5a97`) | GF(p), 224 bits | 224 | 127.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 129.5 | +0.0 | -2.0 | boundary |
| prime239v1 (`icv1-fp239-t506926255561332610517105435906892533-3a5a5c8a`) | GF(p), 239 bits | 239 | 135.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 137.2 | +0.0 | -2.0 | boundary |
| prime239v3 (`icv1-fp239-t543286227487037812274483689551665839-dbcebd3e`) | GF(p), 239 bits | 239 | 135.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 137.2 | +0.0 | -2.0 | boundary |
| prime239v2 (`icv1-fp239-tm1078211060294033492355314871841595491-312285e0`) | GF(p), 239 bits | 239 | 135.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 137.2 | +0.0 | -2.0 | boundary |
| W-25519 (`icv1-fp255-tm221938542218978828286815502327069187962-36f15f4b`) | GF(p), 255 bits | 252 | 141.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 143.8 | +0.0 | -2.0 | boundary |
| Curve25519 (`icv1-fp255-tm221938542218978828286815502327069187962-5eb6d947`) | GF(p), 255 bits | 252 | 141.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 143.8 | +0.0 | -2.0 | boundary |
| Edwards25519 (`icv1-fp255-tm221938542218978828286815502327069187962-f864fecb`) | GF(p), 255 bits | 252 | 141.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 143.8 | +0.0 | -2.0 | boundary |
| id-GostR3410-2001-CryptoPro-C-ParamSet (`icv1-fp256-t158633982881099435961745385353335106019-09cbae2d`) | GF(p), 255 bits | 255 | 143.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.5 | +0.0 | -2.0 | boundary |
| SM2 (`icv1-fp256-t188730266966446886577384576996245946077-80e68e0d`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.8 | +0.0 | -2.0 | boundary |
| id-GostR3410-2001-CryptoPro-A-ParamSet (`icv1-fp256-t196221757141320814079428202834060788997-51c05000`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.8 | +0.0 | -2.0 | boundary |
| brainpoolP256r1 (`icv1-fp256-t300418416528525664980082381967979838673-6ad773ce`) | GF(p), 255 bits | 255 | 143.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.5 | +0.0 | -2.0 | boundary |
| brainpoolP256t1 (`icv1-fp256-t300418416528525664980082381967979838673-bd148c7d`) | GF(p), 255 bits | 255 | 143.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.5 | +0.0 | -2.0 | boundary |
| secp256k1 (`icv1-fp256-t432420386565659656852420866390673177327-76dadd18`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.8 | +0.0 | -2.0 | boundary |
| P-256 (`icv1-fp256-t89188191154553853111372247798585809583-f188c491`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.8 | +0.0 | -2.0 | boundary |
| FRP256v1 (`icv1-fp256-tm35197163533674495870879651530057169373-68fd7b2c`) | GF(p), 256 bits | 256 | 143.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.8 | +0.0 | -2.0 | boundary |
| gost256 (`icv1-fp256-tm447942250905168602652181588181066248577-f512821f`) | GF(p), 255 bits | 255 | 143.3 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.3 | +0.0 | -2.0 | boundary |
| id-GostR3410-2001-CryptoPro-B-ParamSet (`icv1-fp256-tm467140827435544800340052492131167243509-c96d2aed`) | GF(p), 255 bits | 255 | 143.3 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 145.3 | +0.0 | -2.0 | boundary |
| id-tc26-gost-3410-2012-256-paramSetA (`icv1-fp256-tm84256526728449730591029627228991796228-3ff7fff5`) | GF(p), 256 bits | 254 | 142.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 144.8 | +0.0 | -2.0 | boundary |
| brainpoolP320r1 (`icv1-fp320-t1829129012274291271856083956306349739816241699607-60ebe171`) | GF(p), 320 bits | 320 | 176.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 178.2 | +0.0 | -2.0 | boundary |
| brainpoolP320t1 (`icv1-fp320-t1829129012274291271856083956306349739816241699607-fbf419a7`) | GF(p), 320 bits | 320 | 176.2 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 178.2 | +0.0 | -2.0 | boundary |
| P-384 (`icv1-fp384-t1388124618062372383606759648309780106643088307173319169677-cf66c6c9`) | GF(p), 384 bits | 384 | 208.8 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 210.8 | +0.0 | -2.0 | boundary |
| brainpoolP384t1 (`icv1-fp384-t5973228999478432667446284273866865475630836211106699249391-05fe695b`) | GF(p), 383 bits | 383 | 208.3 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 210.3 | +0.0 | -2.0 | boundary |
| brainpoolP384r1 (`icv1-fp384-t5973228999478432667446284273866865475630836211106699249391-975ab6a0`) | GF(p), 383 bits | 383 | 208.3 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 210.3 | +0.0 | -2.0 | boundary |
| Curve448 (`icv1-fp448-t28312320572429821613362531907042076847709625476988141958474579766324-3537278c`) | GF(p), 448 bits | 446 | 240.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 242.1 | +0.0 | -2.0 | boundary |
| E448 (`icv1-fp448-t28312320572429821613362531907042076847709625476988141958474579766324-3d2274d2`) | GF(p), 448 bits | 446 | 240.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 242.1 | +0.0 | -2.0 | boundary |
| Edwards448 (`icv1-fp448-t28312320572429821613362531907042076847709625476988141958474579766324-53658426`) | GF(p), 448 bits | 446 | 240.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 242.1 | +0.0 | -2.0 | boundary |
| W-448 (`icv1-fp448-t28312320572429821613362531907042076847709625476988141958474579766324-649c48bd`) | GF(p), 448 bits | 446 | 240.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 242.1 | +0.0 | -2.0 | boundary |
| gost512 (`icv1-fp511-t33317690176989408428354063686744165122149707869071609401230178415964786752661-6e80cdd8`) | GF(p), 510 bits | 510 | 272.5 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 274.5 | +0.0 | -2.0 | boundary |
| brainpoolP512t1 (`icv1-fp512-t133911538952573548431982907995132016398065354869678696041884562716142492272779-c7f8a111`) | GF(p), 511 bits | 511 | 273.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 275.1 | +0.0 | -2.0 | boundary |
| brainpoolP512r1 (`icv1-fp512-t133911538952573548431982907995132016398065354869678696041884562716142492272779-e69cd968`) | GF(p), 511 bits | 511 | 273.1 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 275.1 | +0.0 | -2.0 | boundary |
| id-tc26-gost-3410-12-512-paramSetA (`icv1-fp512-t97744483583712349266929640403245629889151353128602905529915952558174263790419-b856c63a`) | GF(p), 512 bits | 512 | 273.4 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 275.4 | +0.0 | -2.0 | boundary |
| id-tc26-gost-3410-2012-512-paramSetC (`icv1-fp512-t98513328673960770018330014639034783119468794885417684309186204259789288992276-0206fe7a`) | GF(p), 512 bits | 510 | 272.4 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 274.4 | +0.0 | -2.0 | boundary |
| id-tc26-gost-3410-12-512-paramSetB (`icv1-fp512-tm149097018916599106931901217895018650162542083271457414743815309640053273077069-aac01b1f`) | GF(p), 511 bits | 511 | 272.9 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 274.9 | +0.0 | -2.0 | boundary |
| P-521 (`icv1-fp521-t657877501894328237357444332315020117536923257219387276263472201219398408051703-f767e0d1`) | GF(p), 521 bits | 521 | 278.0 | same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1) | 2 | 280.0 | +0.0 | -2.0 | boundary |
| wap-wsg-idm-ecid-wtls1 (`icv1-f2m113-t1267584991505179-9318c53d`) | GF(2^113) | 112 | 70.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^113) (ordinary_binary_cubic_pullback) | 3 | 73.1 | +0.0 | -3.2 | boundary |
| sect113r1 (`icv1-f2m113-tm122610772499221213-97df4ac6`) | GF(2^113) | 112 | 70.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^113) (ordinary_binary_cubic_pullback) | 3 | 73.1 | +0.0 | -3.2 | boundary |
| sect113r2 (`icv1-f2m113-tm148884003161464613-d77a0def`) | GF(2^113) | 112 | 70.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^113) (ordinary_binary_cubic_pullback) | 3 | 73.1 | +0.0 | -3.2 | boundary |
| ECC2K-130 (`icv1-f2m131-tm22283658519494248867-115e0dc5`) | GF(2^131) | 129 | 78.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^131) (ordinary_binary_cubic_pullback) | 3 | 82.0 | +0.0 | -3.2 | boundary |
| sect131r2 (`icv1-f2m131-tm52073227371480044317-ecd1b011`) | GF(2^131) | 130 | 79.3 | same-field genus-3 cover H -> E of degree 3 over GF(2^131) (ordinary_binary_cubic_pullback) | 3 | 82.5 | +0.0 | -3.2 | boundary |
| sect131r1 (`icv1-f2m131-tm80868651916585429657-35517338`) | GF(2^131) | 130 | 79.3 | same-field genus-3 cover H -> E of degree 3 over GF(2^131) (ordinary_binary_cubic_pullback) | 3 | 82.5 | +0.0 | -3.2 | boundary |
| sect163r1 (`icv1-f2m163-t1731532666194638619521227-c861cf8f`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 99.0 | +0.0 | -3.2 | boundary |
| c2pnb163v3 (`icv1-f2m163-t4581356680384720711965167-d8d162ce`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 99.0 | +0.0 | -3.2 | boundary |
| c2pnb163v2 (`icv1-f2m163-t4927277371730230283198131-5500578d`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 99.0 | +0.0 | -3.2 | boundary |
| wap-wsg-idm-ecid-wtls5 (`icv1-f2m163-tm4590722519955900204932993-eb7e4d1d`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 99.0 | +0.0 | -3.2 | boundary |
| sect163k1 (`icv1-f2m163-tm4845466632539410776804317-e15e581e`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 99.0 | +0.0 | -3.2 | boundary |
| B-163 (`icv1-f2m163-tm6224022517221266226059365-31d88d74`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 99.0 | +0.0 | -3.2 | boundary |
| DSTU4145-163 (`icv1-f2m163-tm6637329309425161870107289-653384de`) | GF(2^163) | 162 | 95.8 | same-field genus-3 cover H -> E of degree 3 over GF(2^163) (ordinary_binary_cubic_pullback) | 3 | 99.0 | +0.0 | -3.2 | boundary |
| DSTU4145-167 (`icv1-f2m167-t744409597681349541568451-8d4c1932`) | GF(2^167) | 166 | 97.9 | same-field genus-3 cover H -> E of degree 3 over GF(2^167) (ordinary_binary_cubic_pullback) | 3 | 101.0 | +0.0 | -3.2 | boundary |
| DSTU4145-173 (`icv1-f2m173-tm118990531082834462588575939-e0f56426`) | GF(2^173) | 171 | 100.4 | same-field genus-3 cover H -> E of degree 3 over GF(2^173) (ordinary_binary_cubic_pullback) | 3 | 103.6 | +0.0 | -3.2 | boundary |
| DSTU4145-179 (`icv1-f2m179-t170443561202618430526296611-54da5c7e`) | GF(2^179) | 178 | 104.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^179) (ordinary_binary_cubic_pullback) | 3 | 107.2 | +0.0 | -3.2 | boundary |
| c2tnb191v3 (`icv1-f2m191-tm21751307834701491275992758225-17fe9713`) | GF(2^191) | 188 | 109.4 | same-field genus-3 cover H -> E of degree 3 over GF(2^191) (ordinary_binary_cubic_pullback) | 3 | 112.5 | +0.0 | -3.2 | boundary |
| c2tnb191v1 (`icv1-f2m191-tm2867709611233468671243809609-a71ef20d`) | GF(2^191) | 190 | 110.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^191) (ordinary_binary_cubic_pullback) | 3 | 113.3 | +0.0 | -3.2 | boundary |
| DSTU4145-191 (`icv1-f2m191-tm65396783579344579503573601949-86ea3d4e`) | GF(2^191) | 190 | 110.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^191) (ordinary_binary_cubic_pullback) | 3 | 113.3 | +0.0 | -3.2 | boundary |
| c2tnb191v2 (`icv1-f2m191-tm99424717553083670833495147979-0036f99e`) | GF(2^191) | 189 | 109.7 | same-field genus-3 cover H -> E of degree 3 over GF(2^191) (ordinary_binary_cubic_pullback) | 3 | 112.8 | +0.0 | -3.2 | boundary |
| sect193r1 (`icv1-f2m193-tm123763275178439483174011368593-bd4a68b5`) | GF(2^193) | 192 | 111.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^193) (ordinary_binary_cubic_pullback) | 3 | 114.4 | +0.0 | -3.2 | boundary |
| sect193r2 (`icv1-f2m193-tm214577892669885555631611261865-7d7aef6d`) | GF(2^193) | 192 | 111.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^193) (ordinary_binary_cubic_pullback) | 3 | 114.4 | +0.0 | -3.2 | boundary |
| sect233k1 (`icv1-f2m233-tm137381546011108235394987299651366779-08b072d7`) | GF(2^233) | 231 | 131.1 | same-field genus-3 cover H -> E of degree 3 over GF(2^233) (ordinary_binary_cubic_pullback) | 3 | 134.3 | +0.0 | -3.2 | boundary |
| B-233 (`icv1-f2m233-tm206777407530349254000433718821372333-2bb37fd4`) | GF(2^233) | 232 | 131.6 | same-field genus-3 cover H -> E of degree 3 over GF(2^233) (ordinary_binary_cubic_pullback) | 3 | 134.8 | +0.0 | -3.2 | boundary |
| DSTU4145-233 (`icv1-f2m233-tm206777407530349254000433718821372333-65002b7a`) | GF(2^233) | 232 | 131.6 | same-field genus-3 cover H -> E of degree 3 over GF(2^233) (ordinary_binary_cubic_pullback) | 3 | 134.8 | +0.0 | -3.2 | boundary |
| c2tnb239v3 (`icv1-f2m239-t1688252296967058035833017058395749899-677552d4`) | GF(2^239) | 236 | 133.5 | same-field genus-3 cover H -> E of degree 3 over GF(2^239) (ordinary_binary_cubic_pullback) | 3 | 136.7 | +0.0 | -3.2 | boundary |
| c2tnb239v2 (`icv1-f2m239-t775701775700413337265223005219721971-012b5b00`) | GF(2^239) | 236 | 133.9 | same-field genus-3 cover H -> E of degree 3 over GF(2^239) (ordinary_binary_cubic_pullback) | 3 | 137.0 | +0.0 | -3.2 | boundary |
| sect239k1 (`icv1-f2m239-tm1879124296849029916489117776215138963-fe13c2c6`) | GF(2^239) | 237 | 134.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^239) (ordinary_binary_cubic_pullback) | 3 | 137.3 | +0.0 | -3.2 | boundary |
| c2tnb239v1 (`icv1-f2m239-tm317806026669043444686057064183337243-1c7f8907`) | GF(2^239) | 237 | 134.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^239) (ordinary_binary_cubic_pullback) | 3 | 137.3 | +0.0 | -3.2 | boundary |
| DSTU4145-257 (`icv1-f2m257-tm549493087903001725763798196569983228979-1541e241`) | GF(2^257) | 255 | 143.3 | same-field genus-3 cover H -> E of degree 3 over GF(2^257) (ordinary_binary_cubic_pullback) | 3 | 146.5 | +0.0 | -3.2 | boundary |
| B-283 (`icv1-f2m283-t2863663306391796106224371145726066910599667-6a731228`) | GF(2^283) | 282 | 157.1 | same-field genus-3 cover H -> E of degree 3 over GF(2^283) (ordinary_binary_cubic_pullback) | 3 | 160.2 | +0.0 | -3.2 | boundary |
| sect283k1 (`icv1-f2m283-t7777244870872830999287791970962823977569917-ddd7db8c`) | GF(2^283) | 281 | 156.6 | same-field genus-3 cover H -> E of degree 3 over GF(2^283) (ordinary_binary_cubic_pullback) | 3 | 159.7 | +0.0 | -3.2 | boundary |
| DSTU4145-307 (`icv1-f2m307-t2833281535299466325038960383492079652813050003-0d180efa`) | GF(2^307) | 306 | 169.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^307) (ordinary_binary_cubic_pullback) | 3 | 172.4 | +0.0 | -3.2 | boundary |
| c2tnb359v1 (`icv1-f2m359-t892492397682618518080660149488378363069332596639535741-98bc4060`) | GF(2^359) | 353 | 193.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^359) (ordinary_binary_cubic_pullback) | 3 | 196.1 | +0.0 | -3.2 | boundary |
| DSTU4145-367 (`icv1-f2m367-tm29919614407377426422224252239009944793385357330657729169-5613f811`) | GF(2^367) | 366 | 199.7 | same-field genus-3 cover H -> E of degree 3 over GF(2^367) (ordinary_binary_cubic_pullback) | 3 | 202.8 | +0.0 | -3.2 | boundary |
| sect409k1 (`icv1-f2m409-t10457288737315625927447685387048320737638796957687575791173829-87b3ee45`) | GF(2^409) | 407 | 220.4 | same-field genus-3 cover H -> E of degree 3 over GF(2^409) (ordinary_binary_cubic_pullback) | 3 | 223.6 | +0.0 | -3.2 | boundary |
| B-409 (`icv1-f2m409-tm6059503967182126918765909026644927652236777310526686418445029-86b19751`) | GF(2^409) | 408 | 220.9 | same-field genus-3 cover H -> E of degree 3 over GF(2^409) (ordinary_binary_cubic_pullback) | 3 | 224.1 | +0.0 | -3.2 | boundary |
| c2tnb431r1 (`icv1-f2m431-t1799562222397931009847908291975708857846270858398803555164615329-2fc909a8`) | GF(2^431) | 418 | 225.9 | same-field genus-3 cover H -> E of degree 3 over GF(2^431) (ordinary_binary_cubic_pullback) | 3 | 229.0 | +0.0 | -3.2 | boundary |
| DSTU4145-431 (`icv1-f2m431-t57433707327821439195582946117660639423082714303598267113724573283-5a9bd0cb`) | GF(2^431) | 430 | 232.0 | same-field genus-3 cover H -> E of degree 3 over GF(2^431) (ordinary_binary_cubic_pullback) | 3 | 235.2 | +0.0 | -3.2 | boundary |
| B-571 (`icv1-f2m571-t99534385013609758659469819150465382236412396452349171016760770327496674607579419075443-58357d9f`) | GF(2^571) | 570 | 302.7 | same-field genus-3 cover H -> E of degree 3 over GF(2^571) (ordinary_binary_cubic_pullback) | 3 | 305.8 | +0.0 | -3.2 | boundary |
| sect571k1 (`icv1-f2m571-tm148380926981691413899619140297051490364542574180493936232912339534208516828973111459843-3b150b80`) | GF(2^571) | 569 | 302.2 | same-field genus-3 cover H -> E of degree 3 over GF(2^571) (ordinary_binary_cubic_pullback) | 3 | 305.3 | +0.0 | -3.2 | boundary |
| ECC2K-95 (`icv1-f2m97-t378251973071909-f27192f3`) | GF(2^97) | 95 | 61.1 | same-field genus-3 cover H -> E of degree 3 over GF(2^97) (ordinary_binary_cubic_pullback) | 3 | 64.3 | +0.0 | -3.2 | boundary |
