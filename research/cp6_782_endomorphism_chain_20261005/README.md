# CP6-782 G1: CM discriminant −339, and the 5·5·7 chain endomorphism verified on points

**Status: order-level facts reproduced natively; two chain endomorphisms constructed
and verified on real CP6-782 points by the autoresearcher's explicit builder; no
timing, no benchmark, no ECDLP claim.** Class label: scalar-multiplication
*construction*, not an attack (same discipline as
`research/cryptopro_b_glv_chain_20261005`).

## Question

A session reported, from a live run of the autoresearcher's endomorphism sweep:
"CP6-782 (the 782-bit curve Zexe used for BLS12-377 proofs) has CM discriminant
−339, class number 6, with its cheapest endomorphism 9+ω of degree 5²·7, matching
CryptoPro-B's chain shape." That claim was not recorded anywhere (not in this repo,
not in `crypto-autoresearcher/research/endosweep_20261005/sweep/`, which has 35
targets and no CP6-782). The question here is whether the sweep's **explicit
construction** verifies that chain on actual points of the curve, or whether it
was only an order-level fact.

## Curve

CP6-782 G1 as shipped in arkworks `ark-cp6-782` (`frozen_inputs/arkworks_cp6_782_g1.rs`,
`arkworks_cp6_782_fq.rs`, fetched 2026-10-05 from `arkworks-rs/algebra`):
`y² = x³ + 5x + b` over `F_p`, `p` 782 bits, scalar field `r` = the BLS12-377 base
field prime (377 bits), cofactor `h` 390 bits. Order pinning: `h·r` lies in the
Hasse interval, `(h·r)·G = O` for the arkworks generator, and `h·G` has prime
order `r`. The harness's own one-point check is skipped because `r < 4√p`; the
Frobenius-discriminant factorisation below is itself a strong consistency check
on `#E = h·r` (a wrong order would not give `−339 · f²`).

## Results

| fact | value | source |
|:--|:--|:--|
| `t² − 4p` | `−339 · f²`, `f` 388 bits | sweep (`cp6_782_sweep.json`, `discriminant_certificate`) |
| fundamental discriminant | **−339** | sweep; native `scan --discriminant=-339` agrees on the order |
| class number `h(−339)` | **6** | sweep; native scan: forms `(1,1,85) (3,3,29) (5,±1,17) (7,±5,13)` |
| geometric units | `±1` only | native scan: `geometric_unit_count: 2` |
| minimum non-scalar degree | **85** (`ω`, norm `5·17`) | native scan `minimum_degree_witness`; sweep `min_nonscalar_degree` |
| cheapest chain by the cost model | `4+ω`, degree `3·5·7 = 105`, 48.8 M | sweep `cheap_endomorphisms` |
| best GLV-2 configuration (total) | **`9+ω`, degree `5²·7 = 175`**: 2432 M vs 3651 M generic, **1.50× modelled** | sweep `best` |
| runner-up | `4+ω`, degree 105: 2440 M, 1.50× | sweep |

So the session's line is right on every count except one nuance: `9+ω` is the
cheapest *total* GLV-2 configuration (its LLL basis is one bit shorter: 189 vs
190 bits), not the cheapest endomorphism. The cheapest single map is `4+ω` of
degree `3·5·7`, and the two configurations are within 0.3% under the model. The
minimum non-scalar *degree* is 85, not 175; `ω` itself is rejected because a
17-isogeny step is priced at 63 M against the 5-5-7 chain's 53 M.

### Explicit construction on points: yes, for both chains

`harness/endosweep/explicit.py` (autoresearcher, Python; frozen output here)
computed the ℓ-division polynomials over `F_p` for ℓ ∈ {3, 5, 7}, found the
`F_p`-rational kernel polynomials, walked the isogeny class in Kohel/Vélu
kernel-polynomial form, closed the walk back to `j(E)`, composed with the
isomorphism to `E`, and checked the composite on a point of prime order `r`:

| element | degree | steps | matched element on points | eigenvalue | walks examined | GLV-2 check | wall |
|:--|--:|:--|:--|:--|--:|:--|--:|
| `9+ω` | 175 | 5, 5, 7 | `−10+ω` (= `−(9+ω̄)`, same norm) | `λ_ω − 10` | 2 | 4 scalars, max coeff 187 bits, Babai bound 189 bits, `n` 377 bits | 24 s |
| `4+ω` | 105 | 3, 5, 7 | `5−ω` (= `4+ω̄`) | `5 − λ_ω` | 2 | 4 scalars, max coeff 189 bits, Babai bound 190 bits | 84 s |

Both records (`frozen_inputs/cp6_782_chain_5_5_7.json`,
`cp6_782_chain_3_5_7.json`) carry the full j-walk, the intermediate curves
`(a, b)`, `λ_ω`, the matched element and the eigenvalue. `found: true` means
step 5 of the builder's contract held: the composite acts on the prime-order
subgroup as one of the predicted scalars, up to units and conjugation. Because
the conductor of `Z[π]` is 388 bits, this is also the evidence that `End(E)` is
the maximal order of `Q(√−339)` at 3, 5 and 7 (the rational ℓ-kernels and the
closed walk exist exactly as the maximal order predicts).

Wall times are Python on this Mac and are **not** performance claims.

## What this does and does not say

* It establishes a constructive 2-dimensional GLV decomposition on CP6-782 G1
  with a degree-175 (or degree-105) endomorphism, verified on points. A native
  projective evaluation at 782 bits, benchmarked against width-w NAF, is the
  measurement that would settle the 1.50×; the CryptoPro-B modules
  (`src/ecc/cryptopro_b_*.rs`) are 256-bit specific and would need a 782-bit
  (13-limb) field to port.
* For ECDLP search it changes nothing: the order has only the units `±1`, so
  there is no equivalence-class rho gain beyond negation. The native scan's
  decision tree records `extra_geometric_automorphisms: excluded_for_this_order`.
  The endomorphism is a *scalar-multiplication* speed-up usable inside a walk's
  step function, not a reduction of the search space.
* Not a registered curve: CP6-782 is not in `src/ecc/curve_zoo.rs` and nothing
  here adds it.

## Provenance and reproduction

Python provenance, frozen per `AGENTS.md` (nothing in this directory is executed
from this repository):

* Driver scripts, run from `/Volumes/SSD990/crypto-autoresearcher` against its
  `harness/endosweep/` (commit `efaa5e0b09`), left untracked in that repo at
  `research/endosweep_20261005/cp6_782/`:
  `cp6_chain.py` (sha256 `2fc1f37e…7017`), `cp6_chain_el.py` (`ea45a955…df45`),
  `cp6_sweep.py` (`216c4080…7ad5`).
* `frozen_inputs/SHA256SUMS` covers every frozen file.

Native, from this repository:

```sh
cargo build --manifest-path tools/endomorphism-search/Cargo.toml --locked --release
tools/endomorphism-search/target/release/endomorphism-search scan --discriminant=-339 --degree-bound 400
```

reproduces `native/scan_D-339_deg400.json`: class number 6, two units, minimum
non-scalar degree 85, the degree-175 elements `(a, b) = (9, 1), (−10, 1)` and
their conjugates, all with `map_status: abstract_ring_element` and
`curve_binding: not_bound_to_a_curve`. The native `probe` cannot bind them to
CP6-782 (field cap `p ≤ 16381`); the binding evidence is the frozen Python
record above.

Coordination: Conductor was unreachable during this work
(`127.0.0.1:8443: connection refused`), so no task card was attached; the
change is confined to this new directory.
