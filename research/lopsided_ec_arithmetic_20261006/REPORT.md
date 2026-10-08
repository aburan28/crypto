# EC-arithmetic speedups from lopsided thin products: screening (2026-10-06)

Question and method in the [protocol](README.md). Rounds 1–4 covered the
index-calculus pipeline, the general envelope, and the refutation audit.
This round turns the technique on elliptic-curve arithmetic itself.
Status labels: **[proposed]** / **[derived]** / **[verified]** /
**[measured]** (none in this round).

## Headline

**No speedup: no EC kernel has the technique's shape.** The method needs
many wanted inner products over one shared thin middle with small-integer
entries at large N. Our kernels are constant-size field formulas,
sequential group operations, single inner products, or single-RHS
transforms — every one fails before size is even reached. The nearest
miss (bulletproof's `<a,b>`) is killed twice over, in closed form, below.

## Screening table [derived]

| # | Kernel (repo site) | Verdict | Blocking gate + reason |
|---|---|---|---|
| 1 | Point add/double (`src/ecc/point.rs`: `add`, `double`, affine variants) | Negative | **Scale + shape**: fixed rational formulas in a handful of field ops; no N, no inner products, nothing shared |
| 2 | Single scalar mul (`scalar_mul`, `_vartime`, `_ct`) | Negative | **Shape**: sequential double-and-add group ops; no shared middle across steps |
| 3 | Batch Schnorr verification (`src/ecc/schnorr.rs`) = MSM | Negative | **Ring**: the "multiplication" is a group action (scalar times point), not a ring product — no bilinear identity can act on it |
| 4 | IPA `<a,G>`, `<b,H>`; KZG commitments (`src/zk/bulletproofs.rs`, `kzg.rs`) | Negative | **Ring**: MSMs again — group ops, however inner-product-flavoured the notation |
| 5 | Bulletproof `<a,b>` (single inner product over field scalars) | Negative — nearest miss | **Shape then ring**: one wanted pair, so no amortization (quantified below); scalars are 256-bit mod-p, not small integers |
| 6 | Extension fields, Miller loop, final exp (`src/bls12_381/`) | Negative | **Scale + shape**: fixed O(1)-size bilinear maps; pairings sequential with no shared middle |
| 7 | Field multiplication (schoolbook and would-be Karatsuba sites) | Negative | **Scale**: single O(1) products; the technique is asymptotic in N and has nothing to be thin in |
| 8 | Batch transforms / NTT-adjacent code in zk | Negative | **Shape**: single-RHS transforms or full outputs with a non-thin middle |

## Why the nearest miss misses, exactly

Two closed-form arguments, either sufficient:

1. **One-shot reading bound.** The data-structure form (Theorem 24)
   quotes query time *after* preprocessing that has already read both
   matrices. For a single wanted entry (`|W| = 1`), total work is bounded
   below by reading the inputs, `Omega(D)` — the naive inner product.
   The `O(D^0.437)` query never beats it one-shot; the method's power is
   *entirely* shared preprocessing amortized over `|W| >> 1` wants. Our
   IPA computes exactly one `<a,b>`. There is nothing to amortize over.
2. **Entry-size bound.** Entries must be small integers (`N^{O(1)}`); the
   recursion's encoding cost explodes otherwise. Our scalars are 256-bit
   field elements for every batch size in view — hundreds of bits where
   the bound allows `O(log N)`. Even a hypothetical many-pair setting
   over `F_p` fails before size (Figure 1, right panel).

The deeper pattern: EC arithmetic's "products" are either field
multiplications inside O(1)-size formulas or group actions. Neither is a
large family of Z-inner-products sharing one thin middle — the only
object the Schoenhage-recursion pruning can grip.

## What would qualify (revisit hook) [proposed]

A qualifying kernel computes many (`|W| >> 1`) inner products over Z with
small entries, all sharing one thin middle `D`, with the wanted pairs
known before the computation. No inventoried kernel has this shape;
multi-pairings share no middle, batch inversion shares a product chain
rather than inner products, and batch transforms have a single RHS. If
such a kernel is ever named — e.g. a future batch proof system whose
verifier is literally a sparse wanted-set of small-integer inner
products — it earns a measurement protocol per round 2's promotion
rules, screened first against `table2.json`'s widest row. Until then,
this screen stands and needs no re-running.

A Sec-6-style computer search for batch-EC bilinear identities would
likewise need such a shape to target first; the search space "faster
formulas for one product" is Karatsuba/Toom-Cook territory, not thin
products. No search is proposed here.

## Boundary, method, ratio

- EC kernels are counted in field/group ops on fixed-size objects; the
  IC unit `S` does not apply and no cross-method claim is made. Any
  future positive would need matched baseline/candidate operation
  counts on identical inputs with identical outputs.
- No measurement rows: gate evaluations are not answers, and the
  one-shot/entry-size arguments are derivations from stated bounds.
- Falsification target: a kernel with the qualifying shape above. The
  table is the search that found none.
- Wall time appears nowhere.

## Graphs checked — no change [verified by inspection]

- Scoreboard, progress-timeline, leaderboard files, browser data, curve
  registry: no measurement landed; nothing references EC-arithmetic
  speedups. Figure 1 is a new visual for this round, not a
  canonical-graph update.

## Open checks [proposed]

1. PRs #1493, #1494, #1501 are still open on other branches; this round
   touches none of their files.
2. If a batch proof system with a thin-product-shaped verifier enters
   the tree, the revisit hook above names the exact screen it faces.

## Sources

- Paper: `https://arxiv.org/pdf/2610.06783` (v1); Theorems 24–25, Table 2
  (p.40); Sec. 6 (computer-search direction).
- Rounds 1–4: `research/lopsided_thin_product_20261006/` (merged PR
  #1485), `research/lopsided_other_speedups_20261006/` (PR #1493),
  `research/lopsided_implications_20261006/` (PR #1494),
  `research/lopsided_general_regime_20261006/` (PR #1501).
- Kernel sites: `src/ecc/point.rs`, `src/ecc/schnorr.rs`,
  `src/zk/bulletproofs.rs` (IPA ll. 230–260), `src/zk/kzg.rs`,
  `src/bls12_381/`.
