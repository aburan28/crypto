# Fixed-core affine batches retain the Boolean F5 row-label signature

The frozen native structural screen passed **all 18 discovery and all 18 fresh
holdout primary cells**. In every n=16/20/24, seed and affine-family batch of
32 systems, the largest exact selected-row signature class contains **32/32**
assignments. All 144 cells and 2,016 system evaluations completed, and the
native verifier regenerated every polynomial, selected-row bit, digest and
criterion counter. The worker and verifier source bytes were identical in
discovery and holdout. No timing or complete F5 batch cache was measured.

| Variables | Multiplier row labels per system | Selected | Selected / labels | Pruned | Koszul / Frobenius prunes | Different quadratic-core signatures across both splits | Primary discovery / holdout |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 948 | 870 | 0.918 | 78 | 66 / 12 | 6 | bridge only |
| 16 | 2,192 | 2,056 | 0.938 | 136 | 120 / 16 | 6 | 6/6 / 6/6 |
| 20 | 4,220 | 4,010 | 0.950 | 210 | 190 / 20 | 6 | 6/6 / 6/6 |
| 24 | 7,224 | 6,924 | 0.958 | 300 | 276 / 24 | 6 | 6/6 / 6/6 |

The table copies exact counts from the sealed raw records. Its units are
**row labels**, not independent relations or row-space ranks. Each split has
72 cells and 1,008 evaluations, including repeated base systems across
families and batch sizes. The 2,016 evaluations are not 2,016 independent
quadratic cores. The three discovery and three disjoint holdout cores at each
size have six distinct signatures in total, while every affine assignment
within one core retains its complete bitstream. Aggregate criterion work is
367,416 GF(2) word XORs in discovery and 320,208 in holdout; individual work
counts can change even when the prune signature does not.

For quadratic `f_j=q_j+a_j` at degree four, admitted multipliers have degree
at most two, so the quartic projection of `t*f_j` depends only on `q_j`.
This study checks the additional necessary condition left open by the merged
[graded-block study](https://github.com/aburan28/crypto/pull/1207): the
repository's Boolean-specific `F5Criterion::new` selects the **same labelled
rows** as affine coefficients change. Both invariants hold on these public
generated fixtures. The selected signatures differ between quadratic cores,
so a cache needs an exact core/signature guard. Equal signatures are not a
certificate that cached elimination, criterion recomputation or output
unpacking will be cheaper than fresh complete F5 calls.

The [discovery bundle](qualified_discovery_01/manifest.json) has manifest
SHA-256 `529255bac0bac72d59d77155fba6f50b476292748431f229eb440f33c0d69564`.
The [holdout bundle](qualified_holdout_01/manifest.json) has manifest SHA-256
`71064b2ee7d48b2dc2dcda7fcbdb98469a6b8a4335688f1430f0a79167dc324f`.
All 21 discovery and 24 holdout member hashes and complete native replays
passed after copying. The holdout binding points to the passing discovery,
and its worker/verifier SHA-256 digests exactly match discovery:
`977c5d14b9991a567f61a5cab87f8dabe52dbb39b32fb804151447046a20d35a`
and `b6fb54267d9c1ab4ddca0cedb42fefd54f4b164230d795d23f2daefe81039713`.
Four optimized Rust tests pass, including generated-core and affine-walk
contracts, selected-bit/count replay, tamper rejection and sealed failure
refusal. The producer and verifier call the same inherited Boolean-safe F5
criterion; this is source-bound replay, not an unaffiliated proof of its
criterion or a complete F5 solver benchmark.

**Decision:** the structural gate is met and a separate matched complete-F5
batch timing experiment is justified. It must charge criterion recomputation,
high-block setup, low-block updates, output materialization, fallbacks,
memory and destruction against the fastest correct same-binary F5 path.
Natural relation yield, complete index-calculus cost and Pollard-rho cost are
still null. No curve or key input, DLP crossover or cryptanalytic
breakthrough is established.
