# Performance campaign results

The results ledger of the performance campaign on this repository: what was
changed, how it was measured, what it bought, what did not work, and what
the numbers cannot show.  It is written for a reader who was not there.  The
index it reports is defined in `PERFORMANCE_INDEX.md`; the class of change
is defined in `AGENTS.md` §3.

## Summary

The campaign performance index is **1.904x**: the per-kernel product, over
the four crypto PRs (#908, #915, #1167, #1172), of callgrind
instruction-count ratios on 133 kernels in eight areas, combined by the
index formula with equal area weights.  It comes from rewrites of hot paths
in existing code that leave the kernels' fingerprinted outputs and counted
units unchanged (one documented contract change, under "Decisions and open
choices", and three differences that no fingerprint sees, under "How it was
measured"): word-sized and fixed-width modular arithmetic (Montgomery,
Shoup, u64 Euclid), Jacobian coordinates with one inversion, keyed merges
and bit-packed state; the gain is concentrated in dlp (10.087x), pdp
(2.987x), relation (1.956x), field_ec (1.757x), sat (1.375x) and bool_gb
(1.210x), while fp_gb (1.000x) and gf2_la (1.004x) did not move.  About half
of the dlp figure (47% of its logarithm) is three #1167 kernel steps
(203.18x, 17.78x, 14.40x), which come predominantly from `ec-vartime-walks`:
an inversion-free Pohlig-Hellman digit search and word-sized variable-time EC
arithmetic for the cryptanalysis walks (see "Largest kernel speedups").  It
is an instruction-count index of fixed kernels in the **engineering** class
(`AGENTS.md` §3), not a wall-clock index and not an end-to-end speedup of any
attack or workload.  Scope: the changes merged through PR #1172; later PRs
append to this file.

## How it was measured

**Per-PR attribution.**  For each of the four PRs, the merge commit is built
and compared with its first parent (main just before the merge), both with
the same `perfbench` harness.  Only that PR's changes are in the comparison;
other people's commits on main are excluded.  All four PRs landed as merge
commits, so the first parent is well defined.  The campaign figure is the
per-kernel product of the four per-PR speedups; the area indices and the
index are then computed from those products.

The alternative, one comparison of the original baseline (`f9c77524` plus
the current harness) against `259ecbd5` (the head of #1167, and the campaign
tip when that comparison was run on 2026-10-01, before #1172 merged), mixes
in every non-campaign commit merged to main since `f9c77524`, and is invalid by
design.  Four kernel fingerprints differ in it, all because of
non-campaign commits: `bool_gb/f4_basis_m2_n17` and
`bool_gb/f4_basis_random_n18` already differ at `92ac774b` (main before
#908), and `fp_gb/tower_f4_kummer_m2_n12` and
`fp_gb/tower_f4_kummer_m3_n9` changed with `00983dfc` ("f4_fp_tower: stop a
step's elimination at a full echelon"), which arrived with a merge of main.
Even a kernel with no campaign change moved in that comparison:
`fp_gb/f4_pkm_kummer_m2_t5` read +1.7% Ir.

**Why instructions.**  `perfindex.py instr` runs each kernel once under
callgrind, restricted to the timed region, at one thread, and applies the
formula of `PERFORMANCE_INDEX.md` to the instruction count `Ir`: kernel
speedup is base `Ir` / candidate `Ir`, an area index is the geometric mean
over its kernels, and the index is the weighted geometric mean of the eight
area indices, equally weighted (`docs/perf/weights.json` gives every area
weight 1, normalised to 0.125 over the eight areas).  The count is
deterministic and immune to host load, which matters because the host was
shared and changed during the campaign (see "Caveats").  It is not perfectly
repeatable for every kernel: the base count of `dlp/ecc2k130_certify_m131_len12`
varies about 1% between runs because of random `HashMap` keys, and
`relation/residual_walk_lmw_b22_x3` varies by about 1.5k Ir (0.0006%).
Each instruction comparison is a single run (one round), so the 95% interval
in the per-PR tables is the point value and the A/A column reads n/a.

**Fingerprints.**  Every kernel returns a 64-bit fingerprint of everything it
computed.  If any kernel's fingerprint differs between the two sides, or
between two runs of one side, the comparison is INVALID and no index is
reported; a changed kernel fingerprint therefore invalidates a comparison,
whoever changed it.  All four per-PR comparisons are valid.

**Class of change.**  Every change is **engineering** in the sense of
`AGENTS.md` §3 (legitimate, bounded, and not a finding).  `AGENTS.md` §3
glosses the class as "`S` falls, ratio to the floor flat"; for these kernels
the operation count `S` does not change (`PERFORMANCE_INDEX.md`), the time
per counted operation does.  The kernels' outputs, orders and counted units
are unchanged as far as their `perfbench` fingerprints see them, except
where the notes record a documented contract change.  There is one:
`ec-vartime-walks` on composite moduli (see "Decisions and open choices").
A fingerprint pins only what its kernel hashes, and the write-ups say where
equivalence rests on tests instead: the `pq_sparse_solve` kernel's
fingerprint pins only column 320, so equivalence on every other column rests
on the differential tests; the `perfbench` fingerprint does not see
`RelationSystem::ops`, and nothing checked the `residual_walk` change against
the original insert until test commit `a42923f4`; for `ecc2k130_guard`,
behavioural equivalence is proven by tests and certificates, not only by the
kernel fingerprint.  The per-change write-ups also record differences that
no fingerprint sees: the derived `Debug` output of `CsrMatrix`
(`koblitz_sparse_la`) and of `FactorBase` (`residual_walk`) changed, and a
debug build of `residual_walk` no longer panics on the overflow in `sub_mod`
for p > 2^63 (release results are identical for all inputs).

**Per-change figures.**  The wall and `Ir` figures in "Accepted changes" are
per-change measurements: paired A/B/A' wall time at one thread unless noted,
and callgrind counts, against the campaign baselines or the tree before the
change.  The campaign baselines are `f9c77524` (the crypto merge-base with
main) with the then-current harness overlaid: `crypto-base-v4` (pdp
kernels), `v5` (adds relation and sat), `v6` (adds dlp, field_ec, fp_gb).
The notes do not state the base for every row.  These figures are a
different measurement from the per-PR index, with different bases, and do
not add or multiply into it.

## Results

| PR | merge commit | base (first parent) | head (second parent) | merged (UTC date) | performance index |
|:--|:--|:--|:--|:--|--:|
| #908 | `524d3877` | `92ac774b` | `83191df0` | 2026-09-29 | 1.201x |
| #915 | `717d0aaf` | `0ec2a0dc` | `c4369263` | 2026-09-29 | 1.135x |
| #1167 | `88ff556f` | `b7d3f9ca` | `259ecbd5` | 2026-10-01 | 1.367x |
| #1172 | `441293ee` | `848f145b` | `31635a73` | 2026-10-01 | 1.022x |
| **campaign (kernel product)** | | | | | **1.904x** |

Callgrind `Ir`, one thread, 133 kernels, every comparison valid.  Equal
weights mean an area contributes one eighth of the logarithm of the index
however many kernels or changes it has: #1172 merges one change (CDCL hot
paths in `sat.rs`), and reads 1.022x because its sat area index, 1.178x, is
one of eight.

### Area indices

| area | kernels | #908 | #915 | #1167 | #1172 | campaign |
|:--|--:|--:|--:|--:|--:|--:|
| bool_gb | 19 | 1.212x | 1.001x | 0.997x | 1.000x | **1.210x** |
| dlp | 10 | 1.021x | 2.405x | 4.122x | 0.997x | **10.087x** |
| field_ec | 24 | 1.763x | 0.997x | 1.000x | 1.000x | **1.757x** |
| fp_gb | 19 | 1.002x | 1.000x | 0.998x | 1.000x | **1.000x** |
| gf2_la | 19 | 1.005x | 1.004x | 0.995x | 1.000x | **1.004x** |
| pdp | 21 | 1.473x | 1.002x | 1.999x | 1.013x | **2.987x** |
| relation | 10 | 1.352x | 1.136x | 1.273x | 1.000x | **1.956x** |
| sat | 11 | 0.993x | 1.003x | 1.172x | 1.178x | **1.375x** |

The accepted changes include none aimed at fp_gb or gf2_la, and those two
areas read 1.000x and 1.004x.  dlp and relation have ten kernels each, so a
single large kernel moves their area index.

### Largest kernel speedups (campaign product)

| kernel | speedup |
|:--|--:|
| `dlp/pohlig_hellman_smooth47` | 238.36x |
| `pdp/build_system_m2_n23` | 124.56x |
| `field_ec/ecdsa_verify_secp256k1_x4` | 31.58x |
| `field_ec/point_scalar_mul_secp256k1_x8` | 31.54x |
| `field_ec/ecdsa_verify_p256_x4` | 30.55x |
| `field_ec/point_scalar_mul_p256_x8` | 30.48x |
| `dlp/aut_folded_rho_j0_p25` | 26.84x |
| `dlp/rho_dp_zp_multi_q34_x3` | 19.30x |
| `dlp/ecc2k130_certify_m131_len12` | 18.84x |
| `dlp/gaudry_schost_negation_demo32` | 18.60x |
| `dlp/rho_floyd_zp_q36` | 17.53x |
| `dlp/collab_walkers_demo32_x8` | 14.48x |
| `pdp/instantiate_m2_n23` | 11.24x |
| `pdp/enumerate_m3_n15` | 10.37x |
| `pdp/weil_descend_s4_sym_n17_l4` | 8.74x |
| `pdp/enumerate_m2_n23` | 8.37x |
| `sat/semaev_s4_encode_n19l6_x4` | 6.17x |
| `pdp/build_system_m3_n15` | 5.97x |
| `relation/pq_sparse_solve_m127_c640` | 4.68x |
| `pdp/subspace_oracle_n20_l7` | 3.59x |

The top row is the product of 1.173x (#908), 1.000x (#915), 203.18x (#1167,
349,778,070 to 1,721,536 `Ir`) and 1.000x (#1172).  In the same #1167
comparison `dlp/gaudry_schost_negation_demo32` reads 17.78x and
`dlp/collab_walkers_demo32_x8` 14.40x.  The figures are the measured
ratios.  These three #1167 steps come predominantly from the
`ec-vartime-walks` change (merge `8403dd16`; see "Accepted changes"), whose
kernel list is exactly these three kernels.  The #1167 window also contains
the `utils::mod_inverse` word path (`8925e497`), which all three kernels call,
so the per-PR figures are not that change alone.  The change measured on its
own, against the tree before it, gives 17.70x (`gaudry_schost_negation_demo32`),
14.31x (`collab_walkers_demo32_x8`) and 201.1x (`pohlig_hellman_smooth47`)
in the implementer's first round, and 17.71x, 14.35x and 201.0x in a later
review of the merged code against `b072fcf5`; the per-PR figures are 17.78x,
14.40x and 203.18x.  The 1.173x that `pohlig_hellman_smooth47`
shows in #908 comes from the one commit in that window that touches
`src/ecc`, the Jacobian `scalar_mul` change (`d4851a0a`).

## Accepted changes

Crypto repository, from the notes' "Accepted and merged" ledger.  `PR` is
the first of the four whose merge commit contains the change (from git
history).  `wall` is the median wall time of the kernel at one thread unless
noted, `Ir` is the callgrind instruction count, `CI` is the 95% interval of
the paired wall ratio, `load` is the host load average.  Where a row gives a
ratio without `wall` or `Ir`, the notes do not state the unit.  Rows are in
the notes' order.

| PR | commit | change | measured |
|:--|:--|:--|:--|
| #908 | `0cce6204` | `koblitz_groebner`: Macaulay readback by set bits | 1.35–1.68x wall, 1.57–1.95x Ir (bool_gb `matrix_f4_*`) |
| #908 | `2f41cedf` | `mq_fes`: blocked Moebius (two layers per pass) | 1.79x wall, 2.28x Ir (`fes_moebius`) |
| #908 | `d4851a0a` | `ecc::point`: Jacobian `scalar_mul`, one inversion | 27.7x P-256, 31.4x secp256k1 |
| #908 | `a5c02b47` | `mq_monica`: bit-packed state | 1.83x wall, 2.26x Ir |
| #908 | `52d7e6ef` | `crossbred`: word-major tables | 1.20x Ir |
| #908 | `fe17f3c1` | `pq_groebner_f2`: `from_monos` keyed sort (>16 terms) | `solve_inherited_m3` 1.09x Ir |
| #908 | `9c0ac2e6` | `pq_descent_symbolic`: split with one sort | descend 3.47x |
| #908 | `3d214fed` | `enumerate_decompose` on `FastCurve` (`FastEnumeration`) | enumerate 5.89x / 6.78x |
| #908 | `994784af` | `semaev_leading_form`: `MonoHasher` (folded multiply) | S5 2.11x wall / 3.10x Ir; S6 1.75x |
| #908 | `248af250` | `koblitz_sparse_la`: division-free (lazy lanes, Shoup) | `block_wiedemann` 2.93x (Ir 2.71x), `sparse_solve` 2.69x (Ir 2.79x); relation area 1.24x |
| #908 | `83191df0` | `m_can_decompose` on `FastCurve`; descent wrap | `koblitz_rung_k0n31` 3.77x wall (review run: CI 3.50-5.45, A/A ±8.9%, load about 9; the implementer's runs give about 3.1-3.3x), Ir 2.99x; descent 1.05x Ir, wall neutral |
| #908 | `b072fcf5` | clippy 1.98 lint fixes (CI) | no behaviour change |
| #915 | `d616c008` (+ `943ccfe8` nits) | `ecc2k130_guard`: fixed-width Montgomery certify search (after #908 merged) | `ecc2k130_certify_m131_len12` 17.5x wall (CI 17.2-19.6), Ir 19.1x; `ecc2k-guard` JSON byte-identical m=2..450 len 12 |
| #915 | `4fc58708` (+ `c4369263` nits) | `pollard_rho`: single-word Z_p^* walks (Floyd + DP multi) | `rho_dp_zp_multi` 19.7x wall (review run, CI 13.36-35.93; the implementer's run gives 24.5x), Ir 19.1x; `rho_floyd_zp_q36` 27.7x wall (review run, CI 21.49-31.92; the implementer's run gives 21.9x), Ir 17.9x. Open note: `dp_walks_word` allocates 520 B/target up front (only caller m=3) |
| #915 | `8f05c77b` | `residual_walk`: u64 Euclid, Shoup rows on support, FxMap, arena states | `residual_walk_lmw_b22_x3` 2.06x wall (CI 1.99-2.14), Ir 3.58x; 2.05-2.08x at 22/33/34 bits. `mul_mod` word branch wall-neutral (kept) |
| #1167 | `14290207` | `SubspaceOracle`: quartic reduction by t^4..t^6 table, Horner u->u^2, pclmul pair loop | `subspace_oracle_n20_l7` 2.82x wall quiet host (CI 2.79-3.09), Ir 3.59x; portable fallback 1.55-1.71x; aarch64 unmeasured |
| #1167 | `abf31579` (merge) | `pq_sparse_la`: pivot-only inverse + Vec counts; `pq_wiedemann`: delayed reduction | `pq_sparse_solve` 5.55x wall (Ir 4.66x); `pq_wiedemann` 2.74x (Ir 2.42x). Review 1 caught a `usize::MAX` column panic (fixed `ffa41f17` + regression test) |
| #1167 | `f6b44065` (merge) | `koblitz_relation_solver`: Shoup mulmod, word invmod, lead-bounded back-substitution | `relation_solver_u1024` 1.61x wall (CI 1.51-1.81), Ir 1.587x; wide path (m>=2^63) 1.50-1.54x standalone. Review 1: needs fix (prefetch dropped, `d14236c1`); review 2: accept |
| #1167 | `8925e497` (merge, + doc-link nit) | `utils::mod_inverse`: u64 extended Euclid; `cga_hnc` word group law; `aut_folded_rho` canonical-form carry | `aut_folded_rho_j0_p25` 19.5x wall (CI 12.3-21.6, load 7-9), Ir 26.85x. The notes also list `koblitz_signed_rho` 1.116x Ir and `pohlig_hellman_smooth47` 1.173x Ir, others neutral; this file does not credit those figures to this change (see "Where the sources disagree"). Review accept first round (+ `51834619` whole-stack differential test) |
| #1167 | `4542f7be` (merge) | `F2BoolPoly::add` keyed merge (popcnt dispatch), `DecompositionTemplate` bitsets, S3/`FieldSum` products over supports | Ir vs pre-change: `build_system_m2_n23` 124.8x, `build_system_m3_n15` 5.91x, `instantiate_m2_n23` 11.25x, `instantiate_m3_n15` 3.33x, `symmetrised_groebner_m3_n15` 1.42x; wall 96.6x / 4.96x / 11.75x / 2.31x / 1.20x. bool_gb ±0.2%. Review 1: needs fix (`57f39078`); review 2: accept. The `koblitz_symmetrised` gate test asserts a wall-clock ratio and failed once under load 13 (pre-existing, unrelated path) |
| #1167 | `4b26b717` | `pq_groebner_f2`: clippy 1.98 (`sort_unstable_by_key` + `Reverse`; same comparator) | no behaviour change |
| #1167 | `883ec267` (merge) | `binary_semaev_s4`: `AnfPoly`/`AnfF2m` xor as ordered merges; toggles kept for small sums | `weil_descend_s4_sym_n17_l4` 8.645x Ir, 9.69x wall (CI 8.18-14.48); `rank_refute_anf` neutral. Review 1: needs fix (small-shape wall), fix `fa6d176b`; review 2: accept |
| #1167 | `8403dd16` (merge, + `f1685a91`) | `ec-vartime-walks`: variable-time EC arithmetic for the cryptanalysis walks (Euclid `inv_vartime`) | Ir against the tree before the change (first round): `gaudry_schost_negation_demo32` 17.70x, `collab_walkers_demo32_x8` 14.31x, `pohlig_hellman_smooth47` 201.1x; against `f9c77524`'s build, final round: 18.54x / 14.43x / 235.4x Ir (the extra 1.17x on `pohlig_hellman_smooth47` is the #908 factor, from `d4851a0a`). Wall at one thread against `f9c77524`'s build, final round, host load 6-10, A/A spreads ±4.6% / ±27.9% / ±11.3%: 11.72x (CI 7.53-13.02) / 9.08x (CI 7.92-10.11) / 181.8x (CI 148.7-215.6); the other seven dlp kernels neutral. Per-step `Ir` against `f9c77524`'s build (Gaudry-Schost / collab / Pohlig-Hellman, first round): 2603M / 1728M / 410M before; 574M / 360M / 29M after `inv_vartime`; 257M / 170M / 3.3M after the one-word path and `scalar_mul_vartime`; the projective digit search takes Pohlig-Hellman to 1.74M; the walk bookkeeping and word helpers give 140.5M / 120.1M. Bit-identical on every prime-p job (reviewer's 338-solve / 185-DP transcripts), changed on composite p (see "Decisions and open choices"). Source: `results-ec-vartime-walks.txt` |
| #1172 | `31635a73` (merge, + `452b0b54`) | `sat-core`: CDCL hot paths in `sat.rs` | sat area 1.180x Ir in the notes (1.178x in the #1172 index); review 1 found an `add_vars` O(n) regression (fixed), review 2's out-of-range literal/trigger issue fixed in `452b0b54` |

The last two rows are merged and counted in the indices above but are not in
the notes' "Accepted and merged" table; they are taken from the notes'
decisions and PR history and from git.  The notes give no merge commit for
`pq-sparse-wiedemann` and `mod-inverse-u64`; `abf31579` and `8925e497` are
the merge commits found in git by branch name.

**C library and suite.**  The notes also record accepted changes to the C
library and the Rust suite in the separate cryptanalysis repository (PR
#168, merged there; later suite ports were merged in #205 and #207, and one
on the branch with no PR); they
are measured against that repository's own baselines and are not part of the
index above.  C library:
REDC subtraction form (Ir index 1.118x); `ic_fbdiv` inverse trial division
(2.57x / 1.92x wall; Ir 0.82x, the divide blind spot); Montgomery powmod +
Miller–Rabin (1.53x / 1.48x); `ca_ec_group_mul` Jacobian (4.16x);
`_Thread_local` factorize memo (`pohlig_zp` 1.98x); Shoup dense solve /
Lanczos (sparse 1.88x, IC precompute 1.51x); `ca_zp_group_mul` (Cheon
1.37x); BSGS 64-lane batch (EC 3.55x, `pohlig_ec` 3.87x); Mestre batched
(`count_points` 6.01x); grumpy 3-wide (1.94x); rho cached hash + fastmod (Ir
1.18x); kangaroo (Ir 1.17x); `ca_htab1` (`pohlig_zp` 1.33x).  Suite:
Jacobian point (Smart 6.9x); Montgomery `brent_u64` 1.29x, `brent_u128`
Mont128 16.4x; pm1 split 4.98x; `weak_curves` batched restarts 2.21x; QS
SWAR scan 1.34x; ports of `mq_monica`, `mq_fes`, `koblitz_groebner` and
`crossbred`; suite instruction index on touched kernels 4.95x.  Later suite
ports, from the notes: the 101 ported `perfbench` kernels (`75ac845`) have
fingerprints equal to the base; `koblitz_sparse_la` (`1f05ba7`, #205) gives
`block_wiedemann` 2.92x wall / 2.713x `Ir`, `sparse_solve` 2.60x / 2.790x and
relation area 1.110x against `suite-base-v3`; the `semaev_leading_form`
hasher (`2514bbc`, #207) gives S5 2.35x wall / 3.217x `Ir` and S6 2.39x /
3.061x.  `port-f2-sorts` is merged on the cryptanalysis branch as `fa9b807`
with no PR: `pdp/descend_s3_n31_np20` 2.36x wall / 3.46x `Ir`,
`descend_s4_n23_np6` 2.43x / 3.30x, pdp area 1.082x, overall 1.040x (CI
1.016-1.044) against `suite-base-v3`.

## Changes after #1172

Changes accepted after the four PRs above.  They are not in the 1.904x, which
covers #908 to #1172 only; each gets a per-PR index when its PR merges.

### `gf2-fold-reduce`: `Gf2` reduced by two carry-less folds

`Gf2::mul`, `sqr`, `sqr_k`, `inv` and `batch_inv` (`semaev_decomp.rs`,
`koblitz_fast.rs`) reduce a product by two carry-less folds with the sparse
modulus tail, when the tail is short (`2·deg(tail) <= n + 1`) and the CPU has
`pclmulqdq`.  The byte-table path stays for dense tails and for other CPUs.
Branch `perf/opt-gf2-fold-reduce`, merged into the crypto branch as
`809a7318`; the review verdict is accept.  All 133 `perfbench` fingerprints
are equal to a rebuilt `b072fcf5` control at one and four threads.  Source for
every figure here: `results-gf2-fold-reduce.txt`.

Review, interleaved min-of-samples wall time against that control, one thread,
host load 2.5-4 (ctl / cand): `gf2_n53_mul_sqr_1m` 2.44x, `gf2_n53_inv_16k`
2.49x, `gf2_n53_batch_inv_4x64k` 1.66x, `koblitz_n53_fast_mul_2k` 2.86x,
`koblitz_curve_mul_n53_1k` 3.61x, `pdp/mitm_m3_n31` 1.59x, `pdp/mitm_m4_n31`
1.67x, `pdp/pair_table_build_n31_l10` 1.64x, `pdp/subspace_oracle_n20_l7`
1.48x, `pdp/enumerate_m3_n15` 1.82x, `pdp/enumerate_m2_n23` 1.87x,
`dlp/koblitz_signed_rho_k0_n41` 1.76x, `relation/koblitz_collect_k0n41` 1.17x,
`relation/koblitz_pair_table_k0n41_folded` 1.28x,
`relation/koblitz_descent_k0n53_m2_x4` 1.11x, `relation/koblitz_collect_k0n53`
1.05x, `relation/koblitz_rung_k0n31` 1.01x; `koblitz_n53_add_many_lazy_64x2k`
0.99x, neutral, because this host runs its `Simd512` path.  Instruction counts
on the same pair: `gf2_n53_mul_sqr` 2.97x, `inv` 5.05x, `batch_inv` 2.53x,
`koblitz_n53_fast_mul` 6.27x, `koblitz_curve_mul_n53_1k` 8.98x.  AArch64
compiles the dispatch out and was not measured.

Measured on the merged tree, after main's own `SubspaceOracle` table change
(`14290207`) and other commits had landed, callgrind `Ir` against the build of
`441293ee`: `pdp/subspace_oracle_n20_l7` 192.3M to 139.6M (1.377x), so the two
changes compose; `pdp/enumerate_m2_n23` 1.92x, `pdp/pair_table_build_n31_l10`
2.15x, `pdp/mitm_m3_n31` 2.12x and `mitm_m4_n31` 2.23x, `dlp/koblitz_signed_rho_k0_n41`
2.44x, `relation/koblitz_collect_k0n53_w477_u1024` 1.39x.  Fingerprints are
equal on every kernel measured except `relation/koblitz_pair_table_k0n41_folded`,
whose fingerprint differs between main at `441293ee` and main at `ff5a31d06`
and is the same in the merged tree as in main at `ff5a31d06`; the difference
arrives with main's commits between those two (`3314760bc` is the one among them
that touches this kernel's code), not with this change.

Caveats the review recorded:

- **The fallback is not unchanged.**  For a field that does not fold (tail
  degree above `(n+1)/2`, or no `pclmulqdq`) the new `folds()` test makes
  independent-multiply throughput 7-49% slower than before (`n = 7` +39%,
  `n = 15` +49%, `n = 31` +24%, `n = 41` +14%, `n = 53` +7%, `n = 63` about
  0; `batch_inv` 5-12% slower; latency unchanged).  No pipeline or `perfbench`
  modulus is affected: `find_irreducible` and `find_irreducible_sparse` give
  short tails, and every `perfbench` kernel folds.  The commit messages on the
  branch say "unchanged"; the merge commit corrects that.
- **Layout shifts on kernels with no `Gf2` code.**  Against the control,
  `sat/random3sat_n150` +2.4% `Ir`, `random3sat_n200` +2.8%,
  `sat/semaev_s4_n15l5_cnf` +3.3%, `sat/koblitz_bool_cnf` +2.0% and
  `dlp/aut_folded_rho` +7.0%.  Interleaved wall time for these is within noise
  (ratios 0.93-1.05, no consistent sign), so this is code layout, the same
  effect as the `ct_scalar_mul` kernels below.

Two things happened when the change was integrated (the first in the merge
commit `809a7318`, the second in a later commit):

- **Merge conflict.**  Main had gained `Gf2::dot2` (the `SubspaceOracle`
  quartic route) at the spot where this change makes `Gf2::sqr` a
  `folds()`-dispatching wrapper.  Both are kept; `dot2` calls `Gf2::reduce`,
  which dispatches to the fold path, so `dot2` and the other new call sites
  pick the fold up unchanged.
- **One test narrowed.**  `general_quartic_matches_994784af_even_for_unreduced_b`
  (from the `SubspaceOracle` change) passes a curve constant `b` wider than the
  field (it sets bits up to 63; the failing case was `n = 20` with `b` about
2^58), so the product exceeds the
  documented domain of `reduce` (below 2^(2n-1)).  The byte table stays
  F2-linear there; the fold path promises nothing, and `Gf2::mul` on such an
  operand is not `reduce(clmul(a, b))` once it folds.  Production never
  passes such a `b`: `ic_boundary` builds it as `1` for a Koblitz instance
  and as `(rng.gen::<u64>() & gf.mask) | 2` for a random binary one, so it is
  below 2^n, and the other callers pass `1`.  (The commit message of
  `e378cae2` says `ic_boundary` draws `b` from `1..short` and that the
  narrowing is disclosed in the merge commit; both are wrong: those draws are
  walk scalars, and the disclosure is here and in the PR description.)  The
  fold review asserted `operand < 2^n` in `mul`,
  `sqr`, `sqr_k`, `inv` and `batch_inv` over the whole library suite with no
  failure.  The off-domain case now runs on table fields in release builds
  only; every in-domain case still runs on every field (commit `e378cae2`).
  This narrows a test that another change had added, so it is stated here and
  in the PR.

## Caveats

### Hosts and wall time

The container moved between two host types during the campaign.

| host | CPU | relevant features | role in the campaign |
|:--|:--|:--|:--|
| H1 | Intel Xeon @ 2.10 GHz, 4 vCPU, 16 GB | AVX-512F/BW/VL/IFMA/VBMI2, VPCLMULQDQ, GFNI | `crypto-base-v3`, `suite-base-v1/v2`, `clib-base-v1` built here |
| H2 | Intel Xeon @ 2.80 GHz (Cascade Lake), 4 vCPU | AVX-512F/BW/VL/DQ/VNNI, pclmulqdq, no VPCLMULQDQ | `crypto-base-v4` to `-v6` built here; the notes also call it the host of most batch A/B measurements (see "Where the sources disagree") |

Binaries carry no `target-cpu` flags, so a baseline built on one host runs
unchanged on the other, and every wall comparison made with
`perfindex compare` is a paired A/B/A' on the host it ran on (the
instruction counts are single runs, and a few per-change figures are
in-process timings or `Ir` against the pre-change tree).  rustc 1.94.1; CI
pins clippy 1.98 for crypto (`rust-lint.yml`) and stable for the suite.  The
per-change table does not say which host produced each figure, and a figure
need not come from the host that built its baseline: the `residual_walk`
comparison ran on the 2.10 GHz host against `crypto-base-v5`, which the notes
say was built on the 2.80 GHz host.  The per-PR instruction indices carry no
host manifest; their counts do not depend on host load.

- **Shared-host load.**  Wall figures were taken on a shared 4-core
  container while other builds ran.  Load averages recorded in the
  write-ups run from 0.03 (a quiet-host review run) to 21.  Spreads from
  the A/A arm, where recorded, are wide: one review run reports ±3.6% for
  `rho_dp_zp_multi` and ±24.4% for `rho_floyd_zp_q36` in the same run, and
  the implementer's run of the same kernels at load 10–20 gives ±16.7%
  (`rho_floyd_zp_q36`) and ±18.7% (`rho_dp_zp_multi`).  A difference inside
  the A/A spread is not a result (`PERFORMANCE_INDEX.md`).
- **Isolation rule not met.**  `AGENTS.md` §10 makes isolated timed runs
  through `tools/isolated_bench.py` mandatory for any wall-clock number, says
  wall time is evidence only from uncontended runs, and asks that the number
  of contended runs be reported.  No source shows that tool in use: the
  write-ups use `perfindex compare`, some pinned with `taskset`, at the load
  averages above, and record no count of contended runs.  The per-change wall
  figures therefore do not meet §10; the deterministic instruction counts
  are the primary figures.
- **Same kernel, different runs.**  `ecc2k130_certify_m131_len12` measured
  23.66x (CI 14.23–27.06, A/A ±6.9%) in the implementation run and 17.49x
  (CI 17.23–19.56, A/A ±9.7%) in review; the notes quote the review figure.
  `subspace_oracle_n20_l7` measured 3.07x at load 11–15 and 2.815x (CI
  2.791–3.087) on a quiet host (load 0.03–0.35); the notes quote the quiet
  figure, 2.82x.  `aut_folded_rho_j0_p25` carries CI 12.3–21.6 around its
  19.5x (load 7–9).  `koblitz_rung_k0n31` measured 3.77x (CI 3.50–5.45, A/A
  ±8.9%, load about 9) in review and 3.12x (CI 2.57–5.10, A/A ±13.0%) in the
  implementer's unpinned run, 3.29x pinned (CI 2.97–3.62) and 3.22x at four
  threads; the notes quote the review figure, the implementer's summary says
  about 3.1–3.3x.  `rho_dp_zp_multi_q34_x3` and `rho_floyd_zp_q36` measured
  19.70x (CI 13.36–35.93) and 27.74x (CI 21.49–31.92) in review, and 24.46x
  (CI 14.01–29.21) and 21.86x (CI 14.53–22.96) in the implementer's run, so
  the order of the two kernels flips between runs; the notes quote the review
  figures, and the Floyd figure lies above the implementer's CI.
- **Kept on instruction counts, wall neutral.**  Some changes were kept
  because their `Ir` gain is deterministic although the wall effect is
  inside the A/A spread: the descent changes in `83191df0` (1.05x Ir, wall
  neutral), the `residual_walk` state arena (Ir 1.216x; wall 1.02–1.06x in
  five paired runs), and the `residual_walk` `mul_mod` word branch (wall
  0.999x, CI 0.936–1.034; Ir 1.032x).
- **Tests under load.**  The `koblitz_symmetrised` gate test asserts a
  wall-clock ratio and failed once under load 13.
- **No campaign-wide wall index.**  The sources contain no wall-time index
  for the campaign.  `PERFORMANCE_INDEX.md` asks that the instruction index
  be reported next to wall time, never instead of it; at campaign level only
  the instruction index exists, and wall figures exist only per change.
- **Multi-core.**  The write-ups for the seven changes that have one
  (`ecc2k130_guard`, `pollard_rho`, `residual_walk`, `koblitz_sparse_la`,
  `pq_sparse_wiedemann`, `SubspaceOracle`, `koblitz_index_calculus`) record
  four-thread runs; six state that there is no multi-core regression and the
  `koblitz_index_calculus` write-up says no kernel's verdict is "slower".
  Five of the seven say the kernel is single-threaded, so for them a
  four-thread run says nothing about parallel scaling.  Two four-thread
  readings fall below 1 and were judged neutral only by the A/A rule:
  `koblitz_rung_k0n31` in the sparse-LA write-up at 0.908x (CI 0.856–0.932;
  a rerun gave 0.983x, CI 0.980–1.218) and `koblitz_rung_k0n41` in the
  index-calculus write-up at 0.928x (CI 0.902–1.126, `Ir` identical).  The
  sources record no four-thread run for the other accepted changes.  The
  index itself is at one thread.

### What callgrind cannot see

- **AVX-512 and VPCLMULQDQ.**  Callgrind never runs these paths (Valgrind
  does not emulate them); it measures the scalar and AVX2 dispatch paths.
  For kernels that dispatch to them, wall time is the headline.
- **A divide is one instruction.**  Replacing a divide by a few
  multiplications raises `Ir` while the time falls: the C library's
  `ic_fbdiv` measured 2.57x wall for 0.82x Ir.  The inverse also happens:
  replacing software 128-bit division by hardware divides lowers `Ir` more
  than the time, and `residual_walk_lmw_b22_x3` reads 3.58x Ir against 2.06x
  wall.
- **Mispredictions and cache misses are invisible.**  The branch-free
  Pollard-rho step forms three products per step to avoid a mispredicted
  branch: on the Floyd walk it reads 1.56x wall for 2.69x more `Ir`, on the DP
  multi-target walk 1.25x wall for 1.83x more `Ir`, each against the branchy
  word step.  A branchy Shoup reduction in `koblitz_sparse_la` cut `Ir` and
  read 0.79x wall.

### Kernels below 0.99 in the campaign product

| kernel | campaign speedup | what the notes say |
|:--|--:|:--|
| `field_ec/ct_scalar_mul_p256_x16` | 0.866x | inlining artifact, verified with callgrind (profiled at `0ec2a0dc` and `717d0aaf` only, the #915 flip) |
| `field_ec/ct_scalar_mul_secp256k1_x16` | 0.936x | recorded with the p256 kernel as an inlining artifact ("same pattern"); no callgrind profile of this kernel is recorded |
| `fp_gb/f4_pkm_kummer_m2_t5` | 0.973x | 1.3% two-state flip with no relevant commit; same pattern; not profiled |

- **What was verified (p256 kernel).**  Callgrind at `0ec2a0dc` (96,443,183
  Ir) and `717d0aaf` (103,443,151 Ir), the #915 flip (0.932x) only:
  `Uint::mont_mul` is identical (79,055,648 Ir, 76-82% of the kernel).  The
  only difference is that `Uint::add_mod` and `Uint::sub_mod` are inlined
  into `P256ProjectivePoint::double`/`add` in one build and compiled as
  separate functions in the other.  The range `0ec2a0dc..717d0aaf` contains
  no commit touching `src/ecc` or `src/ct_bignum.rs`; that statement is
  scoped to this range.  The #908 flip (0.929x, `92ac774b` to `524d3877`)
  was not profiled, and its window does contain a commit touching
  `src/ecc/point.rs` (`d4851a0a`).  Of the product 0.866x = 0.929x * 0.932x,
  only the 0.932x factor is verified; the 0.929x factor rests on the
  notes' attribution of the same pattern.
- **Why the product reads below 1.**  The kernel's `Ir` is 96.1M at
  `92ac774b`, 103.4M at `524d3877`, 96.4M at `0ec2a0dc` and 103.4M from
  `717d0aaf` on: it flips with unrelated changes elsewhere in the crate.
  #908 and #915 each contain an upward flip (0.929x and 0.932x); the
  downward flip lies between them, outside both windows; the product counts
  both upward flips.  The secp256k1 kernel is also two-state (66.0M / 70.5M /
  66.1M) but not with the same pattern: its `Ir` is 65,960,846 at
  `92ac774b`, 70,466,510 at `524d3877`, `0ec2a0dc` and `717d0aaf`,
  66,141,198 at `b7d3f9ca` and 65,960,430 at `88ff556f`.  Its product counts
  one upward flip (#908, 0.936x); #915 reads 1.000x, #1167 and #1172 read
  1.003x and 0.997x, and the downward flip lies between `717d0aaf` and
  `b7d3f9ca`, outside every window.
- **Not applied.**  The notes say that an `#[inline]` hint on
  `add_mod`/`sub_mod` in `ct_bignum.rs` would pin it; no source records it
  applied or measured.  That file belongs to the constant-time ladder task,
  so the hint is left to it.  Other kernels that run through these helpers
  may carry the same few-percent layout noise.
- **`fp_gb/f4_pkm_kummer_m2_t5` is unexplained.**  Its `Ir` by revision is
  1.414G at `92ac774b`, 1.433G at `524d3877`, 1.413G at `0ec2a0dc`, 1.413G at
  `717d0aaf`, 1.413G at `b7d3f9ca` and 1.432G at `88ff556f` (per-PR raw
  data); each flip is 1.3-1.4% (+1.34% in #908, +1.39% in #1167), and no
  commit in any of the compared ranges touches the F4-over-F_p code.  It is
  consistent with the same code-layout effect, but it was not profiled, so
  treat it as unexplained noise of that size, not as a verified artifact.

### Other kernels with a "slower" verdict in the per-PR tables

The product multiplies per-PR ratios, so a loss in one PR can be hidden by a
gain in another.  Besides the `ct_scalar_mul` kernels above, these kernels
carry a "slower" verdict in a per-PR table.  The notes give no cause for any
of them.  Many more kernels read below 1.0 with a "neutral" verdict (40, 24,
58 and 30 kernels in #908, #915, #1167 and #1172, for example #1167
`fp_gb/gaudry_s4_macaulay_p1039` 0.9815x and `gf2_la/rref_macaulay_k5m3d5`
0.9822x); they are not listed.

| PR | kernel | `Ir` base → candidate | speedup |
|:--|:--|:--|--:|
| #908 | `sat/koblitz_bool_cnf_m3_n15_b1500_x2` | 1,000,382,121 → 1,024,335,738 | 0.977x |
| #908 | `sat/random3sat_n150_r426_x12` | 914,630,364 → 938,717,608 | 0.974x |
| #908 | `sat/random3sat_n200_r426_x5` | 876,744,330 → 903,010,298 | 0.971x |
| #908 | `sat/semaev_s4_n15l5_cnf_x2` | 345,339,565 → 357,782,968 | 0.965x |
| #1167 | `bool_gb/fes_moebius_random_n22` | 161,445,807 → 165,640,113 | 0.975x |
| #1167 | `dlp/rho_floyd_zp_q36` | 107,948,870 → 110,164,568 | 0.980x |
| #1167 | `gf2_la/echelon_macaulay_k5m3d5` | 217,198,053 → 222,216,685 | 0.977x |
| #1172 | `dlp/rho_dp_zp_multi_q34_x3` | 35,668,842 → 36,703,717 | 0.972x |

### Hardware class

Every measurement here is on Linux x86-64, on the two Intel Xeon hosts
above.  Arm64 (Apple silicon, Graviton), GPUs and FPGAs were not measured,
and nothing here implies them.  The write-ups say so per change:

- `subspace_oracle_n20_l7`: aarch64 was neither compiled nor measured; only
  x86_64 gets the pclmulqdq-compiled pair loop.  The x86 portable fallback
  was timed in process at 1.55-1.71x, against 2.76-2.82x with pclmulqdq.
- `ecc2k130_certify_m131_len12`: the portable (non-BMI2) path is tested
  against the BigUint search at every width but was never timed on a CPU
  without BMI2; Arm64 was not measured.
- The `koblitz_index_calculus` change (`83191df0`): Arm64 and non-clmul x86
  were not measured.
- The `pollard_rho` change (`4fc58708`): the branch-free step was measured on
  x86-64 only; Arm64 is expected, but not measured, to behave the same.

### What the index does not say

- **Not an end-to-end speedup.**  A faster kernel lowers the time per
  counted operation; it does not change the operation count `S`, a degree of
  regularity, a relation yield or any exponent (`PERFORMANCE_INDEX.md`;
  `AGENTS.md` §8).
- **Not one comparison.**  1.904x is a product of four ratios, each against
  a different base (main just before each PR).  It is not the speedup of the
  final tree over one fixed baseline; that comparison was refused as
  invalid (see "How it was measured").
- **Fixed kernels, equal weights.**  The kernels have fixed seeded inputs;
  areas weigh the same however many kernels they have.  No workload profile
  is in the sources, so no projected workload speedup (the Amdahl form in
  `PERFORMANCE_INDEX.md`) is given.

### Where the sources disagree

- The campaign instruction index (`campaign.md`) lists
  `fp_gb/f4_pkm_kummer_m2_t5` at 1.433G Ir at `717d0aaf`; the per-PR raw data
  for #915 has 1.413G at both its base and `717d0aaf`.  This file uses the
  raw data, which is also what the product 0.973x requires (flips in #908
  and #1167, none in #915).
- The notes give sat-core as "sat area 1.180x Ir"; the per-PR index for
  #1172 gives 1.178x.  This file uses the index.
- The notes record the `f9c77524` baseline harness (`crypto-base-v6`) as 153
  kernels; the per-PR comparisons cover 133.  `PERFORMANCE_INDEX.md` rule 4
  puts the slower cells in `Tier::Full`, the write-ups run
  `relation/koblitz_rung_k0n41` as `Tier::Full` (`--full`), it is not among
  the 133 in the per-PR raw data, and the per-PR runs of `perfindex.py instr`
  pass no `--full`.  So the index excludes at least the `Tier::Full` kernels,
  a scope limit; the sources do not list the other kernels outside the 133.
- The notes' row for `mod-inverse-u64` (`8925e497`) gives
  `koblitz_signed_rho` 1.116x `Ir` and `pohlig_hellman_smooth47` 1.173x `Ir`,
  "others neutral".  The per-PR raw data disagree.  For #1167,
  `koblitz_signed_rho_k0_n41` reads 1.000x (and 1.001x, 1.000x, 1.000x in
  the other three PRs, so no per-PR window holds a 1.116x),
  `pohlig_hellman_smooth47` 203.18x, and
  `gaudry_schost_negation_demo32` and `collab_walkers_demo32_x8` 17.78x and
  14.40x.  The 1.173x is exactly the #908 per-PR factor for
  `pohlig_hellman_smooth47` (410,253,853 to 349,700,822 `Ir`), and the
  `ecc2k130_guard` review attributes the 1.17x on `pohlig_hellman` (and 1.05x
  on `gaudry_schost`) to changes in `src/ecc/point.rs` and similar files
  already between `f9c77524` and `b072fcf5`, the parent of this change's
  commit, not to the `ecc2k130_guard` diff.  Git shows that commit
  (`fb72b23b`) touching only `aut_folded_rho.rs`, `cga_hnc.rs` and
  `utils/mod.rs`.  The notes do not state the base of the
  1.116x and 1.173x figures; this file does not credit them to `8925e497`.
  The #1167 steps of the three dlp kernels come predominantly from the
  `ec-vartime-walks` change, whose own measurements
  (`results-ec-vartime-walks.txt`, first round, against the tree before the
  change) give 17.70x, 14.31x and 201.1x `Ir` for them; the notes'
  accepted-changes ledger did not list that change with figures.
- The notes call the 2.80 GHz host (H2) the host of "most batch A/B
  measurements".  Of the per-change write-ups that name a CPU, four
  (`ecc2k130_guard`, `pollard_rho`, `residual_walk`, `pq_sparse_wiedemann`)
  report the 2.10 GHz host (H1) and one (`koblitz_index_calculus`) the 2.80
  GHz host; the `koblitz_sparse_la` and `SubspaceOracle` write-ups name no
  host.

## Tried and rejected

Listed in the notes:

| what | numbers | note |
|:--|:--|:--|
| LTO fat + `codegen-units=1` (2026-09-26, bool_gb + gf2_la kernels, 3 paired rounds) | `gf2_la/rref_random_2048` 0.77x; area indices 1.005x / 0.769x; no kernel reliably faster; build 3 min | default profile stays (`PERFORMANCE_INDEX.md`, "Build settings tried") |
| `f4_gf2` slices | no figure in the notes | |
| `FastField` Euclid inversion | slower than Fermat on toy40; no figure in the notes | |
| orbit `add_r`/`sub_r` | no figure in the notes | suite |
| relation solver: Möller–Granlund 2-by-1 | 0.82x | |
| relation solver: Shoup with branches | 1.65x slower | |
| relation solver: single 128-bit path | 17% slower | |
| sparse LA: `chunks(K)` split | block_wiedemann `Ir` 728.7M vs 713.3M for the hand split | |
| sparse LA: per-product fold counter | row kernel 266M `Ir` vs 222M (spilled an accumulator) | |
| sparse LA: branchy Shoup | 0.79x wall despite lower `Ir` (0.789x block_wiedemann, CI 0.640–0.871; 0.851x sparse_solve, CI 0.775–0.891); 4.36M mispredicts vs 0.63M | branch-free folds fixed it: 1.458x / 1.389x against the lazy-rows-only build |

Further variants measured in the per-change write-ups.  The `ecc2k130_guard`
variants were measured in a standalone copy of the search (m=131, length 12),
as search-only `Ir` and wall time, the minimum over 15 interleaved rounds, so
they are not the kernel's `Ir`.  The `pollard_rho` variants were measured in a
standalone micro-benchmark of the step on the Floyd q36 instance, best of 80,
under load.  The first two `SubspaceOracle` rows are `Ir` of
`pdp/subspace_oracle_n20_l7` against the inherited first stage (203.2M;
690.0M before the change, 192.3M after).

| change | variant | numbers |
|:--|:--|:--|
| `koblitz_sparse_la` | `std::hint::select_unpredictable` | same as `min()` (`Ir` 476.3M vs 479.5M, same 31.6k mispredicts), but needs Rust 1.88 while the suite declares 1.87; `min()` kept |
| `residual_walk` | u32 divides in the Euclid when p < 2^32 | u64/u32 time ratio 0.94–1.04, no gain |
| `residual_walk` | quotient-one shortcut in the Euclid | 0.68–0.74x |
| `residual_walk` | mask-form corrections instead of `select_unpredictable` | 0.678x wall (CI 0.627–0.692); `Ir` 238.5M vs 234.8M |
| `ecc2k130_guard` | all sibling products before recursing | `Ir` 32.4M vs 29.6M; 2.60 ms vs 2.29 ms |
| `ecc2k130_guard` | product-scanning (FIPS) Montgomery | `Ir` 39.6M vs 29.6M |
| `ecc2k130_guard` | one loop with a loop-invariant `deeper` flag | `Ir` 30.3M vs 29.6M; 2.30 ms vs 2.26 ms |
| `ecc2k130_guard` | adx beside bmi2 | identical code; 25.02M `Ir` both ways |
| `pollard_rho` | exponents through one selected addend (v3) | 3–8 fewer `Ir` per step, no wall gain (19.4–27.6 ns/iter against 17.7–18.9 for the committed step) |
| `pollard_rho` | v3 plus a 256-byte partition table (v4) | 2–6 fewer `Ir` again, no wall gain (8.5 ns/step) |
| `pollard_rho` | select the multiplier first, two products instead of three (v5) | 24.3–32.8 ns/iter, 1.3–1.9x slower |
| `SubspaceOracle` | `assert_unchecked(has_clmul)` | `Ir` 203.2M → 217.7M (+7.2%) |
| `SubspaceOracle` | per-target tables of v²x² and bv² | `Ir` 203.2M → 228.0M (0.891x); wall 1.013x (CI 0.851–1.068, A/A ±5.4%) |
| `SubspaceOracle` | `Gf2::sqr` rewritten as `reduce(sqr_wide(a))` | against the tip, `Ir` `pair_table_build` 0.930x, `mitm_m3` 0.933x, `mitm_m4` 0.954x; reverted, replaced by a private fused squaring used in the quartic route only |
| `m_can_decompose` | `FastCurve` walk with one `add` per pair, no batching | 295.7M `Ir` vs 244.4M on `koblitz_rung_k0n31`; superseded by `add_many` rows |

## Decisions and open choices

### Coordinator decisions

**`ec-vartime-walks`: composite p (review 1, 2026-09-29).**  The Euclid
`inv_vartime` differs from the Fermat `inv` only off the prime domain.  On
composite-p `pollard_collab` jobs, which lie outside `FieldElement`'s
documented prime-only contract (`JobContext::new` never checks primality),
the branch tables and the walker DPs change; on every prime-p job they are
bit-identical (reviewer's 338-solve / 185-DP transcripts).  The change was
accepted as documented, not as a refusal: rejecting non-prime p in
`JobContext::new` would be a behaviour change and belongs in its own PR, and
the reviewer queued it as a follow-up task.  The docs must say "over a prime
p" wherever parity is promised.

**Review 2.**  The reviewer re-raised composite p and offered a
bit-identity-preserving alternative: one primality test per context, with
the Fermat path kept for composite p.  The sign-off was kept: composite p is
outside the documented prime-only contract, the output on it was not a group
computation before either, and the gate would add a primality test to
`JobContext::new`, `EcGroup::new`, `pohlig_hellman_curve` and
`cheon_attack`.  Review 2's other items were fixed by the coordinator (`f1685a91`):
even moduli are stepped affinely in `linear_dlog_vartime` (F_2 panic parity;
the reviewer's exhaustive F_2 test is un-ignored), the `scalar_mul_vartime`
doc excludes p = 1, and `JobContext::half_p` is documented as derived.

### Open choice and open notes

Only the composite-p question is recorded in the notes as an open choice for
the user.  The other items are the notes' hand-offs and open notes.

- **Composite p in the cryptanalysis walkers (open choice).**  Three positions are on
  record: leave the behaviour as merged, with the prime-only contract
  documented (the coordinator's current position); reject non-prime p in
  `JobContext::new` in a separate PR (the reviewer's queued follow-up, a
  behaviour change); or add one primality test per context and keep the
  Fermat path for composite p (the reviewer's alternative, which preserves
  bit-identity).  The notes list it as an open choice for the user and
  record no later decision.
- **`ct_scalar_mul` inlining (hand-off).**  The notes hand an `#[inline]`
  hint on `add_mod`/`sub_mod` in `ct_bignum.rs`, which they say would pin the
  artifact, to the constant-time ladder task (`ecc-ct-ladders`); no source
  records it applied or measured.  Batch D, which holds that task, was
  resumed on 2026-10-02 (18:05 UTC) and the task's notes were extended with
  this finding.
- **`fp_gb/f4_pkm_kummer_m2_t5` (not profiled).**  The 0.973x reading is
  unprofiled.
- **`dp_walks_word` memory (open note).**  It allocates 520 B per target up
  front, where a target cost about 72 B before; at m = 10^6 that is about
  520 MB against about 72 MB.  The only caller in the tree uses m = 3.
- **`koblitz_symmetrised` gate test (note).**  It asserts a wall-clock ratio
  and failed once under load 13 (pre-existing, unrelated path).
- **Not merged, not counted.**  As of the notes' last entries (2026-10-02,
  18:05 UTC): batch C (`gf2la-misc`, `fp-tower`, `f4-fp`, `gaudry-mm`,
  `cantor-fp-inv`, `binary-ecc`, `gf2-fold-reduce`) was resumed; batch D
  (nine tasks, among them `kg-echelon` and `ecc-ct-ladders`), held back at
  16:10 UTC, was resumed at 18:05 UTC with concurrency 1; suite ports to the
  cryptanalysis branch have no PR yet (`port-f2-sorts`, `fa9b807`, is merged
  on that branch).  The notes expect `fp-tower` to conflict with main's
  `f4_fp_tower` (changed by `00983dfc` and `a8b1e7ad`) and to need a
  fingerprint re-pin.

## Reproducing

The commands of `PERFORMANCE_INDEX.md`, for building both sides with the
same harness and comparing:

```bash
# 1. Build both sides with the SAME harness (this tree's examples/perfbench).
python3 scripts/perf/perfindex.py build --ref <base-rev> --out /tmp/pi-base
python3 scripts/perf/perfindex.py build --current      --out /tmp/pi-cand
# 2. Wall time, paired and interleaved, one thread (the headline).
python3 scripts/perf/perfindex.py compare --base /tmp/pi-base/perfbench \
    --cand /tmp/pi-cand/perfbench --rounds 5 --threads 1 --out /tmp/pi-1t
# 3. The same at every core, and the deterministic instruction count.
python3 scripts/perf/perfindex.py compare ... --threads 4 --out /tmp/pi-4t
python3 scripts/perf/perfindex.py instr --base ... --cand ... --out /tmp/pi-ir
```

The per-PR indices above build the first parent and the merge commit
themselves with `--ref`, then run `instr`.  For PR `<pr>` with first parent
`<base>` and merge commit `<merge>` from the table below:

```bash
python3 scripts/perf/perfindex.py build --ref <base>  --out /tmp/pi-<pr>-base
python3 scripts/perf/perfindex.py build --ref <merge> --out /tmp/pi-<pr>-cand
python3 scripts/perf/perfindex.py instr --base /tmp/pi-<pr>-base/perfbench \
    --cand /tmp/pi-<pr>-cand/perfbench --out /tmp/pi-<pr>-ir
```

| PR | base (merge commit's first parent) | candidate (merge commit) |
|:--|:--|:--|
| #908 | `92ac774b` | `524d3877` |
| #915 | `0ec2a0dc` | `717d0aaf` |
| #1167 | `b7d3f9ca` | `88ff556f` |
| #1172 | `848f145b` | `441293ee` |

`git rev-parse <merge>^1` confirms each base.  Use one harness for every
build: the index is only comparable when both sides run the same kernels,
and the per-PR figures here are for a harness with 133 kernels.  A
comparison whose fingerprints differ is invalid and reports no index.  The
campaign figure is the per-kernel product of the four speedups, with the
area indices and the index then computed by the same formula.

To append a later PR: add its row to the PR table and its column to the area
table, add its changes to "Accepted changes" and its pair to the table
above, and recompute the campaign figure as the per-kernel product over all
PRs.
