# Two-word F6 geometric closure on the n=83 K_1 field

The [frozen stage protocol](PROTOCOL.md) ran once in a release build on an
arm64 macOS 26.6 host with Rust 1.93.1. The tested curve is K_1 over the
registered n=83 field, subgroup order 8,569,786,107,849,059. It is distinct
from the K_0 n=83 confidence-gate curve with its 81-bit subgroup.

The exact pair index and batched affine adder passed the n=83 reference-law
tests, including exceptional additions. The probe verified all pair sums
against scalar point additions, found a planted four-summand witness, and
proved the selected outside-range target absent from this small supplied
point list. The point list consists of consecutive generator multiples. It
is an arithmetic control, **not an admissible natural-yield factor base**.

An untimed follow-on correctness test uses the pinned K_1 curve identity
`icv1-f2m83-t6151469093347-cdcc5432` and a genuine polynomial-subspace
base of dimension five. It constructs 37 curve points, projects them by the
public cofactor to 36 distinct nonidentity subgroup-usable points (18 signed
columns), verifies each point's subgroup order, and independently replays
an exact planted four-point closure. This tiny base is a correctness fixture;
its planted hit says nothing about ordinary target coverage.

| Points | Unordered pairs | Pair build reference / batch median ms | Pair ratio | No-witness query scalar / batch median ms | Query ratio |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 136 | 0.544 / 0.181 | 3.00× | 0.572 / 0.078 | 7.36× |
| 32 | 528 | 2.142 / 0.513 | 4.17× | 2.218 / 0.275 | 8.05× |
| 64 | 2,080 | 8.619 / 1.434 | 6.01× | 8.817 / 1.046 | 8.43× |

All three paired timings for each arm are in [raw.jsonl](raw.jsonl).
The measured implementation at commit `32b9601fe` was
`f6_wide_geometry.rs` SHA-256
`e6c2bf0a2f91a1444efcd1a7c7d1010c6edf27d6493d88dc3f6b93f1d63efb32`;
the probe was `f6_wide_n83_probe.rs` SHA-256
`5a6a604616639ddfbce312bdeca53566bbedb0455c1bc0d54a7827ea10274795`;
the protocol was SHA-256
`d3791cdaf253cb0b10c1a48a78bcfcd5e6ef4df9f1b6ce1d257c5c2bc165663b`.

These are exploratory component timings on an unisolated host. They do not
measure a complete ordinary-query PDP, F6 versus F4/F5, factor-base
construction, relation yield, final LA, target recovery, or rho. The current
F6 Boolean backend still cannot encode the four-summand n=83 chain, and the
pair closure is a known meet-in-the-middle component. No n=83 algorithmic or
end-to-end IC speedup is established by this result.

## Exact coverage ceiling for a four-summand tail

For any base of `B` distinct usable subgroup points, at most
`binomial(B+3,4)` unordered four-point multisets exist. Each has only one
group sum. Thus a uniformly drawn subgroup target has four-summand coverage
at most `min(1, binomial(B+3,4)/r)`, with no independence or distribution
assumption. Repeated sums only lower the coverage. On the tested K_1 subgroup,
`B=4,096` gives at most **0.1371%** coverage and requires 8,390,656 pair
entries; reaching a 1% ceiling needs at least `B=6,733` and 22,670,011
pairs. On the K_0 n=83 confidence-gate subgroup (`r =
2,417,851,639,230,796,216,685,689`), the same `B=4,096` ceiling is
`4.86e-12`; even a 1% ceiling needs `B=872,790`, or 380,881,628,445
pair entries. These are ceilings, not predicted yields. The wide pair path
is useful as an exact terminal closure and arithmetic kernel, but a complete
n=83 F6 design needs a higher-arity strategy or a different relation law.

## Actual K_1 base inventory

The separately [registered inventory](BASE_PROTOCOL.md) constructed
standard polynomial-subspace bases and then applied the public cofactor.
These are actual point counts, not nominal `2^ell` sizes. No target solve
or natural-yield measurement occurred.

| Dimension | Curve points | Usable projected points `B` | Signed columns | Exact four-summand coverage ceiling | Pair entries |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 251 | 250 | 125 | 1.945×10⁻⁸ | 31,375 |
| 10 | 997 | 996 | 498 | 4.814×10⁻⁶ | 496,506 |
| 12 | 4,135 | 4,134 | 2,067 | 0.1422% | 8,547,045 |

The [raw inventory](base_inventory.jsonl) fixes every count. At the largest
tested base, a uniformly drawn K_1 target has four-summand success
probability at most 0.1422%, so its expected independent query count is at
least 703 if the draw law is uniform. This is a bound, not a measured rate.
The K_0 gate has a far larger subgroup; this K_1 inventory does not transfer
its factor-base counts or timings to K_0. The inventory used
`koblitz_index_calculus.rs` SHA-256
`a3b184cfd2008a74668598b3003eb643db2e9a3a8492ebfbbf996f696ae731bb`,
`f6_wide_n83_base_inventory.rs` SHA-256
`2ceb4ec04708782a7f17a2432acaf066c334bfe9398f0ac00c219c4a2be44528`,
and `BASE_PROTOCOL.md` SHA-256
`7a929d602fa4b75dd4cd0ed440ff3ae8aef5b9b5dc5c123885f73cdd7f7b904e`.

## Actual K_0 confidence-gate base inventory

The [separate K_0 protocol](K0_BASE_PROTOCOL.md) used the public generator
in `gate-m83-T001.json`. Its 81-bit subgroup order and K_1's 53-bit order
both passed exact Sage `is_prime(proof=True)` through the repository's
checked launcher; see [script](prime_check.sage), [output](prime-check.txt),
and [runtime receipt](sage-runtime-info.json). The Rust constructor also
checked the point, generator order, registry curve ID, and Frobenius
eigenvalue. This is the confidence-gate **curve**, distinct from K_1 above.

| Dimension | Curve points | Usable projected points `B` | Signed columns | Four-summand coverage ceiling | Pair entries |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 261 | 258 | 129 | 7.814×10⁻¹⁷ | 33,411 |
| 10 | 1,051 | 1,048 | 524 | 2.091×10⁻¹⁴ | 549,676 |
| 12 | 4,057 | 4,054 | 2,027 | 4.662×10⁻¹² | 8,219,485 |

The [raw K_0 inventory](k0_base_inventory.jsonl) preserves the counts.
At the largest actual base, an independently uniform target requires at
least **2.145×10¹¹ queries in expectation** for a four-summand hit even
under the most optimistic collision-free mapping. This is a mathematical
lower bound from the coverage ceiling, not a measured run time. At `B=4,054`,
the same ceiling is 0.1484% for seven summands and 75.36% for eight;
the nine-summand count exceeds the group size and gives no nontrivial
coverage guarantee. Actual yield can be lower at every arity. The direct
Semaev chain would need 594 Boolean variables at eight summands or 689 at
nine (`m·12+(m−2)·83`), far beyond the 128-variable backend. A faster
four-summand tail cannot by itself solve the n=83 gate.

The K_0 inventory used `koblitz_index_calculus.rs` SHA-256
`a183fad6ea71841c2dcba6a1f05885101e21590f5dff55664e1bda5eda1f0746`,
`f6_wide_n83_k0_base_inventory.rs` SHA-256
`4471f57411a03c54b6db7d02bd0a1a7befce31fa8740ee2183f77aecba9049d1`,
and `K0_BASE_PROTOCOL.md` SHA-256
`58f03fa4634d37e65ef15224f2a91e91508ff9f0c5d5a89fde90e5ff12480cd0`.

## K_0 paired arithmetic on T001

The [separately frozen K_0 stage protocol](K0_STAGE_PROTOCOL.md) used the
first 64 subgroup-usable points of the actual dimension-8 base and the
public T001 point. The target passed on-curve and subgroup checks. All 2,080
pair sums matched reference addition exactly. Both query methods completed
and agreed that T001 has no four-summand expression in this fixed 64-point
subset. This says nothing about the full 258-point dimension-8 base or
larger bases.

| K_0 component | Scalar median | Batched median | Median ratio | Three paired ratios |
| --- | ---: | ---: | ---: | --- |
| 2,080 pair sums | 8.917 ms | 1.312 ms | 6.80× | 6.56×, 6.95×, 4.02× |
| Exact no-witness query | 8.475 ms | 1.058 ms | 8.01× | 8.01×, 8.42×, 8.44× |

The third batched pair run was slower than the other two; all observations
are preserved in [k0_stage_raw.jsonl](k0_stage_raw.jsonl). This host has no
isolation receipt. The ratios are exploratory component diagnostics and do
not imply faster complete F6, F4/F5, index calculus, or rho on K_0. A full
four-summand query against the actual 4,054-point dimension-12 base would
have 8,219,485 pairs, and its optimistic coverage ceiling remains
4.662×10⁻¹². The stage probe used `f6_wide_geometry.rs` SHA-256
`7d5a9aece64f31fab027612ca004e310fa93f5797cca1ce7e6c1655600f7b762`,
`koblitz_index_calculus.rs` SHA-256
`1c98b492fafe80ae3a9dbea5400e7fceb114dec7b9a548e4d98f33a691f05b5d`,
`f6_wide_n83_k0_stage_probe.rs` SHA-256
`6306f4453a59f09d461123c88703086deca75eff3e3c2e37cde537dfc3ef26ef`,
and `K0_STAGE_PROTOCOL.md` SHA-256
`4287d2dc6b4e6bdd5fe308eccba8a7a1f40b0e60a02fb82c901ebb2a47b1c658`.

## K_0 pair-index representation

The [index protocol](INDEX_PROTOCOL.md) compared the original vector of
pairs per distinct sum with one representative pair per sum. For exact
four-summand search with repetitions allowed, either pair represents the
same group element. The probe built the complete pair index for each actual
base and queried public T001. All baseline and candidate runs passed their
point and subgroup checks, returned the same exact no-witness result, and
had the same pair counts. Three observations per arm, after one warmup, are
preserved in [baseline](index_baseline.jsonl) and
[candidate](index_candidate.jsonl) raw files.

| Usable points | Pairs | Index build baseline / candidate median | Build ratio | No-witness query baseline / candidate median | Query ratio |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 258 | 33,411 | 20.634 / 19.627 ms | 1.05× | 16.077 / 15.780 ms | 1.02× |
| 1,048 | 549,676 | 358.010 / 300.099 ms | 1.19× | 277.368 / 279.922 ms | 0.99× |

The build medians satisfy the registered no-more-than-10%-regression gate,
so the simpler exact representation is retained. Query time is effectively
unchanged at this sample size. The baseline was commit `7ea5f6470`, with
`f6_wide_geometry.rs` SHA-256
`7d5a9aece64f31fab027612ca004e310fa93f5797cca1ce7e6c1655600f7b762`.
The candidate file's SHA-256 is
`a8a88887e925b0416a0995a07c0e24c18d24807f80ee6bcd415418eb8b92c2f7`;
the probe SHA-256 is
`263a0b043c23f890a695370e66aaebbcd8ecd6a053c9413ed821201f5a2c5f39`;
the protocol SHA-256 is
`e7a5f4e2a7395830e7d56a7e98dbc802489a4374ee9cb4f4f6e4b72c303fed5d`.
Both arms used Rust 1.93.1 on the same unisolated arm64 macOS host. These
are exploratory component timings, not controlled CPU speedup claims or a
complete n=83 F6, F4/F5, or one-target IC comparison.

## K_0 preallocated pair lookup

The [capacity protocol](INDEX_CAPACITY_PROTOCOL.md) tested reserving hash
lookup space for the known maximum pair count. All outcomes and pair counts
matched the one-representative baseline. The [baseline](capacity_baseline.jsonl)
and [candidate](capacity_candidate.jsonl) files retain all three builds and
queries per base after warmup.

| Usable points | Build baseline / candidate median | Build ratio | Query baseline / candidate median | Query ratio |
| ---: | ---: | ---: | ---: | ---: |
| 258 | 19.996 / 18.385 ms | 1.09× | 16.194 / 16.026 ms | 1.01× |
| 1,048 | 306.919 / 294.729 ms | 1.04× | 282.504 / 278.408 ms | 1.01× |

The preallocation meets the registered retention gate and is kept. It
allocates the lookup capacity at construction start; actual peak memory was
not measured, so the memory tradeoff remains unknown. The baseline source
SHA-256 was
`a8a88887e925b0416a0995a07c0e24c18d24807f80ee6bcd415418eb8b92c2f7`,
the candidate source SHA-256 is
`7cae9b6c9610bfe17cb0c16e053b616971fdff3a91748e1fc6a41c12b2499a60`,
and the capacity protocol SHA-256 is
`29a52195325552e6e2d8c547e8ae37eaea5e568c5d78ffd41e0f298e8d7d4b00`.
The probe and compiler were the same as above on the same unisolated host;
these are exploratory stage timings only.

## Full K_0 dimension-12 pair-index feasibility

The [full-index protocol](FULL_INDEX_PROTOCOL.md) used all 4,054 actual
subgroup-usable points (2,027 signed columns), indexed all 8,219,485
unordered pairs, and queried public T001. The target passed on-curve and
subgroup checks. The first run built the index in **5.099 s** and completed
the exact no-witness query in **4.501 s**. The program wrote a PASS row in
[the first raw record](full_index_raw.jsonl), but the `/usr/bin/time -l`
wrapper exited 1 after the program because this sandbox denied its
`sysctl kern.clockrate` call; [its stderr](full_index_memory.txt) is
preserved. That wrapper produced no peak-memory measurement. The first
program's source was reconstructed byte-for-byte from the pre-instrumentation
patch as [a source snapshot](full_index_initial_source.rs), SHA-256
`efd09c56620ad868ff355db9f602d87cbdc64f4c25d24b7261973198002e0208`;
it differs from the rerun only in the added peak-RSS reading and JSON field.

The separately [registered memory rerun](FULL_INDEX_MEMORY_PROTOCOL.md)
read `getrusage(RUSAGE_SELF).ru_maxrss` inside the same process after the
query. [Its raw record](full_index_memory_raw.jsonl) reports build **5.099 s**,
exact no-witness query **4.491 s**, and peak process resident memory
**2,678,554,624 bytes (2.49 GiB)**. The rerun exited 0. Its source SHA-256
is `d17d2e1952d101b435295728248be061de04c4bd305db66f52dc6c890e8648ad`;
the retained index source SHA-256 is
`7cae9b6c9610bfe17cb0c16e053b616971fdff3a91748e1fc6a41c12b2499a60`.
The full-index and memory-rerun protocol SHA-256 values are
`402dc7450091edc1c2485702ece2b3fb62e49e7365bccc405b120633063a10f5`
and `3a81f25df1e1b412616d7e307a47c559702e5ee19f93f438a25bce6c087d1985`.

These two runs establish that the exact four-summand index can complete at
the largest tested K_0 base on this host. They do not estimate ordinary
relation yield from one target, recover a logarithm, or establish any
F6/F4/F5 or IC/rho speedup. The mathematical four-summand coverage ceiling
for this base remains `4.662e-12` for a uniform subgroup target, so
full-base feasibility alone does not meet the high-arity admission gate.
