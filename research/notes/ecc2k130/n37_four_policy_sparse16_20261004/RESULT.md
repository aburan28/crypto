# Sparse n37 four-policy PDP: exact misses, no fixed-block seed-selection lead

**Decision: `NO_FIXED_BLOCK_M3_SELECTION_LEAD`.** The K42 four-policy panel
decomposed all 2,048 of its targets, so the preregistered follow-up took the
first K16 signed-Frobenius columns of each already frozen base before creating
new Q. Its formal three-summand multiset capacity is only 1.203 times the
prime subgroup order, and the exact search did find misses. On two fresh
point-only 1,024-target blocks, the original source base hit 1,435/2,048 and
the descendant-native base hit 1,445/2,048. Their difference is **10/2,048
= 0.488 percentage points**, below the frozen 21/2,048 gate, with exact
paired McNemar **p≈0.754762**. Both policies reached rank 16 on the fixed
target-blind probe stream, and every reported logarithm was independently
verified. This is an exact PDP/rank diagnostic, **not** a cold IC/rho result,
a complete recovery of all 2,048 Q, or evidence of an ECC2K-130 speedup.

The [protocol](PROTOCOL.md) was committed as `469126b2` and [PR #1326](https://github.com/aburan28/crypto/pull/1326)
opened before any new target was generated. The native point freeze and
separate general-law replay were committed as `746166b0` before the PDP
source was run. That replay regenerated all 2,048 scalar-to-point mappings,
every accepted/rejected candidate, and all signed-Frobenius orbit keys; it
confirmed 4,096 distinct orbits across the old b03/b04 and new blocks. The
producer/replay source was committed as `7039ab0b` before measurement. The
producer reads only point files, never the verifier-only scalar labels.

Each of the four policies has **16 effective log columns, 592 signed classes
and 1,184 distinct usable physical subgroup points**. Original/transported
and native/pullback are algebraically paired controls through the verified
degree-73 map; only source-selected versus native-selected support is a
different base. The complete m≤3 oracle enumerated 702,705 raw identity,
singleton and nondecreasing-pair candidates per policy, then checked every
residual up to a proved hit or all 1,185 positions. The secondary m≤2 oracle
used the same fixed support.

| Frozen policy | Distinct table sums | m≤2 hits / 2,048 | m≤3 hits / 2,048 | Probes to rank 16 | Proved rank misses | Verified target logs / 2,048 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Original source | 697,673 | 9 | 1,435 | 24 | 5 | 1,435 |
| Transported leaf | 697,673 | 9 | 1,435 | 24 | 5 | 1,435 |
| Descendant-native leaf | 696,415 | 8 | 1,445 | 23 | 7 | 1,445 |
| Pullback source | 696,415 | 8 | 1,445 | 23 | 7 | 1,445 |

The m≤3 hit rates are 70.068% source and 70.557% native. Descriptive
95% Wilson intervals for the deterministic hash-sampled blocks are
[68.049%, 72.013%] and [68.546%, 72.491%]. There are **410 source-only** and
**420 native-only** hits. A paired normal-approximation 95% interval for the
native-minus-source hit-rate difference is approximately
[-2.269, +3.245] percentage points; the exact conditional test, not that
approximation, applies the frozen decision rule. The m≤2 control gives
9 source-only and 8 native-only hits, also no selection lead. Across blocks,
source m≤3 hits were 733/1,024 and 702/1,024; native hits were 704/1,024
and 741/1,024. That reversal is another reason not to choose a base from
one block's apparent rank or yield advantage.

The producer built each table with 701,520 full-point additions and retained
at least 11,262,240 bytes in table vectors. The four-policy run's macOS L0
peak RSS was 66,289,664 bytes. Its 7.10-second descriptive wall interval
includes all four tables, both m≤2 and m≤3 target controls, rank, checks and
JSON output. Source and native rank probes required 8,806 and 11,760 table
lookups, respectively; their m≤3 target queries required 1,046,053 and
1,037,911 lookups. These are stage diagnostics. The archived support excludes
fresh source-base scan, leaf-seed selection, scalar orbit closure, transport
and pullback costs; there is no same-Q rho arm or isolated host. Consequently
complete `S`, cold and online IC/rho ratios, and n131 transfer remain null.

The independent Rust replay used a hash-map table instead of the producer's
sorted packed table and the general binary-curve group law for source
witnesses and scalars. It checked both coupled leaf arms via the degree-73
map, both exact tables, all four rank trajectories, all 8,192 policy-target
decisions and every one of the 5,760 reported target logarithms. The raw
4,669,993-byte manifest is retained losslessly as [`RESULT.json.gz`](RESULT.json.gz):
compressed SHA-256 `dbafb183d4edf27813dd15bc22aed6e1fe9809887aa3cd54bedf589120d43357`,
expanded SHA-256 `a17e6b2cad2f33e6f9df950bed272eb6dcc85bd13594625527250e61e7c0e1f7`.
The [replay receipt](REPLAY_INITIAL.json) has SHA-256
`8b4b4c184cee5541f233a706d432c3f733c30b5562a9feae6e821f691142a05e`.
Paired witness, rank-row and hit-flag mutations, a changed decision, and a
changed target file all failed replay as intended. [EVIDENCE.json](EVIDENCE.json)
pins source, input, binary, result and mutation receipt hashes and the
non-isolated host boundary. The measured [Cargo.lock](Cargo.lock) is archived
because the repository root lockfile is ignored by Git. Linux CI independently rebuilds, replays and
reproduces the canonical non-timing result before this PR can merge.

This K16 result resolves the K42 saturation without finding a stable
source-versus-descendant selection advantage. It does not test a degree-263
descendant or natural high-arity PDP at n131. The next factor-base yield gate
should move to a larger subgroup with equal useful support and a native,
non-eager oracle; a cold single-target comparison must charge base and map
construction and pair each completed target with automorphism-aware rho.
Repeating this exact K16 eager table as a speed claim would conflate its
stage coverage with the required end-to-end attack cost.

Reproduce in a disposable checkout from the repository root with the
committed source lock. `SOURCE_Cargo.toml` preserves the exact manifest from
source commit `7039ab0b2940ac55809da62c91361b8082625716`, SHA-256
`f1b3abc246c4971b48b377d5696bc96594deb729320348c4277307e42ed61721`.
The current root manifest gained later binaries, so the historical replay
restores its measured manifest and lockfile before checking the unchanged
source lock. The old measurement and result remain unchanged.

```sh
cp research/notes/ecc2k130/n37_four_policy_sparse16_20261004/SOURCE_Cargo.toml Cargo.toml
cp research/notes/ecc2k130/n37_four_policy_sparse16_20261004/Cargo.lock Cargo.lock
sha256sum --check research/notes/ecc2k130/n37_four_policy_sparse16_20261004/SOURCE_LOCK.sha256
cargo build --release --locked --example n37_four_policy_sparse16_inputs_replay --example n37_four_policy_sparse16_pdp_replay
target/release/examples/n37_four_policy_sparse16_inputs_replay research/notes/ecc2k130/n37_four_policy_sparse16_20261004/inputs /tmp/n37-sparse16-input-replay.json
gzip -dc research/notes/ecc2k130/n37_four_policy_sparse16_20261004/RESULT.json.gz > /tmp/n37-sparse16-result.json
target/release/examples/n37_four_policy_sparse16_pdp_replay /tmp/n37-sparse16-result.json /tmp/n37-sparse16-result-replay.json
cmp /tmp/n37-sparse16-result-replay.json research/notes/ecc2k130/n37_four_policy_sparse16_20261004/REPLAY_INITIAL.json
```
