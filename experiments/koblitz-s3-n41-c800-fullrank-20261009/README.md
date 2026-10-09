# N41 C800 full-rank S3 holdout

The 800-column low-x factor base solves its complete log system during
reusable preparation. A verified four-point S3 relation then recovers one
public target scalar from signed Frobenius labels, with the 41 inverse
target conjugates charged inside the online interval. The T025 correctness
control reached rank 800, verified the scalar, used direct recovery, and
made 26 target S3 calls. T025 was used to select this configuration and is
excluded from the timed holdout.

The [freeze](FREEZE.json) fixes Koblitz curve
`EC1N41Ckb1hcbf4a37f6e66` over `GF(2^41)` with modulus `x^41+x^3+1`,
cofactor 4, subgroup order 549756390943, and 24 new public targets
T033–T056. The [C800 candidate](identities/candidates/IC1N41Ckb1fb65600PDP4rootRCsampleLAgaussTDdirectISO0hb8875aea0063.json)
has 65,600 subgroup-usable base points, 800 signed-Frobenius columns, a
134,217,728-slot S3 table, and 839,680,000 retained indexed-coordinate
bytes. The baseline has 20,008 usable points, 244 columns, a 16,777,216-slot
table, and a rank-243 preparation threshold. Both bases use the same
natural low-x construction, with exact [C244](local/n41_base_c244.json.gz)
and [C800](local/n41_base_c800.json.gz) encodings retained for Sage replay.

The [runner](run_paired.py) fixes two repetitions per arm on each point in
ABBA/BAAB order, relation seed 20260928, 14 target-independent relation
workers, one target worker, and a 120-second cap per run. Each process
receives one public point. The online interval starts with target-dependent
work after reusable preparation and ends after scalar recovery and check;
its exclusive phases are target query, PDP, relation check, descent, and
recovery check. The [analyzer](analyze_paired.py) retains failures and A/A
repeat spread, and the [checked-Sage replay](replay_sage.py) verifies the
encoded bases, group relations, matrix rank, and recovered scalars.

The target labels were frozen after the C400 pilot and before C800 timing.
The factor base size and full-rank recovery are declared differences from
the baseline. Local wall-time ratios remain exploratory unless the CPU
host-wide isolation and noise gates are met.

## Frozen 24-target result

The [raw ledger](local/raw_holdout96.jsonl.gz), [analysis](local/analysis.json),
[checked Sage runtime](local/sage_runtime_info.json), and
[independent Sage replay](local/sage_replay.json) retain the complete panel.
All 96 native solves verified, and Sage independently recovered all 96
scalars after checking both encoded factor bases, 1,736 distinct relation
witnesses, and the two prepared relation matrices. C800 reached rank 800 on
every run and used direct target recovery with zero fallback attempts.

| One previously unseen N41 target at a time | C244 baseline | C800 full rank |
| --- | ---: | ---: |
| Verified native and Sage runs | 48/48 | 48/48 |
| Usable factor-base points / folded columns | 20,008 / 244 | 65,600 / 800 |
| Prepared rank | 243 | 800 |
| Median online wall time | 0.536 ms | 0.150 ms |
| Median target S3 calls | 1,234 | 71 |
| Median reusable preparation, supplementary | 0.259 s | 2.486 s |
| Root-table slot payload | 268,435,456 bytes | 2,147,483,648 bytes |
| Additional retained indexed coordinates | 0 | 839,680,000 bytes |

Across 24 paired target blocks, the exploratory geometric ratio of baseline
to C800 online wall time is **3.514x**, with a target-bootstrap 95% interval
`[2.661, 4.585]`. C800 was faster on 22 points; 19 point ratios were at
least 2x. The geometric target S3-call ratio was 11.399x. The two repetitions
of each arm had p95 absolute log time ratios of 0.885 and 0.889, respectively
(about 2.42x within-arm spread). The macOS host has no exclusive-core
isolation receipt, and its A/A spread exceeds the 0.10 noise gate. The
controlled wall-time speedup is therefore unknown. The frozen rows remain
useful for the target-work and correctness comparison; a host-isolated paired
run is needed to promote the timing result.
