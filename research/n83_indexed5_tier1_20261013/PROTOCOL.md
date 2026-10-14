# N83 indexed5 cloud-instance capacity and same-point gate

The merged five-byte index makes a larger prefix factor base fit a smaller
root table. This gate measures whether K1798 completes on a fresh standard
16 GB Linux VM, then measures K1500 and K1271 on the same public point and
resource envelope. It also runs the strong-rho reference on that point.

## Frozen inputs and code

- Source parent: `759978250c1ce53d7b87bd8a1c2304781fcc4fc4` (merged PR
  #1787). Producer source SHA-256:
  `64171575e58b97c4ad9041b4b49e5f2b5697d26f71be73861fc36ea879db9161`.
  Cargo.lock SHA-256:
  `ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90`.
- Public point: `research/n83_indexed5_20261013/primary_target.jsonl`,
  SHA-256 `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The producer receives only this point. The separate checker sidecar is
  `answer.jsonl`, SHA-256
  `651eb7d0ef08a1655036746311b19b6989e340ca4752b6d675606055e59c12a1`.
- Curve: Koblitz `a=b=1` over `F_2^83`, polynomial low terms `[0,2,4,7]`,
  subgroup order `8569786107849059`, cofactor `1128547018`, generator
  `(6614244651398184860292015,8868260410029249391539529)`.
- Source-pinned flags: `KIC_INDEXED_ROOT5=1`, `KIC_STATE28=1`,
  `KIC_PACKED_ROOT_TABLE=1`, Bloom/dense/load75/sparse-canonical off,
  rank seed 29, four rank threads, one target thread. Each arm starts with
  empty process caches and `construct:83:1:K`.

## Ordered arms and stop rules

| Order | Arm | K | Actual base if completed | Root-slot allocation bound | Wall cap |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | K1798 indexed5 | 1,798 | 298,468 | 5 GiB | 1,800 s |
| 2 | K1500 indexed5 | 1,500 | 249,000 | 5 GiB | 1,200 s |
| 3 | K1271 indexed5 | 1,271 | 210,986 | 2.5 GiB | 1,200 s |

The table size and usable-point entries are calculations or archived prefix
counts, not new measured outcomes. A process is stopped if sampled RSS
exceeds 14 GiB or its cap expires. The runner records exit status, child CPU,
peak RSS, sampled RSS, cgroup memory/CPU counters, stdout/stderr, source and
binary hashes, base dump, and target row. It attempts every arm even if an
earlier arm fails; an infrastructure stop is a censored row.

The strong-rho reference uses the archived `fold-L64-DP6-E0` policy, seed
`2026101301`, one OS thread, the same public point, and a 120-second cap.
Its config SHA-256 is
`cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b`.

The standard public `ubuntu-24.04` runner is an instance-scoped Tier 1
environment. This first gate uses one run per arm, so its CPU timing ratios
are exploratory. A qualified Tier 1 speedup requires the separately
specified five alternating same-point pairs and noise/interval gates in the
cryptanalysis tiered benchmark service. The primary metric here is each
independently replayed single-target online interval; setup, base/index,
rank, and final matrix costs remain separate. The independent verifier checks
the curve/subgroup, complete point-defined base, four-point relation, scalar,
phase sum, and archived first-witness parity.

The branch-scoped workflow uploads the raw receipt even when a build, arm,
or verifier fails. It does not alter repository-wide Actions settings or
send credentials to the benchmark children.
