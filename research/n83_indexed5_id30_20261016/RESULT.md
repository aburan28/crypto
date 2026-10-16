# N83 exact five-byte root index with 30-bit state IDs

The opt-in `KIC_INDEXED5_ID30=1` layout retains five bytes per slot and exact
canonical-root equality while moving one fingerprint bit to the state ID.
For the N83 symmetric one-root state formula
`83*K*(K-1)/2 + 42*K`, the largest representable folded-column count rises
from K3596 to K5086. The [Lean proof](StateId30Capacity.lean) checks both
boundaries and the all-K upper exclusion; the [proof map](PROOF_MAP.md)
records its exact scope and `#print axioms` results.

The [frozen K1024 control](PROTOCOL.md) ran two fresh processes on one
previously exposed public point. Independent replay of both arms passed:

| Encoding | Actual usable base points | Exact root entries | Base SHA-256 | Four-point witness | Recovered scalar |
| --- | ---: | ---: | --- | --- | ---: |
| Legacy 29-bit ID | 169,984 | 87,028,652 | `f3bb5adfd099df1be8e90edf56c9a764028e060d7e0c3a97759e9b699c423629` | `[25135,120536,27047,45109]` | 3,640,978,018,804,724 |
| Opt-in 30-bit ID | 169,984 | 87,028,652 | `f3bb5adfd099df1be8e90edf56c9a764028e060d7e0c3a97759e9b699c423629` | `[25135,120536,27047,45109]` | 3,640,978,018,804,724 |

The independent checker rebuilt the factor base, summed the four witness
points to the public target, replayed the scalar against the generator,
checked the exact index layout and source/binary/input hashes, and reconciled
the five exclusive target-online phases. The [versioned analysis](analysis_control_20261016.json)
has SHA-256 `69f77eaaae6e08185046c1701450e03e032c5892e6644dc370e2d15c07e58552`.
The full raw control archive is
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-id30-control-20261016-1225.tar.gz`,
SHA-256 `f831bdb50bff9dd27e47c5e21554aed6358d42b5d62859e10cc4ba147c85ddf7`.
The producer source SHA-256 is
`f14d9747cd24627d39d7cae4c9bef69ed2819544125ef15f31306ce5d75f3fc7`.
The [IC1 run rows](ic1_rows_control_20261016.json) bind the two exact
candidate manifests to the same one-point workload `eebf98ee1031`:
legacy `IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0h8eb8e220bff6`
and opt-in `IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0hd0bd0a9c5003`.
The materializer checks the archive and replay hashes and was run twice
without changing the versioned records.

The legacy and opt-in local target-online intervals were 56.980 and 65.056 ms;
their outer cold process intervals were 226.786 and 127.693 s. These are
ordered, shared-host correctness controls, with simultaneous test activity;
`performance_claim` and `selected_factor_base` remain `null` in the analysis.
The focused release example suite passed 21 tests, including ID boundary and
fingerprint-collision cases; log SHA-256
`ebb49ea05f9e0f4f0a218566ab588c4372f6c1205c07239d6610d823576d75f1`.
The full serial release library suite passed 4,148 tests with 117 ignored
and zero failures; log SHA-256
`1d325404c9c978287ec3d1e245a8a9ef12529d902d4fd30f2b7c9f3918b82840`.
`rustfmt --check`, `actionlint` for the dedicated workflow, Python syntax,
and `git diff --check` passed.

The next discriminating run is the pinned K4096 feasibility diagnostic on a
host with sufficient RAM, followed by paired Tier 1 one-target timing if the
exact base, relation and scalar replay pass. The separate current-source
online range screen chooses the factor base for the already prepared five
fresh-point Tier 1 comparison.
