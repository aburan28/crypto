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

## K4096 feasibility diagnostic

The predeclared 900-second local diagnostic completed in 806.573 seconds of
outer cold wall time. Its independent checker reconstructed all 679,936
usable factor-base points, 1,392,477,548 exact root-index entries, the
four-point target witness `[212723,543487,465789,600458]`, and scalar
3,640,978,018,804,724. The peak sampled process RSS was 18,615,615,488
bytes and the lowest sampled system free-memory reading was 59%; neither
stop gate fired. The completed target-online interval was 49.886 ms on
this local host, recorded as a diagnostic rather than a paired speed claim.

The [K4096 analysis](analysis_k4096_20261016.json) has SHA-256
`8da0af76705063e04ef4dd82309fd007ba13df0d9dd1068490028be8fe5b08ca`.
The raw archive at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-id30-k4096-20261016-1238.tar.gz`
has SHA-256
`41b7be597432218a38dc74081f09c38bfa0d26f4d10570f213ed6cbcde5bd289`.
The [IC1 run row](ic1_row_k4096_20261016.json) binds the candidate
`IC1N83Ckb1a1fb679936PDP4rootRCguidedLAgaussTDdirectISO0ha0487e78541e`
to the same exposed one-point workload. The exact source SHA-256 is the
earlier `f14d9747...` control source; the run occurred before the subsequent
lint-only source edit, so this receipt is retained under its executed source.

The final lint-clean source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and rebuilt binary SHA-256
`761b97c2c3b04c464d2c34b469168eae4cba5f26f9a86ab0b5163bd794547eb2`
were used in a separate K4096 confirmation attempt. It reached the frozen
900-second cap at 900.971 seconds outer wall, with exit `-15`, a 16,766,222,336
byte sampled peak RSS and a 42% minimum free-memory reading. The resource
stops did not fire. The [bounded analysis](analysis_k4096_postlint_20261016.json)
has SHA-256 `0cafdd108bdb238438bdf21f2d01625c6ee1b6955c4ddc0daf86b3abe8474ddf`;
its raw archive at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-id30-k4096-postlint-20261016-1256.tar.gz`
has SHA-256
`68fa76d41205226ba3d03b29df80ac9e77da809f57f9f3a7ac83ad47437b2b9c`.
This attempt has `BOUNDED_UNKNOWN` status and no completed base or target
row. Its [bounded diagnostic row](bounded_row_k4096_postlint_20261016.json)
keeps the candidate and online fields null. The earlier exact K4096 result
remains bound to its original source and binary.

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
The same full suite on the final lint-clean source passed 4,148 tests with
117 ignored and zero failures in 1,067.77 seconds; log SHA-256
`1e6908a796ae40de903f9e904a772860bad7068a3de2800dfa075b1f2a7af5f9`.
Its touched release example passed all 21 tests.
`cargo clippy --locked --example koblitz_orbit_dlp_fast -- -D warnings`
exited zero on the final source. The repository's existing command-line
allowance for `clippy::chunks_exact_to_as_chunks` emitted an unknown-lint
warning under the local toolchain.
`rustfmt --check`, `actionlint` for the dedicated workflow, Python syntax,
and `git diff --check` passed.

The next discriminating run for the larger base is a final-source K4096
completion on a sufficiently stable host, followed by a same-source,
same-point Tier 1 A/A-controlled online comparison against the strongest
supported smaller-base finalist, with each exact base and scalar replayed.
The separate current-source online range screen chooses the factor base for
the already prepared five fresh-point Tier 1 comparison; it used the prior
29-bit producer source and is kept as its own workload.
