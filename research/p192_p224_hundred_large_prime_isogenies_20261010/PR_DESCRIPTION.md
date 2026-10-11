Isogeny searches previously required checkout-specific research helpers. This
PR ships a dedicated `isogeny` executable with embedded algorithms, public
curve data, source provenance and independent exact replay. Crypto's ReleaseMe
workflow builds Linux glibc, Linux musl, Apple silicon and Windows archives;
users run the installed binary without Cargo or a checkout.

The continuation found 100 distinct large-prime-degree maps:

| Public isogeny class | Prime degree | New maps |
| --- | ---: | ---: |
| P192 | 2341 | 50 |
| P224 | 1117 | 50 |

These are 100 maps across two prime degrees. Each walk starts on the published
curve; subsequent edges start on independently verified codomains. Every edge
constructs one Frobenius eigenline and reports coverage PARTIAL. Independent
field and polynomial implementations check all kernels, codomains, exact
rational-map identities and subgroup orders, plus 20 public scalar-transport
checks per map (2000 total). Replay occurred on the same host as construction.
The campaign measures operational intervals on a shared host and claims no
matched throughput gain.

`isogeny walk --curve p224 --ell 1117 --count 50 --out ./p224-walk` constructs and
replays a chain. `isogeny verify-walk --input ./p224-walk` independently replays
it from its public anchor. The CLI reports requested and verified counts and
whether the chain is complete. Source/generator mutations remain rejected even
when their manifest hashes are updated. Walk output requires a fresh directory;
the current command does not resume an interrupted directory.

The generalized P192/P224 kernel route checks extension degrees 2 through 16
and is selected automatically above prime degree 1009 when the smaller
Frobenius eigenvalue order is at most 16. Other supported cases retain full
modular-polynomial enumeration. The screen covers primes through 1000000.
The complete screen found 40 P192 and 38 P224 structural candidates with
eigenvalue order at most 64; those counts are not constructed-map counts.

The catalogue has 449 exact models: 304 prime, 72 binary, 60 Koblitz,
4 subfield and 9 extension models. Prime short-Weierstrass construction supports
274 models through 640 bits, and 246 have complete public generator data for
independent replay. Montgomery/Edwards, binary and extension construction
adapters and seven visible standards-source records remain unresolved.
Names, aliases, exact curve/subgroup identifiers, cofactors, covers, route
graphs and browser joins retain their separate capability records. This is the
complete pinned repository inventory. No IC performance ledger row is promoted.

The earlier matched optimization panel remains frozen. Five paired rounds on
identical 349-model inputs compare reference 5cbc3a2dd with candidate ec6240c25
on ARM Linux using Callgrind, counting all traced processes. Screening uses
67.9589% fewer instructions; the complete P192/73 search, including constructor,
receipts and replay, uses 5.2924% fewer. P224/1471 replay uses 0.1478% fewer.
All 15 paired output comparisons pass. Lost counters from two external Docker
interruptions are preserved as unmeasured; no campaign-total gain is claimed.

The current scientific artifacts include all 100 maps and transported public
generators, resource/exit/hash receipts, source and standalone binary archives,
an editable route graph and a visually reviewed three-page PDF. The scientific
producer was clean source 49a280dc6; the archived 449-model executable was clean
source 3c39c5f23. The native audit matches all 891 manifest-bound files against
the preserved canonical-output context. The ten original canonical outputs
are archived from clean 12a846a56, so later global browser updates preserve
the original scientific bindings and manifest bytes.
The archived executable again replays a saved two-step chain outside a checkout.

Validation passed 62 CLI tests, 3 independent field-reference tests,
114 standalone algorithm tests, strict algorithm Clippy, positive chain replay
and source/generator mutation rejection. The complete release library suite
passed **3815 tests, 0 failed, 110 existing tests ignored**, using two concurrent
test threads and localhost access on merged source bff6dcb27. The earlier
pre-merge run passed 3808 tests. The workspace all-target
compile, 30 tests in repaired examples, 14 tests in repaired binaries and the
449-model trait-index freshness check passed. Both staged accounting fixtures
also passed 100 paired repetitions (200 successful fixture executions). The
S4 staged-format fixture passed another 100 isolated repetitions, followed
by the complete merged-tree suite.

The library follow-up repairs one-word truncation of degree-71 linearised
kernels and applying a packed Frobenius orbit walk on models without that
endomorphism. Accounting fixtures now retain enough redundant random rows,
and the calibration fixture checks finite positive measurements while its
exact oracle/work assertions remain. Production solver and benchmark behavior
is unchanged by those fixture corrections. Already merged upstream API and
union-repair files restore the workspace build. Missing sparse-checkout inputs
were restored from their exact tracked bytes. All earlier failed and
interrupted diagnostics remain committed.

The generated traits index now covers all 449 models. Unsupported coordinate
adapters explicitly retain order certification `not_evaluated`; the earlier
121-model guide panel is identified as a historical snapshot. The frozen
scientific manifest continues to match all 891 original bindings; follow-up
integration receipts are additional commit-bound evidence. Main 8b6b1384e is
integrated, with its committed research, ledger and browser additions retained.
The two textual conflicts are resolved; Actionlint and the affected example's
three release tests pass. The original archive audit passes in the extracted
canonical context after integration.

**Publication status:** this ready PR includes the hundred-map source,
catalogue, maps, standalone archive and passing integration receipts. Hosted
checks for the continuation are required before merge. The earlier head
850e77374 passed all four dedicated isogeny platform jobs; those results do
not validate this newer head. An additional strict CLI Clippy diagnostic
retains 35 inherited/generated-module errors and is not presented as passing.
Merge and ReleaseMe publication remain pending; PR dry runs publish build
artifacts rather than GitHub Releases.

See `research/p192_p224_hundred_large_prime_isogenies_20261010/RESULTS.md`,
`FOLLOWUP_VALIDATION.md`, `SUMMARY.json`, `MANIFEST.json` and
`tools/isogeny-cli/README.md` for exact records and binary usage.
