Isogeny searches previously required checkout-specific research helpers. This
PR ships a dedicated `isogeny` executable with embedded algorithms, public
curve data, source provenance and independent exact replay. Crypto's ReleaseMe
workflow builds Linux glibc, Linux musl, Apple silicon and Windows archives;
users run the installed binary without Cargo or a checkout.

The completed local continuation found 100 distinct large-prime-degree maps:

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

The local catalogue has 449 exact models: 304 prime, 72 binary, 60 Koblitz,
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
source 3c39c5f23. A fresh native audit matches all 891 manifest-bound files.
The archived executable again replays a saved two-step chain outside a checkout.

Scoped validation passed 62 CLI tests, 3 independent field-reference tests,
114 standalone algorithm tests, strict algorithm Clippy, positive chain replay
and source/generator mutation rejection. The library follow-up also passes
57 wide-arithmetic tests, the degree-71 constructor, both isogenous-model
tests and isolated cover/Edwards/gate checks. It fixes two reproduced library
integration defects: one-word truncation of degree-71 linearised kernels and
applying a packed Frobenius orbit walk on models without that endomorphism.
The j=0 accounting fixture now has 32 extra random rows; solver behavior is
unchanged. Earlier failing and interrupted diagnostics remain committed.

**Publication status:** the hundred-map continuation commits 3c39c5f23 and
2a3f7f2ac are local and have not yet been pushed. The current remote head is
850e77374. The required `cargo test --release --lib` rerun uses two concurrent
test threads and localhost access and is still running. The owner's rule
requires fixing failures before pushing; that gate is in place. An additional
strict CLI Clippy diagnostic retains 35 inherited/generated-module errors and
is not presented as passing. The earlier remote head passed all four dedicated
isogeny platform jobs; those results do not validate the unpushed continuation.
Hosted checks, merge and ReleaseMe publication remain pending.

See `research/p192_p224_hundred_large_prime_isogenies_20261010/RESULTS.md`,
`FOLLOWUP_VALIDATION.md`, `SUMMARY.json`, `MANIFEST.json` and
`tools/isogeny-cli/README.md` for exact records and binary usage.
