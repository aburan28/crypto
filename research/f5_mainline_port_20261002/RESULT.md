# Exact row-basis rank port on current mainline

Ported the opt-in selected-column rank certificate and eight-table
rank-only kernel from `236205132` onto mainline `51926e7d5`. The newer
`echelon_prefix_counted` and `echelon_resume_counted` functions and tests
remain present. Release tests passed: 9 GF(2) tests, including the new
rank comparison and mainline prefix/resume checks, plus 12 matrix-F5
tests. The normal Cargo dependency fetch failed DNS resolution; the
same release tests and example built with `cargo --offline` from cached
dependencies. The default F5 route remains RREF.

One release binary, SHA-256
`6460a57ad655ae66a5c81ac135b0cd3326180e5b1d7942e6cd5a32a399e16c7e`,
ran the selective-echelon reference and row-basis candidate on the
physical Apple ARM64 host, one Rayon thread. Each seed has one warmup
per mode, five A/A pairs and five alternating A/B pairs. All 44 processes
completed, each with all seven F5 cases. Their raw receipts are
[`FROZEN.json`](FROZEN.json), SHA-256
`1e70621879f750eb6c7a05d0c0a80347a5064ea3fd6b4e7ab74b38b03f52e625`,
and [`HOLDOUT_A.json`](HOLDOUT_A.json), SHA-256
`dda5af0edff44f4569d2f10d773a047137fcb003353cc2ffc78d110093517f00`.

| n24 degree-4 seed XOR | Selective / candidate paired complete-call median | Exact five-pair bootstrap 95% interval | Selective/selective A/A range | Selective marginal median | Candidate marginal median |
| --- | ---: | ---: | ---: | ---: | ---: |
| `0` | **1.450×** | 0.582–1.664× | 0.924–1.728× | 99.863 ms | 71.328 ms |
| `badc0de1` | **1.322×** | 1.253–1.482× | 0.943–1.105× | 72.406 ms | 54.230 ms |

All seven cases matched exact rank, canonical row-space fingerprint,
criterion work, row and column counts, and pruning between modes. The
primary selected 7,436 columns certified all 6,924 rows in every call.
The candidate returns original independent rows: raw row fingerprints
and term counts therefore differ from selective echelon, while row
space remains exact. Counted primary reduction work was about 100.2
million word XORs for selective echelon and 98.7 million for the
candidate. All smaller-case complete-call medians were above their
respective A/A minima. No default route changed.

The local screen misses the further-2× complete-call gate. It is noisy,
especially on the frozen seed, and does not qualify any x86-64 speed
claim. The candidate's marginal reduction medians were 56.453 and
43.456 ms, versus 65.713 and 47.462 ms for selective echelon;
unpacking fell from 24.576 and 20.294 ms to 4.252 and 3.605 ms.
Further improvement must substantially reduce rank work to reach 2×
against selective echelon on this workload. This is a solver-stage
diagnostic, not an IC online or matched-rho result.

The measured source SHA-256 digests are `59462bb330021747b51fc88acf240a918e314b69b8fca34b8c4055e750c29571`
for `gf2_elim.rs`, `6bf4b2c565d572923ac4907e78ced5051889a9251cf7b01bdb9483a5f766a075`
for `matrix_f5_f2.rs`, and `1b7dd569ad3908b84f86aa7b850ebaa8695e54f40e0712ff806d051887c23c9b`
for `koblitz_groebner.rs`. The complete commands, environment,
process outputs and host details are in the raw receipts.

The local `origin/main` ancestry includes PR #801 at merge commit
`35f22378c`; its historical source snapshots remove the earlier
source-hash dependency for new row-building work. The isolated Linux
one- and two-thread replay is prepared in
[`f5-mainline-rank-segments.yml`](../../.github/workflows/f5-mainline-rank-segments.yml)
at `ef4a1515d`, but has not run. It selects the first uncontended
exact attempt on each of four seeds and fails the one-thread gate
unless both the median and bootstrap lower bound exceed 2.00× against
selective echelon.
