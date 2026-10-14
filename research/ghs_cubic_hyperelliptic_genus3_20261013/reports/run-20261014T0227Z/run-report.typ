#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show heading: set text(size: 11pt)

= Run report `20261014T0227Z` (addendum): conjugate-arm transfer with real-model arithmetic, census at $q = 11, 13$

Descriptive record of one run. Study: `research/ghs_cubic_hyperelliptic_genus3_20261013/` (Conductor T-2463). Protocol: `PROTOCOL.md` incl. the 2026-10-14 addendum and its correction (sha256 `fbc83114…`).

== Command and environment
- `hyper_genus3 q <census jsonl>` for $q = 5, 7, 11, 13$ (single-threaded) → `results/run-20261014T0227Z.jsonl`.
- Gate: `census_list 5 0 6`, `census_list 7 0 6` (exhaustive, the merged census code path), compared byte for byte with the q = 5, 7 lines of the merged census run `20261013T0403Z`.
- `census_list q list 6 <A_idx,B_idx> …` for $q = 11$ (15 curves) and $q = 13$ (12 curves) → `results/run-20261014T0227Z.census_q11.jsonl`, `…census_q13.jsonl`.
- Base commit `703e0da8d2a471458b8e1be386e0baec9d811c00`; `hyper_genus3.rs` sha256 `fbe21a4d…`, `census_list.rs` sha256 `61860c2e…` (fork of `cover_census.rs` `26f6997d…`).
- rustc 1.93.1; Darwin 25.6.0 arm64, Apple M4 Pro, 14 cpus, shared host (load average 150–230 during the run), no isolation; census with 6 threads.
- Wall time (practicality note only): hyper_genus3 0, 1, 6, 27 s; gate 178 s; census q = 11: 113 s, q = 13: 393 s. All exited 0.

== Recorded values
- 20/20 instances `all_pass = true`. Split-arm records identical to run `20261014T0029Z`. Conjugate-arm P1–P3/P5 fields are identical too.
- V (8 conjugate instances): enumeration = $\#E_0(K)$ 8/8; associativity 2400/2400; order 800/800; inverse/identity 800/800.
- P4c (8 conjugate instances): homomorphism 3200/3200, doublings 800/800; generator tests 8/8; kernel order 2 (3 instances) or 4 (5); push-forward $= 4P$ on all points of $E_0(K)$ (8,052 points over the 8 instances), 0 failures, 0 not evaluable.
- Gate: BYTE_IDENTICAL, 976 lines, sha256 `65021870…`.
- P5c: structures 1,786,323 (q = 11) and 4,855,539 (q = 13), matching the expected count. 27/27 curves have minimal genus 3, include the type $m = 3, r = 5$ hyperelliptic, match the order, and have 0 mismatches.

#figure(image("addendum.svg", width: 100%), caption: [Recorded values of run `20261014T0227Z` (plotted by `reports/plot_addendum.py`, no computation).])

== Failures
The attempt `20261014T0223Z` (same day, before this run) panicked at q = 11 and 13 on a conjugate-arm point lying over the original $z = infinity$. Its records and annotated env file are kept. The operator stopped it during the census gate. This run is the first after the fix. Development runs were not kept.

Raw records: `results/run-20261014T0227Z.jsonl`, `results/run-20261014T0227Z.census_q11.jsonl`, `results/run-20261014T0227Z.census_q13.jsonl`, `results/run-20261014T0227Z.env.txt`.
