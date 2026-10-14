#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show heading: set text(size: 11pt)

= Run report `20261014T0029Z`: hyperelliptic genus-3 GHS cover build and transfer checks

Descriptive record of one run. Study: `research/ghs_cubic_hyperelliptic_genus3_20261013/` (Conductor T-2463). Protocol: `PROTOCOL.md` incl. the 2026-10-13 amendment (sha256 `6bc81636…`).

== Command and environment
- `rustc --edition 2021 -O hyper_genus3.rs -o hyper_genus3`; `./hyper_genus3 q ../ghs_cubic_cover_census_20261012/results/run-20261013T0403Z.jsonl` for $q = 5, 7, 11, 13$, appended to `results/run-20261014T0029Z.jsonl`.
- Base commit `12331a9a7b15781b3970c6220636951f5b56119f` (pre-rebase; the same tree is `f69cec000` on the PR branch); source sha256 `be4c3374cb64504e23da4ea0f48120dd6608d6a6c31d2c3bf798def5e139b68a`; census input sha256 `5f087a57…`.
- rustc 1.93.1; Darwin 25.6.0 arm64, Apple M4 Pro, 14 cpus, shared host (load average about 70 at start), single-threaded tool, no isolation. Seed 20261013, deterministic.
- Wall time (practicality note only): q = 5: 0 s, 7: 1 s, 11: 8 s, 13: 24 s. All four exited 0.

== Arms and workloads
Split arm ($Q = (x-a)(x-b)$): 3, 4, 3, 2 instances at $q = 5, 7, 11, 13$. Conjugate arm ($Q$ irreducible): 2 per field; at $q = 7$ both chosen by the amendment rule (no candidate with cofactor 4 or 8). Exhaustive candidate scan for the $4 | \#E_0(K)$ check: 400, 2352, 24200, 56784 candidates per arm.

== Recorded values
- 20 of 20 instances `all_pass = true`; totals records report `P5b_all_candidates_0mod4 = true` at every $q$.
- P1: 200 random points per instance, all five point checks 200/200. P2: degree 7 (split) or 8 (conjugate), squarefree.
- P3: $N_k$ on the $FF_q$ model = $K$-model count = prediction on 20/20; Mumford enumeration = $\#E_0(K)$ on 12/12.
- P4 (split, 12 instances): homomorphism 4800/4800 pairs, doublings 1200/1200; generator tests 12/12; kernel order 2 (8 instances) or 4 (4 instances); push-forward $= 4P$ on all of $E_0(K)$, 0 failures (preregistered subset: 0 failures, the rest not evaluable there).
- P5a: census records found and consistent for all 11 instances at $q = 5, 7$; none found at $q = 11, 13$ (outside the census sample).

#figure(image("checks.svg", width: 100%), caption: [Recorded values of run `20261014T0029Z` (plotted by `reports/plot.py`, no computation). Left: Jacobian orders against $\#E_0(K)$, 20 instances. Right: per split instance, $|E_0(K)|$ and the number of points with $pi_(0*) phi(P) = 4P$ verified on the preregistered set and on all of $E_0(K)$.])

== Failures, superseded runs
No failures and no interrupted runs. Run `20261014T0027Z` (same code except the record key) is kept and superseded: its instance records contain the key `"P4"` twice. Records of the two runs agree after renaming that key (diff of everything but wall times). Development runs before the amendment were not kept.

Raw records: `results/run-20261014T0029Z.jsonl`, `results/run-20261014T0029Z.env.txt`.
