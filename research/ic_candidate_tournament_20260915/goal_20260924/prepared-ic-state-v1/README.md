# Certified reusable preparation for the accepted F5 and SAT controls

Both accepted complete development controls independently reconstruct the same
mathematical preparation, `ICP1hedbff76da644`. Their 61 F5 and 37 SAT ordinary
relations produce full rank 29/29 and identical, independently replayed column
logs on the exact curve `EC1N17Ckb1hbbe2b5b6b1e6`.

| Source preparation | Ordinary queries | Verified relation rows | Rank | Geometric points | Usable IC points | Folded columns |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| Accepted F5 v2 | 216 | 61 | 29/29 | 63 | 62 | 29 |
| Accepted SAT v3 | 149 | 37 | 29/29 | 63 | 62 | 29 |

This is a deterministic preparation certificate, with no new solver or target
execution and no comparative timing. The query counts describe these two
rank-stopped historical preparations; they do not estimate a new fixed-sample
yield or rank the solver families. Original negative-query, failed-attempt,
source and timing evidence remains in the accepted archives linked by each
certificate. The [protocol](PROTOCOL.md) specifies the scope and stop rule.

The mathematical state has full digest
`edbff76da6442b9f2e5e8235682c9bf1052465f310c765ba8a37d018a60bf107`.
The distinct whole-certificate seals are:

- [F5 preparation](f5-preparation.json):
  `d8c5d6679fe89561606154785fdfba763835bf261dceb08fae193eff9f5636ae`.
- [SAT preparation](sat-preparation.json):
  `91856ab78550436d3f668367f9aebd9e2c0604bd64b1472d9d19ec318e2b144e`.

From this repository, verify the small certificates without extraction,
subprocesses, native binaries, Sage or target generation:

```sh
python3.12 research/ic_candidate_tournament_20260915/prepared_ic_state_v1.py verify \
  --certificate research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json \
  --expected-sha256 d8c5d6679fe89561606154785fdfba763835bf261dceb08fae193eff9f5636ae
python3.12 research/ic_candidate_tournament_20260915/prepared_ic_state_v1.py verify \
  --certificate research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/sat-preparation.json \
  --expected-sha256 91856ab78550436d3f668367f9aebd9e2c0604bd64b1472d9d19ec318e2b144e
```

To derive either certificate again from all accepted archive bytes, use
`prepared_ic_state_v1.py freeze --family f5|sat --out NEW_FILE`. It verifies the
externally pinned archive and every member, retains only whitelisted ordinary
inputs, and produces the same corresponding certificate. Immutable writes reject
different content at an existing path. Full certificate seals use the existing
sorted-key compact canonical JSON identity rules; raw formatted file hashes are
separate. The reader itself must be source-bound by a later execution protocol;
these offline receipts do not claim a new preexecution source attestation.
In that protocol, candidate identity binds exact executed source/policy and the
mathematical state digest; the run/preparation manifest binds the full evidence
certificate. Preparation-run seeds and measurements stay outside candidate
identity, including hashes of such run evidence.

The tests reject target-bearing inputs and native `b!=0` collection rows,
torsion removal, duplicate/reordered geometry, false group relations, rank loss,
changed log/column order, inconclusive rows entering the matrix, wrong external
seals, changed archive/provenance, source swaps and measurement/promotion claims.
Changing old target answers or timings has no effect on reusable state.

The next goal gate is a reviewed source-bound target-only adapter and a newly
frozen fresh paired protocol. It must measure one supplied point at a time from
the first target-dependent computation through scalar replay, charge all failed
attempts, retain preparation costs separately, and use the selected incumbent
and same-point rho under identical resources. All such comparative results
remain **not run**; no historical confirmation budget has been reopened.
