# n83 flat x-index: rejected after the dimension-10 gate

The [preregistered gate](PROTOCOL.md) **rejects** the flat 32-bit slot
index. The attempted runtime change is retained as [`REJECTED.patch`](REJECTED.patch),
and `src/cryptanalysis/f6_wide_geometry.rs` was restored to the exact
parent source. This is a four-summand F6-IC component diagnostic; it is
not a complete IC candidate or an IC/rho speedup (`candidate_id: null`,
`IC_online_ms: null`, `rho_online_ms: null`).

The curve is `icv1-f2m83-tm6151469093347-debefd74`. The full standard
dimension-12 cofactor-projected factor base has 4,054 distinct usable
points, 2,027 sign-folded columns, 8,219,485 unordered pairs and
4,108,723 signed-sum representatives. Every ordinary probe queried the
same public T001 point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`. Every process
exited zero, built exactly the same number of representatives, returned
the same exact `no_witness`, and reported correctness `PASS`.
The separate full-base planted `[0,2,4,6]` control found that same
witness through both portable and PMULL paths and replayed it in the
group. All nine focused release geometry tests passed.

## Frozen full-base A/B panel

All times below are **exploratory on a contended physical Apple M4 Pro**.
The panel used the frozen #1421 full binary and the frozen candidate
binary, with no build between arms. The table retains every execution.
Index build is reusable target-independent setup; the exact query is
target-dependent work. Peak RSS includes both phases.

| Order | Arm | Build (ms) | Query (ms) | Peak RSS (MB) |
| ---: | --- | ---: | ---: | ---: |
| 1 | #1421 baseline | 5,494.900 | 1,917.701 | 545.767 |
| 2 | flat slots | 3,756.022 | 2,096.038 | 300.958 |
| 3 | flat slots | 4,035.176 | 1,849.063 | 305.676 |
| 4 | #1421 baseline | 4,268.223 | 1,747.760 | 549.143 |
| 5 | #1421 baseline | 4,295.978 | 1,322.653 | 549.159 |
| 6 | flat slots | 3,303.453 | 1,284.663 | 305.824 |
| 7 | flat slots | 2,435.870 | 1,303.071 | 305.660 |
| 8 | #1421 baseline | 3,495.098 | 1,149.570 | 548.995 |
| 9 | #1421 baseline | 3,368.010 | 1,726.655 | 548.995 |
| 10 | flat slots | 2,117.043 | 1,051.328 | 305.644 |

The five-process query medians are 1,726.655 ms baseline and
1,303.071 ms candidate (24.5% lower); the build medians are
4,268.223 and 3,303.453 ms (22.6% lower). Maximum RSS is
549,158,912 versus 305,823,744 bytes (44.3% lower). Yet adjacent
same-target query pairs are mixed: the candidate is slower in three
of five pairs. The baseline's own query range, 1,149.570–1,917.701
ms, is much wider than the median difference. The median alone is
therefore not controlled evidence of a query speedup.

## Small-base rejection

The preregistered small screen ran baseline/candidate/candidate/baseline.
Each process internally repeated its query three times; the entries
below are the middle of those three times. All cases were exact misses
with matching representative counts and `PASS` correctness.

| Dimension | Baseline process medians (ms) | Candidate process medians (ms) | Two-process medians (ms) | Decision |
| ---: | --- | --- | --- | --- |
| 8 | 1.616, 4.012 | 4.859, 1.660 | 2.814 vs 3.259 | 15.8% slower; within 20% gate |
| 10 | 55.416, 56.910 | 195.628, 98.188 | 56.163 vs 146.908 | 161.6% slower; **fails 20% gate** |

The dimension-10 candidate also had larger build times in both process
comparisons, and the query regression persisted in its second process.
Under the frozen rule, this rejects the change even though the noisy
full-base medians and memory favored it. No additional replay was used
to redefine the gate after observing the failure.

The candidate geometry source SHA-256 was
`8d18c52a4554d0b3814621777235e7e71ea5696a48b00626738818cf9b53da18`.
The candidate full, small and planted release binaries have SHA-256
`f9e4380b6ba025b2465cd0e9e9c22b5fea2946fae7d86d3092084b874c82a629`,
`8152bd95262f56c4caa2a956d46163f11b33f2cb0b3df3242fb6d053930bea56`
and `5fe4ce6480aa062518b513038dd8496d7fc939572d0b357d513eda69c8987bcc`.
The baseline hashes and parent head are pinned in the protocol.
The [runner](run.sh), `status.tsv`, all 14 ordinary JSONL outputs and
stderr files, planted output, build log and focused test log are in this
directory. [`SHA256SUMS`](SHA256SUMS) pins the committed evidence.

No end-to-end IC online interval, complete n83 F6 decomposition or
same-point rho solve was measured here. The K0 four-summand natural-target
coverage ceiling remains `4.662e-12`; a faster exact miss would not
establish a useful relation rate. The next lookup experiment should
screen the smaller base before spending full-base runs, while the
higher-arity ordinary-relation solver remains the main algorithmic gap.
