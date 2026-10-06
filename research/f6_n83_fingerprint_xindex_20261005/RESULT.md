# n83 fingerprinted x-index: exact, smaller, but rejected by query gate

The [preregistered protocol](PROTOCOL.md) rejects this second lookup
variant. Its attempted runtime change is preserved as
[`REJECTED.patch`](REJECTED.patch); the source file was restored to the
retained #1421 hash. This was an unisolated four-summand F6-IC component
experiment, with no complete IC candidate or end-to-end speedup
(`candidate_id: null`, `IC_online_ms: null`, `rho_online_ms: null`).

The exact curve is `icv1-f2m83-tm6151469093347-debefd74`. The full
standard dimension-12 cofactor-projected base has 4,054 actual usable
points, 2,027 signed columns, 8,219,485 unordered pairs and 4,108,723
signed-sum representatives. The public T001 point was
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`. All 14 ordinary
probe processes exited zero, returned the same exact `no_witness`, and
matched the baseline representative counts. Nine focused release tests
passed. The full-base planted `[0,2,4,6]` control exited zero, found
the same witness via portable and PMULL queries, and replayed it in the
group. Planted correctness does not estimate natural relation yield.

## Small-base gate passed

The A/B/B/A screen used the frozen #1421 binary and fingerprint candidate.
Each process contained three internal repeats; each value below is its
middle query time. All outcomes were exact misses with `PASS` correctness.

| Dimension | Baseline process medians (ms) | Candidate process medians (ms) | Two-process median, baseline / candidate (ms) |
| ---: | --- | --- | --- |
| 8 | 9.325, 5.644 | 2.894, 4.887 | 7.484 / 3.891 |
| 10 | 196.568, 226.111 | 125.701, 116.463 | 211.340 / 121.082 |

Neither dimension violated the frozen 20% regression limit, so the
full-base panel was run. These are exploratory timings on a busy host.

## Full-base gate failed

The ordered A/B/B/A/A/B/B/A/A/B panel used the same public target and
frozen binaries, with no build between arms. Index build is reusable
target-independent preparation; query is target-dependent. Peak RSS
includes both phases. All values are **exploratory**.

| Order | Arm | Build (ms) | Exact query (ms) | Peak RSS (MB) |
| ---: | --- | ---: | ---: | ---: |
| 1 | #1421 baseline | 12,598.266 | 6,142.704 | 473.104 |
| 2 | fingerprint | 8,761.214 | 3,983.899 | 333.431 |
| 3 | fingerprint | 5,122.465 | 3,002.131 | 333.103 |
| 4 | #1421 baseline | 7,202.013 | 3,779.223 | 542.261 |
| 5 | #1421 baseline | 13,646.855 | 4,886.307 | 535.347 |
| 6 | fingerprint | 8,725.750 | 3,989.577 | 289.620 |
| 7 | fingerprint | 6,708.425 | 5,601.472 | 243.909 |
| 8 | #1421 baseline | 9,203.513 | 2,190.464 | 544.309 |
| 9 | #1421 baseline | 3,537.290 | 1,863.665 | 544.834 |
| 10 | fingerprint | 5,660.341 | 2,720.736 | 334.840 |

The five-process exact-query medians are **3,779.223 ms baseline** and
**3,983.899 ms candidate**, so the candidate is 5.4% slower and misses
the frozen requirement of at least 10% lower query time. The build
medians are 9,203.513 and 6,708.425 ms (27.1% lower candidate), and
maximum RSS is 544,833,536 versus 334,839,808 bytes (38.5% lower).
Adjacent query pairs favor the candidate in three of five comparisons,
but the last two reverse direction. The baseline's own query interval
ranges from 1.864 to 6.143 seconds, a 3.3-fold spread. This host does
not meet the isolation or noise gate for a controlled CPU ratio.

The candidate source SHA-256 was
`d870b3c67f5f6b745fbf0ca679dcd8e1b3706a7f55e557b41ac81d76ac4491d8`.
Its full, small and planted binary hashes were
`a3105b8e6f284fbddb714d6cdda0a6c1434e5492670580a42cae3ceda57f77d8`,
`eca153b4fec91f730c3c7d3632360aa6da571fbfea556ca68c2efcec2f273725`,
and `5a1e522ff3eb745ff476c1f4d5f83e80ddaef1426bff287aabdcfec7e7a84803`.
The protocol pins the baseline source and binaries. The two runners,
all ordinary JSONL/stderr files, exit status tables, planted control,
build and test logs, rejected patch, and [`SHA256SUMS`](SHA256SUMS) are
retained here.

Decision: **reject** under the frozen full-base query gate. The rejected
change is not part of the runtime branch. No new n83 ordinary relation
solver, full F4/F5/F6 comparison, one-target IC online interval,
same-point rho run or IC speedup was established. The four-summand
uniform-target coverage ceiling remains `4.662e-12`, so lookup
constants cannot by themselves solve the higher-arity relation gap.
