# Pair-swap quotient preserves n53 rank rows while halving index S3 calls

The exact transformation
`(left,right,relative) -> (right,left,-relative mod 53)` reduces the
23,320-point K220 compact-orbit index from **2,565,200 to 1,282,710 S3
construction calls**, a 49.9957% reduction. It restores the swapped root
aliases in the original scan positions. The complete n53 root-table test
finds the same 2,564,528 key/value entries, and the frozen six-seed run
produced the same 220 novel rank rows in each paired arm. All 12 IC cells and
six same-point strong-rho cells independently replayed the scalar
`12641099082383` for public Q `[6017759979451042,2886969156937]`.

The field is `GF(2)[x]/(x^53+x^6+x^2+x+1)`, curve `EC1N53Ce0hb097de99be9a`,
and prime subgroup order `21044858204113`. The base has 23,320 usable points
and 220 signed-Frobenius columns; its BLAKE3 digest is
`7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973`.
The [control](candidate_control.json) and [pair-symmetric](candidate_symmetry.json)
manifests fix the exact source and algorithm identities. `KIC_PAIR_SYMMETRY`
was the only algorithmic variable; the rank probe cap was absent. The Q,
candidate IDs, six rank/rho seed pairs, binaries, and 60-second/16-GiB
per-cell envelope were committed before the timed run in
[HELDOUT_FROZEN.json](HELDOUT_FROZEN.json). The source SHA-256 was
`b7b290961b0137ce728b80b9ff92934013328c930dbf32fe0be00c01ae788471`.

| Rank seed | Saved rank S3 calls | Same 220 rows | Paired index-time ratio | Paired complete-cold ratio | Paired target-online ratio | Observed rho / pair-online |
| ---: | ---: | :---: | ---: | ---: | ---: | ---: |
| 532053 | 12,296 | yes | 0.330 | 0.783 | 1.087 | 3.289 |
| 532054 | 11,554 | yes | 0.566 | 1.007 | 1.103 | 2.809 |
| 532055 | 8,268 | yes | 0.525 | 0.931 | 1.021 | 3.939 |
| 532056 | 19,292 | yes | 1.848 | 1.021 | 1.391 | 4.620 |
| 532057 | 19,292 | yes | 1.187 | 1.910 | 3.337 | 0.725 |
| 532058 | 5,512 | yes | 0.503 | 0.862 | 0.982 | 2.775 |
| Paired median | 11,925 | 6/6 | 0.545 | 0.969 | 1.095 | 3.049 |

The target-online interval starts with the first target-dependent query and
ends after scalar recovery and its in-process group check. It sums target
query, PDP, relation check, descent, and recovery-check phases. The rho
interval starts its target-dependent walk and ends after recovery/check on
the same Q. Base construction, index construction, rank relation collection,
matrix construction, and final elimination are in the separate complete-cold
IC interval; every failed attempt and its cost remains charged there.

The pair quotient saved 5,512–19,292 rank-stage S3 calls per seed, a median
0.1150% of rank calls. Target search generally ends before a swapped state
would be revisited, so the exact index construction reduction did not produce
a corresponding target-query reduction. Across this shared-host run, median
observed peak IC RSS was 332.2 MiB (ordered) versus 303.1 MiB (paired). Median
target-online times were 45.202 ms (ordered) and 49.864 ms (paired), while
same-point rho's median was 151.670 ms. The median paired target-online
ratio was 1.095, with a descriptive six-draw bootstrap interval of
`[1.002,2.364]`. The median index-time ratio was 0.545
`[0.416,1.518]`, and complete-cold ratio was 0.969 `[0.822,1.465]`.
Seed 532057 had a 139.769-ms paired target query versus 41.883 ms in its
ordered arm, and its same-point rho took 101.384 ms. These host timings are
exploratory under the repository CPU-isolation gate; the controlled online
ratio remains unknown. Exact operation counts and row equality pass the
preregistered gate independently of that timing variance.

The next fixed-base comparison should target the hot second S3 root call in
rank and target extraction, where each rank run made 8.3–11.0 million calls.
First count field multiplications and inversions in that call under the same
rank trace, then test a fixed-target batched or shared-subexpression kernel
against the byte-identical 220-row and single-Q scalar controls. Promote it
only after a charged rank/target/rho replay with an isolated-host receipt.
This directs work toward per-query cost, while the pair quotient remains a
verified index-memory and construction optimization for larger bases.

All 198 raw cell files, including stderr and timeout/status receipts, are
losslessly retained in [`heldout_runs.tar.gz`](heldout_runs.tar.gz). The
[manifest](heldout_runs_manifest.json) names every file, its SHA-256 and
byte count; duplicate base dumps are stored once and reconstructed during
verification. `python3 archive_runs.py` checks the archive, independently
replays every rank, target, and rho scalar from reconstructed raw cells, and
runs `analyze_heldout.py --check`. The computed
[analysis](HELDOUT_ANALYSIS.json) contains every phase, resource measurement,
paired ratio, and the exact-operation gate. Archive SHA-256:
`7f38f6a451adb59b4ab3041b1662a9143fb2f85944d0ecb63ee25d857bc0f797`.
