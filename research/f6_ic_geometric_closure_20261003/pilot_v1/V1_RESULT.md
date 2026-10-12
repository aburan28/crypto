# F6-IC v1 pilot: fewer reductions, limited full-call gain

The preregistered v1 gate enumerated the exact usable base when one of three
summand x-codes was fixed. Source commit
`30bc3faa8ccdca79b05d64902ed9a70aa4bbc430`, binary SHA-256
`c0ec41f10d5acde9c33dbce2dd42e42d4752866ff55e25f773cdbafc37b72fdf`,
workload ID `f57e562ed2e0`, and both candidate IDs are in [`freeze.txt`](freeze.txt).
The new public point was `[130473,87558]`; no target scalar was constructed
by the fixture. This was an IC-variant pilot, with no new rho run.

| Repetition | Inherited F4 online ms | F6-IC v1 online ms | F4 / F6 total reductions | F4 / F6 total splits | F6 residual lookups | F6 point additions |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 48.356 | 46.330 | 379 / 163 | 159 / 87 | 8,928 | 17,859 |
| 2 | 47.810 | 46.099 | 379 / 163 | 159 / 87 | 8,928 | 17,859 |
| 3 | 47.813 | 46.246 | 379 / 163 | 159 / 87 | 8,928 | 17,859 |

All six fresh processes exited zero, completed three target queries, recovered
the same scalar `64024` and independently replayed it. The first two queries
were proved to have no usable relation; the third yielded a verified witness.
All five exclusive online phase costs summed exactly to each charged online
interval. F6 made 72 exact geometric branch refutations and one geometric
witness per call. Its reductions fell by 57% across the three attempts, but
its 17,859 general-representation point additions consumed most of that gain.
This is an engineering diagnostic, not an algorithmic crossover.

The Mac is not host-isolated. These CPU wall times are exploratory; the
controlled speedup and IC-versus-rho ratio remain **unknown**. Raw reports,
stderr, exit status and timestamps are in [`runs/`](runs/). The keyed
[`v1-measurements.jsonl`](v1-measurements.jsonl) retains each full attempt
ledger, five online phases, correctness and raw-output SHA-256.

Decision: preserve v1 and test packed exact curve arithmetic inside the same
branch enumeration. The one-fixed geometry is useful, but needs a cheaper
point kernel before it can plausibly reach the requested full-call target.
