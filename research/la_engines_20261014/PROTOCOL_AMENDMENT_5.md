# Amendment 5 (2026-10-15, before any data for these arms)

## Why

Amendment 4 measured the ported CADO-style filter in one session against
CADO-NFS: at D = 262,144 CADO-NFS was still about 1.26× faster than ours
(16.0 s against 20.2 s, seed 1), and our filter dominated our pipeline.
Phase timers (development instrumentation, removed before this measurement)
put the filter's time in duplicate removal (a hash set of cloned rows), the
column-list build, and the purge's full sort of every component at every
step. The change:

- duplicate removal hashes rows in parallel and compares only rows whose
  hashes collide (the same rows survive as with the set);
- column lists are reserved to their final length before the build;
- the purge computes per-row Ω₃₁ weights in parallel (summed per component
  in row order, as before), selects the heaviest `2·need + 1,024`
  components with a partial selection and sorts only that prefix (the rest
  is sorted only if the prefix runs out).

The filter's output is unchanged: on D = 4,096 and 65,536, seeds 1–2, for
`legacy-bw`, `cado-bw` and `cado-bw-pos`, the new binary's filter report,
core dimension, core nonzeros and verification equal the amendment-4
records exactly (12/12 cells; a development check, not data).

## Arms (m = 3, seeds 1–3, four threads, 1,800 s timeout, interleaved by
arm within each (D, seed), one fresh session)

| arm | as in amendment 4 |
|:--|:--|
| `legacy-bw`, `cado-bw`, `legacy-bw-pos`, `cado-bw-pos`, `cado-nfs` | same definitions, new binary |

D = 1,024 … 262,144 in powers of two. Headline ratios per (D, seed) in this
session: `legacy-bw / cado-bw`, `legacy-bw-pos / cado-bw-pos` and
`cado-nfs / cado-bw-pos`.
