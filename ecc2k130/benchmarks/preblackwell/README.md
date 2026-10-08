# Fastest packed clients for g6, g6e, and g4dn

Generate one thin binary per instance family and time it on the matching
Modal GPU. This is walk-rate engineering. The campaign default remains
the sm_120 Blackwell client. `campaign.json` still carries one
`binaryKey`; a mixed-architecture fleet is a further change.

## Boundaries, written before the timings

| family | GPU | SMs | arch | CLMAD default | scaled prior (not a measurement) |
|---|---|---:|---|---|---|
| g7e | RTX PRO 6000 | 188 | sm_120 | 1 | 15.116 B/s measured |
| g6e | L40S | 142 | sm_89 | 1 (after this receipt) | ~5.5 B/s ([ADA-L4-L40S.md](../../ADA-L4-L40S.md)) |
| g6 | L4 | 58 | sm_89 | 1 (after this receipt) | ~2.3 B/s before the 72 W clock cut |
| g4dn | T4 | 40 | sm_75 | 0 | well under 2 B/s ([T4-G4DN.md](../../T4-G4DN.md)) |

Unit: B completed scalar updates/s. Class: engineering. The one-add
floor of 23–25 B/s is a 188-SM number and is not the T4/L4 target.

Success for this thread is a verified median on each named GPU, with
identity matching the requested build, and a thin `make gpu-*` recipe
that emits that build. Promoting `CLMAD=1` on Ada required the Ada
receipt to pick it; it did.

## Measured table

Frozen in this directory. Every sample valid. Field products remain
5.3125 / update. Software Ada rows are the superseded before-mark.

| variant | median B/s | / prior | / 6000 15.116 | / matched software | class | correctness |
|---|---:|---:|---:|---:|---|---|
| RTX PRO 6000 shipping | 15.115792 | — | 1.000 | — | reference | [top-clmad](../top-clmad/summary.json) |
| L40S software | 4.879118 | 0.887 vs 5.5 | 0.323 | 1.000 | engineering | 3/3, CLMAD 0, 128 regs |
| L40S CLMAD=1 | 8.838416 | 1.607 vs 5.5 | 0.585 | 1.811 | engineering | 3/3, CLMAD 1, 94 regs |
| L4 software | 1.323708 | 0.576 vs 2.3 | 0.088 | 1.000 | engineering | 3/3, CLMAD 0, 128 regs |
| L4 CLMAD=1 | 2.490035 | 1.083 vs 2.3 | 0.165 | 1.881 | engineering | 3/3, CLMAD 1, 94 regs |
| T4 software | 0.533211 | 0.267 vs 2.0 upper | 0.035 | — | engineering | 3/3, CLMAD 0, 128 regs |

Cite [summary.json](summary.json). Do not recompute.

## Generate the clients

Needs nvcc 13.3 (the campaign Docker image). No GPU is required to
compile.

```sh
make gpu-g6e          # sm_89, CLMAD=1  →  ecc2k130
make gpu-g6           # same binary as g6e
make gpu-g4dn         # sm_75, software product
```

Campaign bucket (from `ecc2k130/aws/`):

```sh
ARCHES=89 ./build.sh /path/to/ecc2k130            # g6 + g6e, CLMAD defaults to 1
ARCHES=75 CLMAD=0 ./build.sh /path/to/ecc2k130     # g4dn
```

These keys coexist because the published prefix includes arches and knobs.

## Measure

```sh
make bench-g6e-modal          # Modal L40S, CLMAD=1
make bench-g6-modal           # Modal L4, CLMAD=1
make bench-g4dn-modal         # Modal T4, CLMAD=0
ADA_CLMAD=0 make bench-g6e-modal   # Ada software before-arm
# on the instance:
make bench-ada                # L4 or L40S, both CLMAD arms
make bench-g4dn               # T4
```
