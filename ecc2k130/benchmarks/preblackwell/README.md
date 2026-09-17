# Fastest packed clients for g6, g6e, and g4dn

Generate one thin binary per instance family and time it on the matching
Modal GPU. This is walk-rate engineering. The campaign default remains
the sm_120 Blackwell client.

## Boundaries, written before the timings

| family | GPU | SMs | arch | CLMAD default | scaled prior (not a measurement) |
|---|---|---:|---|---|---|
| g7e | RTX PRO 6000 | 188 | sm_120 | 1 | 15.116 B/s measured |
| g6e | L40S | 142 | sm_89 | 0 | ~5.5 B/s ([ADA-L4-L40S.md](../../ADA-L4-L40S.md)) |
| g6 | L4 | 58 | sm_89 | 0 | ~2.3 B/s before the 72 W clock cut |
| g4dn | T4 | 40 | sm_75 | 0 | well under 2 B/s ([T4-G4DN.md](../../T4-G4DN.md)) |

Unit: B completed scalar updates/s. Class: engineering. The one-add
floor of 23–25 B/s is a 188-SM number and is not the T4/L4 target.

Success for this thread is a verified median on each named GPU, with
identity matching the requested build, and a thin `make gpu-*` recipe
that emits that build. Promoting `CLMAD=1` on Ada requires the Ada
receipt to pick it.

## Generate the clients

Needs nvcc 13.3 (the campaign Docker image). No GPU is required to
compile.

```sh
make gpu-g6e          # sm_89, software product  →  ecc2k130
make gpu-g6           # same binary as g6e
make gpu-g4dn         # sm_75, software product
```

Campaign bucket (from `ecc2k130/aws/`):

```sh
ARCHES=89 CLMAD=0 ./build.sh /path/to/ecc2k130    # g6 + g6e
ARCHES=75 CLMAD=0 ./build.sh /path/to/ecc2k130    # g4dn
```

A mixed-architecture fleet is still a further campaign-record change.
These keys coexist because the published prefix includes arches and knobs.

## Measure

```sh
make bench-g6e-modal          # Modal L40S
make bench-g6-modal           # Modal L4
make bench-g4dn-modal         # Modal T4
# on the instance:
make bench-ada                # L4 or L40S, both CLMAD arms
make bench-g4dn               # T4
```
