# Turing (EC2 g4dn) — software packed walk on sm_75

EC2 **g4dn** carries the **T4**, `sm_75`. This tree had no client recipe
and no rate for it. `clmad` is not a knob here: the instruction is a
CUDA 13.3 / PTX 9.3 carryless op bought for a Blackwell pipe balance,
and `aws/build.sh` already drops it the moment a pre-Blackwell arch is
in `ARCHES`. The T4 client is the shipping packed arithmetic and
storage preset with `PACKED_CLMAD=0`, compiled for `sm_75` only.

## Boundaries, written before the run

- **Floor:** one-addition-per-step walk. On the 188-SM RTX PRO 6000 that
  is ≈ 23–25 B/s. A 40-SM T4 cannot cross a number that scales with SM
  count. The figure that ranks the part is B/s on the T4 itself, and
  iterations per dollar against g6/g6e/g7/g7e.
- **Reference:** RTX PRO 6000 shipping 15.115792 B/s
  ([benchmarks/top-clmad](benchmarks/top-clmad/summary.json)). Scaled by
  SM count alone (40/188) that is 3.22 B/s; Ada's per-SM issue is
  predicted at half of Blackwell ([ADA-L4-L40S.md](ADA-L4-L40S.md)), and
  Turing is not faster than Ada, so the honest prior is **well under
  2 B/s**.
- **Unit:** billions of completed scalar updates per second.
- **Class:** engineering. Same walk, same 5.3125 products/update.

## Acceptance

A T4 receipt is a verified median on a GPU whose name is T4, identity
matching `PACKED_CLMAD=0`, automatic workers. It does not have to beat
18 B/s; it has to be the fastest client this tree can emit for
`sm_75`. Inadmissible: filing a T4 number under Ada, enabling `clmad`
without a compile-and-run receipt, or quoting the 6000 rate as if it
transferred.

## Recipe

On the instance:

```sh
make bench-g4dn
```

Through Modal:

```sh
make bench-g4dn-modal
```

Thin binary, no GPU needed to compile:

```sh
make gpu-g4dn
# or, for the campaign bucket:
# ARCHES=75 CLMAD=0 ./aws/build.sh /path/to/ecc2k130
```
