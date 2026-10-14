# N83 four-byte exact root index capacity and lookup gate

The retained N83 S3 state is the 11-byte first root from the preceding
experiment. This candidate changes only its open-addressing root slots:
each four-byte word stores a nonzero 29-bit state ordinal, one root selector,
and two fingerprint bits. Empty is zero. A fingerprint match is followed by
reconstruction and comparison of the full canonical 83-bit root, so table
lookup remains exact and retains the first inserted witness. The table has
`2^31` slots and an 8-GiB payload. The lower occupancy may reduce target
lookup time even though each slot carries a shorter fingerprint.

The exact symmetric state count is
`S(K)=83*K*(K-1)/2+42*K`. At K3550, `S=523,005,525`, so the complete
two-root upper bound fits `2^31` slots and all ordinals fit 29 bits.
The 11-byte state plus four-byte table payload is below 14 GiB. This is a
payload bound; the observed process peak and 14-GiB guard decide admission.
`RootFour.lean` checks the bound and slot encode/decode arithmetic.

Run six preselected arms on the same previously used public N83 point:

| Order | Arm | Purpose | Child cap |
| ---: | --- | --- | ---: |
| 1 | K3400 root4 | lower-load four-byte table | 1500 s |
| 2 | K3400 root5 | identical-base five-byte control | 1500 s |
| 3 | K3500 root4 | larger factor base | 1500 s |
| 4 | K3550 root4 | observed memory boundary | 1500 s |
| 5 | K2543 root4 | low-load size control | 900 s |
| 6 | K2543 root5 | identical-base five-byte control | 900 s |

`root4` and `root5` both use the same 11-byte state, symmetric S3 traversal,
rank seed 29, four rank threads and one target thread. Every arm runs in a
fresh process with a 14-GiB observed RSS guard. A same-point strong-rho arm
uses the frozen `fold-L64-DP6-E0` configuration. The cloud verifier replays
the complete factor base, subgroup points, four-point target relation,
scalar, exclusive online phases and source layout. Paired controls require
equal base hash, exact root-key count, scalar and target probes. Preserve
timeouts, guard stops and failures with their raw receipts.

These are single-point correctness, capacity and exploratory timing runs.
The selection comparison and fresh-point Tier 1 paired study remain separate
experiments. The root4 arms enter that selection only after this one-point
gate supplies a verified candidate manifest and complete resource receipt.
