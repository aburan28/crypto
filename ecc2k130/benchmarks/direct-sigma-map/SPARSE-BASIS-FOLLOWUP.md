# Sparse polynomial basis plus direct maps: follow-up assessment

Status: **defer before implementation.**  The direct-sigma panel measured the
shared table at 0.826916x control.  A sparse basis remains a mathematically
coherent compound route, but it is no longer a clear GPU candidate without a
new native/SASS model showing that its reduction cut repays that measured loss
and the additional X-conversion lookups.

## Why the old rejection changes

`THROUGHPUT-30B.md` and `TOP-CLMAD.md` reject a sparse irreducible
pentanomial basis because its polynomial-to-normal conversion is a dense
131-by-131 matrix.  The old implementation model charged a general nibble
table for that conversion while retaining the other basis crossings.

The direct-map work changes the composition.  In a sparse polynomial basis:

- the eight maps `toPolynomial((I + sigma^j) fromPolynomial(p))` remain eight
  131-bit linear maps and occupy the same **56,320 bytes** in the admitted
  three-bit shared layout;
- only X must be mapped to the normal basis on the hot path, for the unchanged
  Hamming-weight branch rule; Y needs the map only on the rare reporting path;
  and
- every polynomial product and square uses the sparse pentanomial reduction.

The dense X conversion can use a six-input-bit table with 22 chunks.  Storing
each row as an aligned low `uint4` and one top word costs
`22 * 64 * 20 = 28,160` bytes and 22 vector plus 22 scalar shared loads.  With
the direct `L_j` table the total dynamic shared allocation is **84,480 bytes**.
The measured device reports 1,024 driver-reserved bytes, so this remains below
the roughly 100 KiB per-SM limit with zero static shared allocation and one
512-thread block (16 warps) per SM.  A seven-bit table would raise the simple
layout to 104,960 bytes before driver reservation and is therefore excluded.

## Static envelope

The repository's SASS accounting prices the current dense-modulus reduction
at about 78 ALU instructions per call and about 6.3 calls per update, roughly
490 slots.  The existing sparse-pentanomial estimate is about 40 per call,
roughly 252 slots, for a provisional **238-slot reduction**.

A six-bit dense X transform performs 110 output-word XORs, 44 basic
shift/mask operations, and six cross-word shift/OR operations: about 160 data
ALU before address arithmetic.  The current structured reduced-input transform
is about 108 operations in the same source model.  The provisional net is
therefore about **186 fewer data-ALU operations per update** before compiler
fusion, at the cost of 44 more shared-load instructions and their bank
conflicts.  These are screening values, not SASS counts or a throughput
prediction.

The resulting hot lookup traffic would be 176 shared instructions for the X/Y
direct maps plus 44 for the X normal map, **220 per update**.  Whether the
reduction cut pays for that traffic is now constrained by measurement: the
56,320-byte direct table fell from 15.893066 to 13.143106 B/s, so the compound
candidate would need a **1.2092x** improvement merely to tie the composed
control while adding another 44 shared instructions.  The provisional
238-slot reduction is not enough evidence of that transfer.  Reduction and
lookup issue on different pipes, but the direct result shows the lookup pipe
is already material.

## Required experiment

Before any implementation or GPU run:

1. choose and record one irreducible modulus, initially
   `x^131 + x^8 + x^3 + x^2 + 1`;
2. generate the sparse reduction, normal conversion and all eight direct maps
   in native C++ from independent field arithmetic;
3. prove every map on all basis vectors, prove reduction by long division, and
   test dense products, squares, inverses, point updates and external normal-
   basis checkpoint round trips;
4. compile native `sm_120` SASS and require a modelled end-to-end gain greater
   than the measured 1.2092x break-even, not merely a reduction-local ALU cut,
   with an 84,480-byte dynamic-shared request, zero static shared, zero walk
   spills and one 512-thread block per SM; and
5. freeze a matched control/candidate replay, corpus, checkpoint, A/A and A/B
   protocol before allocating a GPU.

This is a representation change, not a new iteration function.  If the maps
and checkpoint boundary are exact, seeds, branch decisions, point states and
distinguished-point corpus remain identical.  It does not reduce the number
of carry-less products, so the repository's roughly 22--25 B/s one-addition
pipe floor still applies; this assessment does not claim the 26 B/s objective.
