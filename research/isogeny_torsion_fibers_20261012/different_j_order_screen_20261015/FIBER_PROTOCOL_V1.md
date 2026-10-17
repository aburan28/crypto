# Frozen full rational-fiber audit of the selected order-229 pair

**Status: preregistered before implementation and execution.** The
separately frozen [20-prime screen](RUN_SPEC_V1.json) selected
`p=211`, `r=229`, `F_(p²)=F_p[s]/(s²−2)`, source
`y²=x³−21x+37`, destination `y²=x³−81x−243`, and
`T=(3,1)` for the rational order-three kernel. The independent
[direct point count](CENSUS_V1.jsonl) found 44,655 points on both
curves. This audit tests the full map and fibers; it does not run a DLP
solver or measure a method speedup.

Enumerate every source point over `F_(211²)` by constructing the full
square-to-roots table and visiting every field `x`; include the point
at infinity. Require exactly 44,655 distinct valid source points.
Compute the explicit degree-three Vélu image of each point, reject any
off-destination image, sort source points by the public tuple
`(infinity,x0,x1,y0,y1)`, and SHA-256 the line-by-line source/image
tuple stream in that order. Require exactly 14,885 distinct image
points, each with three source preimages; the identity fiber must be
`{O,T,-T}`. Record the source/image digest, count and fiber histogram.

Select the first point in that sorted source list whose multiplication
by the fixed cofactor 195 yields a nonidentity point `G` with `[229]G=O`.
Verify its image has exact order 229, enumerate all 229 source and image
subgroup elements, require both sets to be distinct, and check
`φ(P+Q)=φ(P)+φ(Q)` for all 229² ordered subgroup pairs. SHA-256 the
ordered source/image subgroup tuple stream and record both generators,
the digest and all subgroup checks. A separate Rust checker must
recompute the full enumeration, map, fiber counts, subgroup selection
and digests with independent field and point-add formulas. It must
receive only the public frozen pair and the first runner's compact
receipt, never a precomputed map table.

Build both release binaries, then commit their source and binary
SHA-256 hashes, exact command/script, the preceding screen/census
input hashes, 120-second per-process wall and 512 MiB RSS ceilings
before execution. Run on the shared host with no output overwrite;
retain stdout, stderr and exit codes. An exit, timeout, digest mismatch
or resource breach censors the corresponding audit, not a negative
mathematical result. A passing finite audit establishes these counts
and explicit map checks only for this selected field and models. It
leaves relation yield, target recovery, cold counted `S`, isolation,
second-host performance replay and a general family theorem open.
