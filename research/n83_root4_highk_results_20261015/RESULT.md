# The high-K root4 panel leaves K3550 and K3596 unresolved

The four-process high-K comparison completed and independently replayed all
16,000 N83 target scalars and four-point relations. The balanced per-target
K3596-minus-K3550 mean was **−0.988686 ms**, with a 95% contiguous
50-target block bootstrap interval of **[−2.080616, +0.104552] ms**. The
frozen selection rule records `BOUNDED_UNKNOWN` for these two bases.

| Fresh process order | K | Mean online ms per point | Setup s | Full process s | Peak RSS GiB | Replay |
| :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| First | 3550 | 26.659630 | 495.309 | 642.422 | 13.813 | `PASS_INDEPENDENT_REPLAY` |
| Second | 3596 | 24.262629 | 457.393 | 594.882 | 13.960 | `PASS_INDEPENDENT_REPLAY` |
| Third | 3596 | 23.808320 | 441.544 | 575.844 | 13.959 | `PASS_INDEPENDENT_REPLAY` |
| Fourth | 3550 | 23.388691 | 428.152 | 558.847 | 13.814 | `PASS_INDEPENDENT_REPLAY` |

K3550's first-minus-second process contrast was **+3.270939 ms**, with a
95% interval of **[+3.017855, +3.521506] ms**. K3596's corresponding
contrast was **+0.454309 ms**, interval **[+0.382508, +0.529026] ms**.
The K3550 process shift exceeds the balanced base difference. All four
arms recorded zero cgroup CPU-throttling and OOM-kill events. K3596's
largest observed peak was 14,989,873,152 bytes, only **42,512,384 bytes**
below the 14 GiB RSS guard; the runner's host exposed four AMD EPYC 9V45
CPUs. Online intervals exclude the target-independent setup shown above.

The common workload contains 4,000 fresh public subgroup points reused in
all four arms. K3550 has 589,300 distinct usable base points and 3,550
folded columns; K3596 has 596,936 and 3,596. Both use the same four-byte
exact root table with 2^31 slots, 11-byte root states, sparse canonical
selector, and no Bloom filter. K3596 is the largest integer K admitted by
the current 29-bit state-index capacity check: the construction's regular
state bound is 536,647,262 at K3596 and 536,945,772 at K3597, while the
encoding requires fewer than 536,870,912 states. This is a bound on the
current root4 implementation's range.

The capacity arithmetic is checked in [Capacity.lean](Capacity.lean), pinned
to Lean 4 v4.34.1. `stateBound_3596` and `stateBound_3597` establish the
exact counts; `max_k_current_layout` proves that K3596 fits and every
integer K greater than 3596 exceeds the 29-bit state-index bound. The proof
assumes the symmetric N83 state-count formula and the producer's strict
29-bit admission rule; those links are the explicit Rust construction and
assertion in the executed source. [LEAN_CHECK.txt](LEAN_CHECK.txt) retains
the build command and `#print axioms` output. The only axioms reported for
the general theorem are Lean's standard `propext` and `Quot.sound`.

The common workload ID is `56376712e771`. The K3550 and K3596 candidate
IDs are
`IC1N83Ckb1a1fb589300PDP4rootRCguidedLAgaussTDdirectISO0h2136bb02d218`
and
`IC1N83Ckb1a1fb596936PDP4rootRCguidedLAgaussTDdirectISO0h3dce48cbfe34`.

The producer executed source commit
`8a37349f3d3cc2bc2960b03167adcf57da83ba77`, with Rust source SHA-256
`c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1`
and executable SHA-256
`318d2310ea4e42d21606159dbaf2a8afed38de8d4e2c7948fc3d969c612b7f71`.
[Actions run 38719753677](https://github.com/aburan28/crypto/actions/runs/38719753677)
passed source preflight, four bounded arms, independent replay, and raw upload.
Artifact `11967981375` has ZIP SHA-256
`1927c1dcbbdc2ccabaa708896624e5e8735359462fc2edfec12464442cbc158b`;
its `analysis.json` SHA-256 is
`a8f7ec556a66212400b629ff9ce85980735b609f0802b9c6e8e3cf2ec811152e`.
The raw ZIP is archived on SSD2 under
`/Volumes/SSD990-2/n83-archive-20261015/raw/`. A separate local verifier
replayed all 16,000 answers and reproduced the analysis digest byte for byte.
`records_index.json` pins the compressed workload and measurements; the
current-head CI workflow regenerates both and requires a clean diff.

The next comparison should put K3550/root4, K3596/root4, and the selected
K3300/root11/wide3 pipeline on one allocation with the same fresh points.
That directly measures their full online target paths before fresh same-point
Tier 1 IC/rho pairs for the selected end-to-end candidate.
