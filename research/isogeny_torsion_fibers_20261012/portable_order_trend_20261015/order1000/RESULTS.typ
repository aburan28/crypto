#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)

= Order 1,201: complete recovery and counted trend

The frozen first eligible different-j cell is $(p,r,h)=(1237,1201,1275)$.
The independent full-fiber audit checked 1,531,275 source points,
510,425 three-point image fibers and all 1,442,401 ordered subgroup
pair identities. The seven-arm batch then completed 84/84 cold
children, 1,785/1,785 scalar replays and 12,486 independently
checked relation witnesses.

$S_("group")=("adds"+"doubles")/(N sqrt(r))$ includes common
validation, setup, every target query and final verification. It is
the same counted unit for all seven arms and a *lower bound* on complete
portable S. At $N=64$, medians of three rotated rounds:

#table(
  columns: (2.2fr, 1fr), stroke: 0.5pt + rgb("cbd9de"), inset: 4pt,
  [*Arm*], [*Median $S_("group")$ per target*],
  [Source projected relations], [153.842],
  [Isogeny-image relations], [162.213],
  [Destination-direct relations], [139.075],
  [Full rational-image relations], [1,532.302],
  [Reusable source BSGS], [1.315],
  [Source negation rho], [4.071],
  [Rho A/A repeat], [4.071],
)

Paired full-image/rho median group ratio: 376.364; against BSGS:
1,165.090. Full image exceeds both in group counts, Fp
multiplications and inversions in every prefix and round. Fp axes
overlap group calls and are not added to the table.

#image("TREND_FIVE_V2.svg", width: 100%)

#pagebreak()
= Relation support and remaining gates

The five-order derivation rechecks 420 children and 8,925 scalar
outputs. Orders 229, 421, 433, 457 and 1,201 span 5.24x, but the
field prime/cofactor vary and three cells cluster near 430. The
fitted exponent is unset. At order 1,201 the full-image/destination
support was 122/128 keys. Full image returned 2,316 relations in
23,413 trials; direct destination returned 2,328 in 21,653.
Full image scans all 1,531,275 rational source points in each cold
child. Its largest N64 peak RSS was 206,864,384 bytes below 512 MiB.

Shared-host N64 median CPU was 2.724 seconds for full image,
3.784 milliseconds for rho and 0.954 milliseconds for BSGS. These
are L0 descriptions, not admitted isolated wall results.

Complete portable S still requires disjoint native counts and pinned
cross-host conversion for field work outside group calls,
square-root branches, hash/table probes, row algebra, allocation and
I/O. Physical-core isolated wall with A/A, independent-host
performance replay, and a broader order trend remain separate.
The direct chi-coordinate resultant/root oracle is the structural
next test: match exhaustive decompositions and recover the same
targets without the full rational-source census.

Source-linked evidence and exact commands: `RESULTS.md`,
`BATCH_RUN.md`, `RUN_SPEC_V2.json`, `TREND_FIVE_V2.json`, the run
JSONL and independent replay JSON. The failed V1 target freeze and
the V1 figure with its overlap are preserved.
