#set page(paper: "a4", margin: 18mm)
#set text(size: 9.5pt)
= Full rational-image factor base: exact support and local cost
*Run:* 2026-10-15 19:16 UTC · *Host:* Darwin arm64, macOS 26.6, shared L0 · *Input:* public order-61 target over F67², source j=7, destination j=11, degree-three map.

The protocol was committed before measurement. The native constructor enumerated 4,575 source points and 1,525 rational image points, with three preimages per image point. It selected eight projected sign orbits after examining 21 destination candidates. Source, isogeny-image, destination-direct and full-image relation arms reached rank nine in each of three rotated rounds. Together with direct negation rho and rho A/A, all *18/18* arm-round rows recovered scalar *41*. A separate Rust checker verified those 18 scalars, eight source-preimage witnesses and 178 relation witnesses.

Each factor base covered 56 of 60 nonzero subgroup two-sums. Source and projected-image support sets were identical. Full-image versus projected-image support had 52 shared points and four unique points on each side; full-image versus destination-direct had the same 52/4/4 split. Equal support counts therefore concealed different support members. Full-image relation (trials, relations) were (11,10), (23,19) and (11,11); destination-direct had (15,14), (11,11) and (21,19).

The one comparison unit is microseconds of current-thread CPU, median of three shared-host rounds. Setup is included in total. The timer excludes common JSON reading, model construction and process startup. Rho setup was not separated.

#table(
  columns: (2.2fr, 1fr, 1fr),
  stroke: 0.5pt + rgb("cbd9de"),
  inset: 5pt,
  [*Arm*], [*Total CPU (µs)*], [*Setup CPU (µs)*],
  [Direct negation rho], [58.084], [—],
  [Rho A/A], [56.041], [—],
  [Source projected], [198.333], [110.458],
  [Isogeny-image projected], [258.042], [94.542],
  [Destination direct], [134.584], [81.666],
  [Full rational image], [3962.083], [3842.500],
)

Paired full-image/rho ratios were 73.826, 70.281 and 22.743 (median 70.281). Paired full-image/destination-direct ratios were 31.862, 29.959 and 8.510 (median 29.959). The full-image setup share was 98.872%, 96.982% and 98.532% by round; its median was 98.532%. That setup counted 9,343 group additions, 390,029 field multiplications and 13,819 inversions per round. Peak full-image RSS was 7,077,888 bytes.

#image("GRAPH.svg", width: 90%)

These are complete local timer-interval costs for one target. Lean checks the finite residue sets without axioms; the curve-level point-to-scalar and map proof remains open. Portable counted S, isolated wall timing, increasing-order scaling and independent-host cost replay remain open. A reuse test must charge image construction once and compare direct source BSGS and rho. The Markdown report links the frozen spec, proof map, hashes and raw receipts.
