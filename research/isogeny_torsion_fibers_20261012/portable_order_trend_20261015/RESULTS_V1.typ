#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)

= Portable cost axes and bounded order screen

*V1 evidence, 2026-10-15.* The native analysis re-derived 84 cold order-229
children and 1,785 independently checked scalar outputs. A separate bounded
screen evaluated 30 eligible primes, accepted eight different-$j$ rows and
passed 512 independent map-witness replays. Source and binary digests, exact
commands and raw records are in `FREEZE_V1.json` and `RESULTS_V1.md`.

$S_("group") = ("adds" + "doubles") / (N sqrt(229))$ uses the repository's
counted group-operation convention and includes common work, setup, all
queries and final scalar verification. It remains a *lower bound on complete
portable cost*: the unpriced native work is named in the source report. Field
multiplication and inversion counts are separate, overlapping axes.

#table(
  columns: (2.1fr, 0.7fr, 0.7fr, 0.7fr, 0.7fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 3.5pt,
  [*Arm*], [*N=1*], [*N=4*], [*N=16*], [*N=64*],
  [Source projected], [91.39], [72.00], [61.55], [58.28],
  [Isogeny-image], [101.50], [79.79], [68.15], [64.56],
  [Destination-direct], [77.98], [73.83], [64.84], [69.56],
  [Full rational image], [5981.66], [1534.12], [428.12], [148.49],
  [Reusable source BSGS], [4.10], [2.45], [1.97], [1.76],
  [Source negation rho], [8.06], [6.74], [6.54], [6.64],
  [Rho A/A], [8.06], [6.74], [6.54], [6.64],
)

Values are median $S_("group")$ per target across three cold rounds. At
N=64, paired full-image/rho ratios are 23.402, 21.544 and 22.355; paired
full-image/BSGS ratios are 85.200, 83.618 and 84.147. The full-image arm
exceeded both direct controls in each counted arithmetic axis on every
matched child. The generic two-class collision floor is 0.8862 per target;
the present lower bounds are not full method costs.

#image("PORTABLE_V1.svg", width: 100%)

#pagebreak()
= Larger-order candidates and next gate

The V1 selection picked $r=p$ at p=277, 313 and 349. Those rows passed their
arithmetic checks, but have a different subgroup regime from the order-229
reference. The versioned V2 protocol explicitly selects the first three
distinct prime-to-characteristic orders from the independently checked
candidate rows, before any new target generation or timing.

#table(
  columns: (0.65fr, 0.65fr, 0.85fr, 1.1fr, 1.6fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 4pt,
  [*$p$*], [*$r$*], [*cofactor*], [*$j$ source / dest*], [*panel state*],
  [409], [433], [387], [84 / 368], [V2 cell A],
  [421], [457], [387], [55 / 176], [V2 cell B],
  [457], [421], [495], [119 / 57], [V2 cell C],
  [211], [229], [195], [77 / 160], [checked reference],
)

#image("SCREEN_V1.svg", width: 100%)

The exact next step is a full rational-fiber census, independent subgroup-map
replay and matched 64-target relation/recovery panel at each V2 cell. A
separate L2-or-better timing session and a second-host performance replay
remain required. The broader goal is an end-to-end, cold structural ECDLP
method whose complete portable cost falls below same-target rho across
comparable orders. Direct chi-coordinate roots and pairing transfer remain
separate possible mechanisms, each requiring its own full cost and scalar
replay. The frontier graph is unchanged.
