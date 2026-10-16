#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)

= Four-order complete-solve panel

*Status:* 252/252 new cold children and 5,355/5,355 new scalar outputs
passed independent Rust replay. Across the checked order-229 reference
and three new prime-to-characteristic cells, 336 children and 7,140
scalar outputs are in the native count derivation. New witnesses:
37,061; all returned relation rows and scalar identities checked.
The run spec and all source/binary/input digests were committed before
measurement. Host: Darwin arm64, shared L0.

$S_("group") = ("adds" + "doubles") / (N sqrt(r))$ is the cold
counted group-operation total per target, with N=64 here. It includes
common/setup, all queries and final verification. It is a *lower bound
on complete portable S*: field work outside group calls and
hash/table/row/alloc/process work have no pinned common conversion.
The direct source negation-rho reference is on identical public targets;
source BSGS is the memory control. The two-class generic collision
floor is $sqrt(pi / 4) approx 0.8862$ per target.

#table(
  columns: (2.0fr, 0.8fr, 0.8fr, 0.8fr, 0.8fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 3.5pt,
  [*Arm*], [*$r=229$*], [*$r=421$*], [*$r=433$*], [*$r=457$*],
  [Source projected], [58.28], [80.04], [78.84], [64.15],
  [Isogeny-image], [64.56], [86.17], [85.39], [69.73],
  [Destination-direct], [69.56], [81.83], [68.04], [72.28],
  [Full rational image], [148.49], [405.89], [332.16], [336.69],
  [Reusable source BSGS], [1.76], [1.60], [1.52], [1.54],
  [Source negation rho], [6.64], [5.59], [5.14], [5.31],
  [Rho A/A], [6.64], [5.59], [5.14], [5.31],
)

Medians of three complete cold children per cell. Paired median
full-image/rho group-count ratios are 22.355, 72.401, 64.275 and
63.389. In every cell, prefix and round, full-image exceeds BSGS
and rho in each counted arithmetic axis: group work, Fp multiplications
and Fp inversions. Native categories remain unpriced.

The frozen native missing-price sensitivity rechecks 48 matched
full-image/rho child pairs. At N64, the median extra *omitted rho*
charge per observed walk step needed to tie the full-image group-count
bound is 112.973, 303.573, 334.844 and 288.595 GAE across the four
orders, setting all omitted full-image work to zero. These are
thresholds, not measured native-operation prices or complete portable
S. See `sensitivity/RESULTS.pdf` and `SENSITIVITY_V1.json`.

#image("TREND_V1.svg", width: 100%)

#pagebreak()
= Relation geometry and remaining gates

The full rational-image table changes two-sum support but does not
consistently increase it. The larger three cells have full-image
support of 108, 106 and 112 keys versus destination-direct 112,
120 and 122. Each relation arm reached rank nine and recovered its
requested scalar; the independent checker replayed 37,061 new
relation witnesses. The full-image construction enumerates the
entire rational source in every cold child, so its setup grows with
the field. Its counted cost is 332–406 per target on the new cells,
against rho's 5.14–5.59. The four orders span only a factor of two
and use different cofactors; a stable exponent fit is not admitted.

#table(
  columns: (0.7fr, 1.1fr, 1.2fr, 1.2fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 4pt,
  [*$r$*], [*full/direct keys*], [*full relations/trials*], [*direct relations/trials*],
  [229], [98 / 90], [2,329 / 5,389], [2,628 / 6,776],
  [421], [108 / 112], [2,457 / 9,475], [2,348 / 8,853],
  [433], [106 / 120], [2,393 / 9,595], [2,234 / 8,113],
  [457], [112 / 122], [2,376 / 9,556], [2,339 / 8,830],
)

The direct chi-coordinate resultant/root oracle is the next
structural test: find useful destination columns without a full
rational-source census, first prove exact decomposition coverage
against independent enumeration, then price complete target solves.
The pairing transfer has its own finite-field DLP and input-search
cost gate. Complete portable pricing, L2-or-better interleaved timing
with A/A and a second-host performance replay remain open. This
run's CPU and wall meters are shared-host L0 descriptions; the
canonical frontier and scoreboard graphs are unchanged.

The isolated-host access check found no configured isolab and a zero
Standard On-Demand EC2 vCPU quota in five checked regions; AWS denied
the quota-increase request with an AccessDeniedException. A frozen
Linux hosted-runner protocol is ready for independent-host replay;
that virtual machine cannot admit L2 wall timing.
The full local release suite was red under sandbox socket denial and
shared-host timeouts; focused reruns passed. See the validation receipt.

Source-linked records, run PDFs, exact hashes and accepted limits:
`RESULTS_V2.md`, `BATCH_RUN_SPEC_V1.json`, and `TREND_V1.json`.
