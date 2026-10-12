# n83 six-summand one-prefix pruning gate

Registered before code changes and execution. This is a focused follow-on
to [the direct-system result](../f6_n83_sixsum_wide_20261005/RESULT.md),
stacked on PR #1411 at commit
`7e760719eb92cbf7603932d5e3b39a3a99f5924c`. That experiment
admitted a 428-variable, 415-equation n83 K0 six-summand chain, but its
ordinary T001 own-degree root reduction produced only
`v80 ⊕ v345 ⊕ v426 = 0`: one source bit can be absorbed by two free
intermediate bits. This gate tests whether fixing a *whole valid first
source summand* makes the same reduction reject or restrict other
summands. It is not a relation-collection or speedup experiment.

Freeze curve `icv1-f2m83-tm6151469093347-debefd74`, standard
dimension-16 source subspace (64,907 geometric points), public
subgroup target T001
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`, source preimage
`U=[4^{-1} mod r]T001` with rational 4-torsion offset 0, and source
point indices
`[0,2,4,6,8,10,64,128,512,1024,4096,8192,16384,32768,49152,64906]`.
These 16 deterministic entries are a finite probe panel, not a random
sample or natural relation-yield estimator. Record each point's exact
abscissa and flag any duplicate abscissae.

Implement Boolean constant substitution of the first 16 coordinate
bits into the complete direct chain; independently evaluate the
substituted planted system on the known six-point witness with first
source index 0. For each frozen ordinary prefix, run one own-degree
echelon reduction only. Record term occurrences, distinct columns,
rank, contradiction, exact affine rows, rows involving only the
remaining source-coordinate variables, peak RSS, elapsed bounded time,
and any failure. The structural success criterion is at least one
sound contradiction or at least one new affine row involving only
remaining source bits. A row containing free intermediate bits does
not meet it. Any contradiction must be checked against the planted
control logic, and no candidate witness may count without exact full
curve/group replay.

Run one native worker, with a 120-second process time guard and sampled
7-GiB RSS termination as in the prior gate. A memory/time/column cap is
an inconclusive row, not a refutation. If none of the 16 branches meets
the structural criterion, stop this one-prefix reduction strategy under
the frozen panel; do not extrapolate a global pruning rate. The
combinatorial capacity ceiling, full IC cost, candidate ID, one-target
online interval and matched rho remain unmeasured. No isolated CPU
partition is available; report timings only as exploratory resource
diagnostics and no wall-time speedup.
