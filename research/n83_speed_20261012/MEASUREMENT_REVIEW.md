# Measurement review for the further rho iteration

## What changed in the measurement approach

The control is the previously selected Fold64/DP6/off rho, remeasured inside the
same executable as the candidate. Fresh points exclude the original target,
K600, all previous optimization targets and each other's signed Frobenius
orbits. Three training and five holdout points are distinct. The full finite
grid, seeds, selection rule and resource envelope are prospective; selections
are written once before holdout. Native code performs generation, monitoring,
selection, statistics and replay. Rendering consumes recorded numbers only.

Five identical A/A pairs quantify local noise before the comparison. The held-out
arms alternate order across 15 matched seed pairs. All 139 solver launches,
including every losing scalar/policy configuration, are retained. There were no
process timeouts or global solve-cap failures. Internal per-lane capped-walk
restarts occur and remain fully charged, with their counters preserved.
Exact Barrett passes correctness but is unselected:
its mean process cost did not improve over reference. This negative candidate
is visible in the report, rather than silently removed from the history.

These changes make the comparison fairer than comparing a tuned candidate with
an untuned default, choosing a favorable seed, or using new rho times against
old IC times. They do not establish exhaustive parameter optimality.

## Cost interval and checks

The primary practical cost is the actual fresh one-target process, including
launch, parsing, target-independent curve/basis/jump preparation, target-dependent
walk, collision recovery, built-in scalar replay and receipt output. The online
clock is secondary. This uses one OS thread; lanes are lockstep walk states,
not additional CPU threads. Each child has CPU3 affinity,12 GiB address-space
limit and60 s cap. Native wait4 polling is nominally5 ms plus scheduling and
poll-body overhead, applied to both arms. A fixed polling correction is not
subtracted. Preparation and phase residuals remain visible in raw receipts.

Independent C++ verification is outside both solver intervals and separately
recorded:0.886574336 s wall for the combined fixture-and-139-answer check. The
native file audit checks commands, sources, binaries, protocol, ordered matrices,
inputs, outputs, clean environments and resource records. These checks do not
reconstruct every intermediate walk state or certify measurement on another
machine. No rank transcript or relation solver is part of this rho-only study.

The first external-verification wrapper failed because /usr/bin/time is absent.
No checker child or answer receipt was produced. Its error is retained. The
direct native monitor then completed successfully. Compile/proof failures and
scratch corruption-rejection controls are retained too; they are not additional
solver launches and do not alter the139-run campaign count.

## Statistical unit and noise interpretation

Five targets are the units for the primary estimate; three seed draws are
nested within each. Each target contributes its three-seed arithmetic cost mean.
The geometric mean of the five paired control/selected ratios is1.069651413981749.
The native20,000-resample target-block percentile interval is
[1.0431783781923283,1.0892379529034781], conditional on selection, exclusions and
this host. Arithmetic mean cost ratio1.0742524986298727 is also shown, rather
than assuming the two estimands are interchangeable. Every target mean improved;
individual seed outcomes and all30 wall times remain visible.

A/A ratio median0.995030985672282, minimum0.9625076546310807 and
maximum1.0465408131500493 describe five same-point/same-seed process pairs.
That observed range is not a confidence interval or a proved noise ceiling.
The primary point estimate lies above its maximum; the primary lower confidence
edge lies inside it. There is evidence of a conditional runtime gain on this
panel, with limited resolution near the lower edge. No noise-adjusted ratio,
independent-host claim or host-wide guarantee is inferred. Failure/cap handling
is explicit and would withhold an unqualified primary win if an arm were incomplete.

Pooled medians/minima and all15 seed-pair ratios are descriptive checks, generated
by the additional native summary in analysis. They do not replace the five-target
estimand, become15 independent target clusters, or select a best seed.
Control/selected pooled wall medians are 2.833963098 s / 2.850071739 s, a slight
reverse ordering. Minima are 1.152991396 s / 0.963512532 s. The 15 paired cost-ratio
median is 1.0757736750707487, minimum 0.9943479875332358 and maximum 1.1973251922223191.
Thus all five target averages improve without every individual launch improving.
All these summaries are shown; no summary is selected after observing its sign.

## What future comparisons should preserve

| Concern | Required approach | Current status |
| --- | --- | --- |
| Strong rho control | Use the best verified curve-specific candidate on separate training, then freeze it | Current fold/reference256/DP6/off is selected within this declared grid |
| End-to-end ranking | Same one-target inputs/resources, whole cost, every completion independently checked | Matched rho comparison complete; no new IC comparison |
| IC/rho comparison | Rerun both on the same new points; include IC preparation and complete online phases | Earlier matched study preserved; new IC/256-rho ratio not established |
| Setup amortization | Report preparation separately; state exact number of targets sharing it | Each rho child is cold; no shared distinguished-point table |
| Tail behavior | Keep arithmetic expected cost, geometric typical ratio and individual targets | Both estimands and all cases retained; five targets limit tail inference |
| Timing variability | A/A controls, alternating arms, native reservation and resource receipts | Complete on one shared VM; physical host neighbors unobservable |
| Common operation unit | Calibrate field/scalar/curve/setup costs and count all phases before S claims | Counters retained; calibrated S remains null |
| Admission and replication | Sealed exact-fixture adapter and independent-host replay for ledger promotion | Not established; this is a scoped engineering runtime record |
| Formal correctness | Match mathematical statements; pin toolchain, build and fresh-check proofs | Four mathematical Lean lemmas pass; Rust word/carry bridges written/native |
| Source-linked reporting | Commit raw evidence, losing candidates, visual source, PDF and benchmark copies | New versioned addendum; earlier60PDF pages preserved |

## Canonical graph review

The benchmark index receives this rho-only runtime entry. The new study plots
the scalar candidate, all nine policy means, A/A noise and every held-out case.
The earlier study-local graphs remain historical records. The canonical
index-calculus scoreboard, IC leaderboard, progress timeline and curve-identity
graphs were checked and left unchanged: no admitted IC/rho ratio, calibrated S,
new curve/model identity, endomorphism or transfer route was produced. Combining
the old1.756671 ratio and the new1.069651 ratio across different panels would be
an unmeasured derived claim and is not done.

The fastest tested setting in this follow-up is an engineering choice within
the declared grid and host. It changes constants and resource tradeoffs; it
does not change the field/basis identity,166-phase quotient or generic reference
scale sqrt(r/166), and it supplies no exponent or deployed-curve conclusion.
