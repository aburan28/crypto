# What the follow-up establishes

All 373 attempts completed: 336 verified recoveries, 36 exhausted searches and 1 retained legacy-width error. Of the 96 attempts in the two rank panels, 56 are IC recoveries and 40 are paired rho/A-A controls.

Both IC arms recovered all four held-out degree 37 public points, eight IC executions in total. Increasing the base from four to eight orbit columns changed mean counted S from 8.566326 to 2.618529; the four/eight count quotient is 3.271426. The eight-column/rho count quotient remains 9.012705. These are incomplete-cost diagnostics, not speedups.

The exact E0 degree 83 method now executes ordinary IC queries at both base sizes, but all four runs reached 256 queries without a point relation or rank gain. This discharges an implementation refusal; it does not pass the N83 target-recovery gate.

The zero-hit N83 result is weak evidence about mathematics: under the existing uniform-query multiset-counting model, these tiny bases cover a vanishing fraction of the subgroup. The model is not a calibrated probability for the fixed hash stream. A scalable algebraic base/witness improvement is still required.

Prime base balance changes across models: icv1-fp16-t295-8d3c3165: abs32-mitm2; icv1-fp18-tm337-d28d5e09: abs32-mitm2; icv1-fp20-t727-cd198a38: abs64-mitm2; icv1-fp22-tm1385-475dcb5f: abs64-mitm2. These are descriptive minima on this panel, not a validated tuning rule or exponent.

Next discriminating work: price the new method's native operations and compare the column-count choice on further unseen points; retain all relative phases in any proposed algebraic correspondence; demand growing base dimension and successful exact E0 N83 recovery before any N131 inference. Recent main also has versioned full-price methods; their host-specific calibration is a follow-up input, not a price silently imported into this AMD run.

## Prime balance on fresh targets

| Exact curve | Descriptive minimum | Mean counted S | Counted / rho | Targets |
|---|---|---:|---:|---:|
|icv1-fp16-t295-8d3c3165|abs32-mitm2|10.342938|4.470230|8|
|icv1-fp18-tm337-d28d5e09|abs32-mitm2|10.506144|6.886933|8|
|icv1-fp20-t727-cd198a38|abs64-mitm2|9.743966|6.635641|8|
|icv1-fp22-tm1385-475dcb5f|abs64-mitm2|10.811440|11.527583|8|

## Direct versus folded tables, identical bases

| Exact curve | Direct S | Folded S | Direct/folded counted quotient | Folded/rho counted quotient |
|---|---:|---:|---:|---:|
|icv1-f2m13-t181-515ee569|127.717893|14.786094|8.637703|1.031459|
|icv1-f2m19-t797-b6cf2467|33.991107|3.541676|9.597463|1.099422|
|icv1-f2m23-t5197-69e76b73|13.711346|2.463207|5.566461|2.164225|
|icv1-f2m31-tm90707-c95f16f5|26.775227|1.866643|14.344050|1.560362|

## N83 coverage model, separated from observations

The already-used multiset count C(B+2,3)/(r-1) bounds the number of reachable sums per uniformly randomized nonzero a, conditional on b. Distinct tuples can collide. The following union bounds are model calculations, not measured rates or calibrated probabilities for the frozen deterministic hash stream.

| Columns | Signed points | Queries across two attempts | Checked relations | 256-query model union bound |
|---|---:|---:|---:|---:|
|16|2656|512|0|3.31004325e-13|
|4|664|512|0|5.18947017e-15|

Sources: [sealed-record summary](summary.json), [per-target rows](per-target.json), [counting-ceiling implementation](../../src/cryptanalysis/ic_boundary.rs), [Galbraith et al., On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf), [Caminata and Gorla, Solving degree, last fall degree, and related invariants](https://arxiv.org/abs/2112.05579), [recent main full-price implementation](https://github.com/aburan28/crypto/blob/fd8bcacf4cc097cc0dcb55ba068ad94a3d3202ae/src/cryptanalysis/ecbench/full_price.rs), and the [parent29-reference literature review](../ecbench_open_questions_20261014/literature/REPORT.md). No new mathematical theorem is claimed; the counting formula and its assumptions are inherited, and no fixed-stream probability bound is asserted.
