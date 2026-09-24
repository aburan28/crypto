#!/usr/bin/env python3
"""Render the accounting diagnostic from hash-bound result files."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parent.parent

def main():
    summary=json.loads((HERE/'SUMMARY.json').read_text())
    result=json.loads((HERE/'analysis_01/results.json').read_text())
    census=json.loads((HERE/'quotient_census_01/results.json').read_text())
    assert summary['all_complete_verified'] and summary['cells']==24 and summary['observations']==10560
    cost_md=[];cost_html=[]
    for row in summary['n24_cost_table']:
        values=[row['variant']]+[f"{row['families'][f]/1e6:.6f}" for f in ['planted','cross_planted','unplanted']]+[f"{row['dispatcher_over_arm_planted']:.3f}",'PASS']
        cost_md.append('| '+' | '.join(values)+' |')
        cost_html.append('      <tr>'+''.join('<td>'+v+'</td>' for v in values)+'</tr>')
    bound_rows=[];bound_html=[]
    for row in result['bounds']:
        if row['n']<16:continue
        c=row['optimistic_ceiling_ci95'];ceiling='null' if c is None else f'[{c[0]:.4f}, {c[1]:.4f}]'
        values=[str(row['n']),row['family'],row['kind'],ceiling,row['decision']]
        bound_rows.append('| '+' | '.join(values)+' |')
        bound_html.append('      <tr>'+''.join('<td>'+v+'</td>' for v in values)+'</tr>')
    quotient_rows=[]
    for n in [12,16,20,24]:
        for k in [4,5,6]:
            rows=[r for r in census['cells'] if r['n']==n and r['low_variables']==k]
            ranks=sorted(set(r['annihilated_rank'] for r in rows));dimensions=sorted(set(r['quotient_dimension'] for r in rows))
            quotient_rows.append(f"| {n} | {k} | {', '.join(map(str,ranks))} | {', '.join(map(str,dimensions))} |")
    note='''# Construction-only optimization is not supported by this cost diagnostic

The cost split does not justify implementing a constructor-only candidate for a
universal dramatic claim. Under the measured unchanged-scan assumption, six of
nine relevant 16-point groups and three of nine dispatcher groups have guarded
optimistic ceilings below 2.0. The other groups remain **INCONCLUSIVE**, including
every 64-point group: extracting that scan changed timing beyond the declared
representativeness tolerance. The guards and thresholds were not relaxed.

| Policy | Guarded groups below 2x | Not ruled out | Inconclusive |
|---|---:|---:|---:|
| 16-point scan | 6/9 | 0/9 | 3/9 |
| 64-point scan | 0/9 | 0/9 | 9/9 |
| Fixed dispatcher | 3/9 | 0/9 | 6/9 |

This is an **accounting diagnostic**, not a measured speedup or a general
impossibility theorem. No constructor optimization is implemented. The frozen
protocol's classification string is retained as `engineering diagnostic`; this
report classifies the measurement work as accounting under the repository rule.
The original dramatic-gain objective remains open.

## What was measured

All **10,560 observations over 24 discovery systems and 55 methods** complete with
verified outcomes. Forty-nine retained methods are joined by profiled and
unprofiled rebound arms for the two compiled kernels and their dispatcher.
Profile and rebound share exactly the same non-inlined scan function and builder.
Every new result matches the original policy's model, semantic work and trace.

Encoding, construction, scan and release are separate timed phases. Complete
totals also charge wrapper and result-validation work. Under a fixed scan cost,
even free construction leaves `T_new >= T_scan`, so the maximum ratio is
`T_reference/T_scan`. The optimistic diagnostic subtracts the whole positive paired
profile/rebound total difference and the largest observed clock pair from the scan
timer. Nonpositive estimates yield null ceilings. This is a deliberately favorable
conditional estimate, not an absolute timing-error guarantee.

Both profile/rebound and rebound/original 95% intervals must fit in [0.8,1.25]
before a ceiling is interpreted. The original and extracted scans can differ in
code layout and data placement. A constructor change that also changes scan cost
would be a different mechanism, outside this unchanged-scan bound.

## Relevant conditional ceilings

Intervals below describe the optimistic ratio to the pointwise fastest complete
reference. An interval below 2.0 is usable only when its comparability guards pass.
Inconclusive rows stay inconclusive even when their displayed ceiling is small.
All n12 rows and every phase value remain in the machine-readable result.

| Variables | Family | Policy | Optimistic ceiling, 95% interval | Decision |
|---:|---|---|---|---|
'''+ '\n'.join(bound_rows)+'''

The valid dispatcher groups are all at n20; their upper ceilings are below 0.973.
The n16 comparisons are too sensitive to profiling and extraction to support a
bound. The 64-point kernel's rebound/original intervals exceed the tolerance,
including approximately 1.21–1.33 in several groups. Those observations do not
permit declaring the original 64-point constructor unhelpful.

## Complete costs on n24 discovery fixtures

Units are milliseconds per cold solve plus result validation, pooling the two
discovery seeds and eight random/reverse repetitions. The descriptive ratio uses
the retained dispatcher's planted cost. It is not a promotion statistic. Profiled
methods include their observer work; no setup phase is deducted from these totals.
All methods are shown and correctness is PASS.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Dispatcher / method, planted | Correctness |
|---|---:|---:|---:|---:|---|
'''+ '\n'.join(cost_md)+'''

The campaign took **201.581 seconds**. Whole-worker peak RSS was **19,873,792 bytes**,
including every method, reference preparation and clock calibration. Candidate-
specific memory and calibrated-operation costs remain unmeasured.

## Preserved analysis failure and validation

The original analyzer stopped on a dispatch-arm naming mismatch after every worker
completed. Its source and failure remain in the sealed `run_01` bundle. The additive
`analysis_01` correction binds that input manifest and fixes only the dispatch
name and output location. It changes no timed Rust source, observations, formulas,
guards, thresholds or cohorts. [ANALYSIS_CORRECTION.md](ANALYSIS_CORRECTION.md)
documents the execution and replay paths.

The producer passes **63 Rust tests**. Twenty Python evidence/census tests pass,
including exact replay, original-source retention, source/work corruption, null
ceilings, failed comparability guards, censored outcomes and the preserved original
analysis failure. These are producer checks, not external review. The experiment
uses discovery inputs only and has no holdout or performance promotion.

## A different mathematical degree of freedom

An exact census projects equation-coordinate words along the span of the quadratic
coefficients internal to a chosen low-variable block. The projected system is affine
in that block for every outside assignment, even when the original interaction graph
has edges inside it. Every original solution survives the projection.

The converse is false: projecting `x*y+1` can produce the zero system. A complete
method must enumerate any remaining affine freedom and verify every candidate on
all original equations. The tests explicitly preserve this necessary/sufficient
distinction. The census contains 72 configurations on the same discovery inputs:

| Variables | Low-block size | Observed annihilated ranks | Equation-quotient dimensions |
|---:|---:|---|---|
'''+ '\n'.join(quotient_rows)+'''

These are exact quotient dimensions, not independent-relation counts or predicted
success rates. The complete projected-fiber solver and its cost are unimplemented.
[NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) records the required original-equation
checks, rank-deficient cases, setup accounting and matched controls.

[RUN_LEDGER.json](RUN_LEDGER.json) binds the original execution, corrected analysis,
structural census and report. Production solver, full index-calculus, calibrated-
operation and rho costs remain **null**. Neither this cost bound nor the structural
census establishes an asymptotic or cryptanalytic result.
'''
    (HERE/'CONCLUSION.md').write_text(note)
    base='https://github.com/aburan28/crypto/blob/main/research/boolean_schedule_construction_20260923/'
    section='''<section class="panel" id="boolean-schedule-construction-20260923">
  <div class="panel-head">
    <h2>Construction cost diagnostic: limited bounds, no optimized solver <span class="chip">accounting</span></h2>
    <p>Even free construction misses 2&times; in six guarded 16-point groups and three
      dispatcher groups under the measured unchanged-scan assumption. Every 64-point
      group remains inconclusive because extraction exceeded the declared timing tolerance.
      No constructor optimization or performance promotion is implemented.</p>
    <p>Sources: <a href="'''+base+'''CONCLUSION.md">decisions, phases and limitations</a>,
      <a href="'''+base+'''analysis_01/results.json">corrected diagnostic results</a>,
      <a href="'''+base+'''ANALYSIS_CORRECTION.md">preserved analysis failure</a>,
      <a href="'''+base+'''SUMMARY.json">recorded figures</a>,
      <a href="'''+base+'''RUN_LEDGER.json">artifact ledger</a>.</p>
  </div>
  <div class="table-scroll"><table>
    <caption>Cold solve plus validation milliseconds on two n24 discovery seeds and
      eight random/reverse repetitions. Profiled arms include observer work. Display
      ratios use pooled planted costs relative to the retained dispatcher; no speedup
      or holdout promotion is claimed. All 55 methods are included.</caption>
    <thead><tr><th>Method</th><th>Planted (ms)</th><th>Cross-planted (ms)</th><th>Unplanted (ms)</th><th>Dispatcher / method</th><th>Correctness</th></tr></thead>
    <tbody>
'''+ '\n'.join(cost_html)+'''
    </tbody>
  </table></div>
  <p>All 10,560 observations across 24 systems complete with verified results. The
    original analysis-only naming failure is preserved; an additive correction changes
    neither measurements nor thresholds. Validation includes 63 Rust tests and 20 Python
    evidence/census tests. Whole-worker peak RSS is 19,873,792 bytes.</p>
  <p>Relevant guarded outcomes: 16-point scans 6/9 below the target and 3/9 inconclusive;
    64-point scans 9/9 inconclusive; dispatcher 3/9 below the target and 6/9 inconclusive.
    These conditional timing ceilings do not rule out changing the scan itself.</p>
  <p>An exact discovery census annihilates low-block quadratic coefficient vectors in
    equation space, leaving a necessary affine subsystem. At n24 the quotient dimensions
    are 21&ndash;23 for four low variables, 18&ndash;20 for five, and 14&ndash;17 for six.
    Original-equation checks remain mandatory. Its complete solver and speed are unmeasured;
    production, full-IC, calibrated-operation and rho costs remain null.</p>
</section>
'''
    page=REPO/'docs/index-calculus-scoreboard.html';text=page.read_text()
    start='<section class="panel" id="boolean-schedule-construction-20260923">'
    marker='<section class="panel" id="boolean-byte-sieve-20260923">'
    if start in text:
        a=text.index(start);b=text.index(marker,a);text=text[:a]+section+text[b:]
    else:
        assert text.count(marker)==1;text=text.replace(marker,section+marker)
    old='''Direct
    schedule construction is a prospective experiment requiring a measured cost split
    first. Historical figures and verdicts remain below.</p>'''
    new='''Direct schedule construction was prospective at that run's close. The
    <a href="#boolean-schedule-construction-20260923">subsequent cost diagnostic above</a>
    records guarded bounds and inconclusive cases without implementing the optimization.
    Historical figures and verdicts remain below.</p>'''
    assert old in text or new in text
    page.write_text(text.replace(old,new))
    print('Wrote diagnostic conclusion and 55-method canonical scoreboard panel.')

if __name__=='__main__':main()
