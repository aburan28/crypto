#!/usr/bin/env python3
"""Copy verified frozen figures into the note and canonical static scoreboard."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parent.parent

def main():
    summary=json.loads((HERE/'SUMMARY.json').read_text())
    result=json.loads((HERE/'run_01/results.json').read_text())
    assert summary['all_complete_verified'] and summary['distinct_systems']==192
    assert not summary['confirmation_eligible']
    labels=['Search-only','Degree-3 flat kernel','Degree-3 sparse bucket','Degree-3 hybrid kernel',
        'Selective degree-2 flat','Selective degree-2 one-word','Ordered-specialization search',
        'Fixed-quadratic state','Retained packed state','Full RREF, list','Full RREF, wide',
        'Echelon/affine tail, list','Echelon/affine tail, wide','Packed without diagnostic hashing',
        'Recursive affine/products, list','Recursive affine/products, compact','Direct Gray, scalar',
        'Direct Gray, SIMD','12-variable leaves, scalar','12-variable leaves, SIMD',
        '16-variable leaves, scalar','16-variable leaves, SIMD','Transported blocks, scalar',
        'Transported blocks, SIMD','Initial restriction, list/scalar','Initial restriction, packed/SIMD',
        'Transported 16-variable leaves, scalar','Transported 16-variable leaves, SIMD',
        'Fibers with row elimination','Fibers with scalar column membership','Fibers with SIMD screening',
        'Fibers with zero-only SIMD screen','Direct 16-point, quiet','Full-word 64-point, quiet',
        'Projected bytes, scalar','Projected bytes, SIMD','Leaf16 full words, quiet',
        'Leaf16 projected bytes, scalar','Leaf16 projected bytes, SIMD','Single-stage bytes, quiet',
        'Leaf16 single-stage bytes, quiet','Projected bit planes','Leaf16 projected bit planes',
        'Compiled projected bytes','Compiled full-word 64-point','Leaf16 compiled projected bytes',
        'Compiled full-word 16-point','Fixed compiled dispatcher','Leaf16 compiled full words']
    rows=[r for r in result['costs'] if r['n']==24]
    assert len(rows)==len(labels)==49
    md=[];html=[]
    for row,label in zip(rows,labels):
        values=[f"{row['families'][f]['completion_ns']/1e6:.6f}" for f in ['planted','cross_planted','unplanted']]
        values += [f"{row['display_direct_simd_ratio']:.3f}",'PASS']
        md.append('| '+' | '.join([label]+values)+' |')
        html.append('      <tr>'+''.join('<td>'+x+'</td>' for x in [label]+values)+'</tr>')
    gate_rows=[]
    for arm in json.loads((HERE/'protocol.json').read_text())['new_candidates']:
        g=summary['gate_counts'][arm]
        matched=sum(d['pass'] for d in result['comparisons'] if d['candidate']==arm and d['reference']=='same_policy')
        gate_rows.append(f"| `{arm}` | {g['dramatic_passes']}/18 | {g['incremental_passes']}/18 | {matched}/18 |")
    successes=[]
    for d in result['comparisons']:
        if d['reference']=='retained_frontier' and d['pass']:
            lo,hi=d['ci95_paired_median']
            successes.append(f"| `{d['candidate']}` | {d['split']} | {d['n']} | {d['family']} | {d['paired_ratio_median']:.6f} | [{lo:.6f}, {hi:.6f}] |")
    conclusion='''# Compiled schedules show limited 2x gains; the universal objective is unmet

Every candidate fails the preregistered universal dramatic and incremental gates.
The strongest family is compiled full-word scanning: the 64-point kernel and its
fixed dispatcher each pass **3/18** >2x comparisons. These include two fresh n24
holdout groups. This is a bounded engineering result, not the required universal
gain or a cryptanalytic breakthrough. The primary is ineligible for confirmation.

The complete comparison covers **192 distinct systems, 49 methods and 75,264
observations**. Every result is verified: 147 SAT and 45 UNSAT fixtures per arm.
All 168 predecessor inputs and all 32 prior methods remain, with 24 new holdouts.
No prior figure, exception or frozen verdict is superseded by regrouping.

## Fixed acceptance results

The reference has 39 methods, including quiet, full-word, scalar-projection and
single-stage controls. A candidate removes only itself when it belongs to that
roster. All eighteen group comparisons and every completion are required. The
matched-control threshold is a separate 1.05 lower bound; it cannot replace the
strongest-reference gate. Every universal dramatic and incremental result is
REJECTED, even when some individual groups pass.

| Candidate | >2x strongest-reference groups | >1x strongest-reference groups | >1.05x matched-control groups |
|---|---:|---:|---:|
'''+ '\n'.join(gate_rows)+'''

The complete >2x subgroup results are retained below; none is presented as a pass
of the full protocol:

| Candidate | Split | Variables | Family | Paired median ratio | 95% interval |
|---|---|---:|---|---:|---|
'''+ '\n'.join(successes)+'''

The dispatcher also misses four incremental groups: regression n16 planted and
n20 cross-planted, and holdout n16 planted and cross-planted. Their lower bounds
are respectively 1.000000, 0.999973, 0.925926 and 0.779613. The previous n24
cross-planted seed2097153 exception remains at ratio **0.194210** for the dispatcher
and **0.093432** for direct SIMD byte filtering against the current reference.

## Projection work and measurement order

Across one deterministic solve per fixture, two-stage direct SIMD screening
processes **606,446,080 partial points**, makes **2,368,803 second-stage checks**,
and requires **9,516 complete syndrome checks**. Its single-stage control makes
2,368,803 complete checks on the same ordered points. This large work-count
reduction does not establish a complete-cost gain: setup, coefficient maintenance,
mask extraction and all failed filtering remain charged. These are single-grid
counts, not all eight timing repetitions and not collected independent relations.

The randomized/reverse protocol addresses timing differences noticed between
equivalent wrappers during cyclic discovery. Under the primary protocol, the
median component/dispatcher ratios across fixture medians are **1.002699** at
n16, **0.998049** at n20 and **1.001462** at n24. Their mathematical work and first
models match exactly. These results are consistent with measurement-context or
code-layout effects in the discovery differences; they do not prove a specific
hardware cause. No instruction-cache-cold claim is made.

## Cold complete standalone solve costs

Units are milliseconds including encoding, every image and schedule, filtering,
complete checks, recovery, destruction and result validation. The table pools
two fresh n24 holdout seeds and eight repetitions. Its ratio is a descriptive
ratio of pooled planted costs relative to retained direct SIMD. Acceptance uses
paired pointwise fastest-reference costs, a different statistic. Every row is
class **engineering** and has correctness PASS.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Direct SIMD / method, planted | Correctness |
|---|---:|---:|---:|---:|---|
'''+ '\n'.join(md)+'''

[run_01/results.json](run_01/results.json) retains every interval, model, work count
and exclusive phase cost. [SUMMARY.json](SUMMARY.json) adds individual exceptions,
projection totals and same-kernel alias ratios. Confidence intervals concern
these fixed fixtures, not a population-wide or worst-case guarantee.

## Validation and custody

The producer passes **60 Rust tests**. Thirty-two Python evidence tests pass,
including exact frozen replay, original-equation checks, retained inputs and
methods, compiled-schedule equivalence, random/reverse ordering, null traces,
filter conservation, source/counter corruption and censored costs. The mathematical
contract is in [README.md](README.md), with the producer argument in
[CORRECTNESS.md](CORRECTNESS.md). These are not independent external review.

The first local evidence-test attempt could not copy fixtures because the system
temporary volume was full. The test harness now creates scratch copies on the
checkout volume and registers cleanup before copying. The rerun passed; no timed
source or frozen measurement changed. [VALIDATION_ATTEMPTS.json](VALIDATION_ATTEMPTS.json)
retains both observed outcomes and explicitly leaves unavailable execution
timestamps null.

The complete campaign took **1,556.231 seconds**. Whole-worker peak RSS was
**29,851,648 bytes**, including all methods and reference preparation. Candidate-
specific allocation and calibrated operation counts remain unmeasured. No measured
cell is censored, and no confirmation was launched. [RUN_LEDGER.json](RUN_LEDGER.json)
binds the seven discovery bundles, completed primary, source lineage, report and
validation attempts. Executed artifacts remain immutable. Production, full
index-calculus and rho costs remain **null**. The broad dramatic-gain goal is open.

## Next bounded question

The most useful next question is whether direct construction of the compiled
full-word schedule can materially lower complete cost. The current timing fields
do not separate schedule construction from scanning. A discovery-only cost split
must first establish whether a construction-only change could possibly meet the
strongest-reference target; a counting argument alone is insufficient.
[NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) records that prospective cost-bound gate
and the coefficient-cancellation contract. No implementation or gain is claimed.
'''
    (HERE/'CONCLUSION.md').write_text(conclusion)
    base='https://github.com/aburan28/crypto/blob/main/research/boolean_byte_sieve_20260923/'
    section='''<section class="panel" id="boolean-byte-sieve-20260923">
  <div class="panel-head">
    <h2>Compiled Gray schedules: limited 2&times; groups, universal gain rejected <span class="chip">engineering</span></h2>
    <p>The compiled 64-point kernel and fixed dispatcher each pass 3/18 dramatic
      comparisons, including two fresh n24 holdout groups. Every universal dramatic
      and incremental gate remains rejected. All 192 systems complete with verified
      results; no treatment is eligible for confirmation.</p>
    <p>Sources: <a href="'''+base+'''CONCLUSION.md">decision and all subgroup claims</a>,
      <a href="'''+base+'''CORRECTNESS.md">projection, recovery and schedule arguments</a>,
      <a href="'''+base+'''run_01/results.json">frozen costs and intervals</a>,
      <a href="'''+base+'''SUMMARY.json">case ratios and alias diagnostics</a>,
      <a href="'''+base+'''RUN_LEDGER.json">immutable run ledger</a>.</p>
  </div>
  <div class="table-scroll"><table>
    <caption>Cold complete standalone solve plus validation milliseconds on two fresh n24
      holdout seeds and eight random/reverse repetitions. Display ratios use pooled planted
      costs relative to retained direct SIMD; acceptance uses paired fastest-reference costs.
      All 49 methods are shown. Calibrated IC-operation and rho comparisons are absent.</caption>
    <thead><tr><th>Method</th><th>Planted (ms)</th><th>Cross-planted (ms)</th><th>Unplanted (ms)</th><th>Direct SIMD / method</th><th>Correctness</th></tr></thead>
    <tbody>
'''+ '\n'.join(html)+'''
    </tbody>
  </table></div>
  <p>The run contains 75,264 observations: 147 SAT and 45 UNSAT fixtures per arm.
    Two-stage SIMD filtering makes 9,516 complete checks across 606,446,080 partial
    points in one solve per fixture; reduced checks do not remove its other costs.
    The prior seed2097153 exception remains at dispatcher ratio 0.194210.</p>
  <p>Equivalent component/dispatcher median ratios are 1.002699 at n16, 0.998049 at
    n20 and 1.001462 at n24 under random/reverse orders. No mathematical gain is
    attributed to alias timing. Validation includes 60 Rust and 32 Python evidence
    tests. Whole-worker peak RSS is 29,851,648 bytes; calibrated operations, production,
    full-IC and rho costs remain unmeasured.</p>
  <p>The evidence-test scratch copies were moved to the checkout volume after system
    temporary storage was exhausted; frozen measurements were unchanged. Direct
    schedule construction is a prospective experiment requiring a measured cost split
    first. Historical figures and verdicts remain below.</p>
</section>
'''
    page=REPO/'docs/index-calculus-scoreboard.html';text=page.read_text()
    start='<section class="panel" id="boolean-byte-sieve-20260923">'
    marker='<section class="panel" id="boolean-linear-fibers-20260923">'
    if start in text:
        a=text.index(start);b=text.index(marker,a);text=text[:a]+section+text[b:]
    else:
        assert text.count(marker)==1;text=text.replace(marker,section+marker)
    old='''Partial byte-syndrome filtering is the next
    documented representation hypothesis; its implementation and complete cost remain
    unmeasured. Prior measurements and rejections remain below.</p>'''
    new='''Partial byte-syndrome filtering was the next documented
    hypothesis at that run's close. Its <a href="#boolean-byte-sieve-20260923">subsequent
    implementation and measurement above</a> retain the historical rejection.
    Prior measurements and rejections remain below.</p>'''
    assert old in text or new in text
    page.write_text(text.replace(old,new))
    print('Wrote conclusion and canonical scoreboard with 49 measured rows.')

if __name__=='__main__':main()
