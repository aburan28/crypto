"""Render committed numbers; the scoreboard itself never computes results."""
import json
from collections import defaultdict
from pathlib import Path

HERE=Path(__file__).resolve().parent
START='<!-- BEGIN toy-quotient-controls-20260924 -->'
END='<!-- END toy-quotient-controls-20260924 -->'


def render():
    s=json.loads((HERE/'results/run-002/summary.json').read_text())
    groups=defaultdict(list)
    for r in s['rows']:
        groups[r['split'],r['summands'],r['variant']].append(r)
    totals={key:sum(r['total_opcodes'] for r in rows) for key,rows in groups.items()}
    lines=[]
    htmlrows=[]
    for (split,m,v),rows in sorted(groups.items()):
        total=totals[split,m,v];ratio=total/totals[split,m,'ordered'];verified=sum(r['verified'] for r in rows);count=sum(r['cases'] for r in rows)
        lines.append(f'| {split} | {m} | {v} | {total} | {ratio:.6f} | {verified}/{count} |')
        htmlrows.append(f'<tr><td>{split}</td><td>{m}</td><td>{v}</td><td>{total}</td><td>{ratio:.6f}</td><td>{verified}/{count}</td><td>accounting</td></tr>')
    fresh=[r for r in s['rows'] if r['variant']=='invariant' and r['split'].startswith('fresh')]
    bounds={m:(min(r['opcode_ratio_to_ordered'] for r in fresh if r['summands']==m),max(r['opcode_ratio_to_ordered'] for r in fresh if r['summands']==m)) for m in (2,3)}
    scale=s['paired_controls']['scale']; inv=s['neighbor_degree_pairs_by_variant']['invariant']
    verdict=f'''The finite quotient helps the three-summand Python reference, but regresses for
two summands. On fresh supports, invariant/ordered complete-call opcode ratios
are {bounds[3][0]:.3f}–{bounds[3][1]:.3f} for three summands and
{bounds[2][0]:.3f}–{bounds[2][1]:.3f} for two summands across the five models.
The predeclared requirement to improve every model and both summand counts by
at least 20% is **{str(s['frozen_success']).lower()}**. The three-summand result
is a bounded Python opcode-proxy improvement, not a native runtime or full-DLP
speedup. Classification: accounting; full-DLP S and rho/floor ratios remain null.
'''
    intro=f'''# Results: quotient benefit depends on decomposition size

{verdict}
All {s['verified']}/{s['systems']} systems preserve the exact unordered signed-point
decomposition set. Three uninstrumented repetitions agree, and an additional
instrumented replay agrees on each input. Independent Buchberger audits verify
{s['independent_verified']}/{s['independent_audits']} selected systems. Actual
signature-based SymPy F5B completes {s['f5b_verified']}/{s['independent_audits']};
{s['f5b_budget_exceeded']} exhaust the declared work budget. Those are explicitly
incomplete F5B runs, not successful F5B validations.

## Representation effects

Coordinate rescaling alone lowers the matrix ideal-completion budget in
{scale['lower']} cases, leaves it equal in {scale['equal']}, and raises it in
{scale['higher']}. Thus degree changes of this kind are not unique to isogeny
descent. Selector relabeling and equation reversal preserve every completion
degree while changing cost, as expected from their unchanged degree filtration.

Within the invariant representation, neighbor/source degrees are lower in
{inv['lower']}, equal in {inv['equal']}, and higher in {inv['higher']} comparisons.
This finite-image encoding does not establish that descent makes the algebra
intrinsically easier. Its rank labeling and exhaustive preprocessing remain
representation choices, fully described in README.md.

## Complete-call proxy ledger

Unit: executed CPython opcode dispatches through setup, encoding, matrix solve
and root extraction, reconstruction, and point verification. Each row aggregates
40 inputs (five curve models, eight targets). Ratios use ordered on the same
split and summand count. This excludes C-internal instructions, arithmetic-table
and fixture generation, and independent audit work; no calibrated machine-operation
total is available. Per-model totals and exclusive phase counts are in the raw
and summary artifacts. Every row is classified as accounting.

| Split | Summands | Variant | Python opcodes | Cost/ordered | Verified |
|---|---:|---|---:|---:|---:|
'''
    intro+='\n'.join(lines)+'\n\n## Fresh quotient comparisons by model\n\n| Split | Model | Summands | Opcode cost/ordered |\n|---|---|---:|---:|\n'
    intro+='\n'.join(f"| {r['split']} | {r['model']} | {r['summands']} | {r['opcode_ratio_to_ordered']:.6f} |" for r in fresh)
    intro+='''

## Decision and reproducibility

Retain the finite-quotient implementation as a tested toy representation with
an explicit setup/lifting cost. Reject the blanket improvement claim. A future
three-summand-only hypothesis would require a new predeclared protocol and new
holdouts; the present two-summand regressions must remain visible.

- `results/run-002/raw.json.gz`: canonical raw inputs, phase counts, per-degree
  traces, repetitions, exact-set hashes, independent bases and budget receipts.
- `results/run-002/summary.json`: full per-model aggregates and decisions.
- `results/run-001/`: superseded audit instrumentation run, retained with its
  original independent-audit source and contract. It is not accepted F5B evidence.
- `contract.json`, `protocol-v1.json`, `protocol-v2.json`: protocol and revision history.
- `dependency-replay.json`: exact replay receipt for the preceding full corpus.
- `verify_replay.py`: hashes, exclusive accounting, correctness and exact replay gates.

The independent observer reports only degrees seen at signature-reduction
boundaries and final reduced-basis degree. Neither these observations nor the
Boolean matrix ideal-completion budget are intrinsic degree of regularity.
'''
    (HERE/'RESULTS.md').write_text(intro)
    section=START+'''\n<section class="panel" id="toy-quotient-controls-20260924">
<div class="panel-head"><h2>Finite symmetry quotient: three-summand proxy benefit, two-summand regression <span class="chip">accounting</span></h2></div>
'''
    section+=f'<p>{s["verified"]}/{s["systems"]} systems preserve exact signed-point multiset coverage. Fresh invariant/ordered Python opcode ratios: {bounds[3][0]:.3f}–{bounds[3][1]:.3f} for three summands; {bounds[2][0]:.3f}–{bounds[2][1]:.3f} for two. The frozen uniform-improvement criterion is {str(s["frozen_success"]).lower()}.</p>\n'
    section+=f'<p>Coordinate rescaling alone changes completion degree: {scale["lower"]} lower, {scale["equal"]} equal, {scale["higher"]} higher. In invariant coordinates, all {inv["equal"]} neighbor/source degree comparisons are equal. These are representation-dependent completion budgets, not intrinsic regularity.</p>\n'
    section+='''<p>Unit: complete decomposition-call CPython opcode dispatches, including setup,
encoding, solving, reconstruction and verification; 40 inputs per row. This is
a Python proxy, not calibrated machine operations. C-internal work, arithmetic-table
and fixture generation, and independent audit costs are outside this metric.
Full-DLP S, rho/floor ratios and speedup are null for every row.</p>
<table><thead><tr><th>Split</th><th>Summands</th><th>Variant</th><th>Python opcodes</th><th>Cost/ordered</th><th>Verified</th><th>Class</th></tr></thead><tbody>
'''+ '\n'.join(htmlrows)+'</tbody></table>\n'
    section+=f'<p>Independent Buchberger audits: {s["independent_verified"]}/{s["independent_audits"]} verified. Signature-based F5B: {s["f5b_verified"]} completed, {s["f5b_budget_exceeded"]} budget-exceeded and retained as incomplete.</p>\n'
    prefix='https://github.com/aburan28/crypto/blob/main/research/toy_quotient_controls_20260924/'
    section+=f'<p>Evidence: <a href="{prefix}results/run-002/raw.json.gz">raw run</a>, <a href="{prefix}results/run-002/summary.json">per-model summary</a>, <a href="{prefix}RESULTS.md">results and limitations</a>. Earlier pilot and neighbor figures above remain unchanged.</p>\n</section>\n'+END+'\n'
    page=HERE.parents[1]/'docs/index-calculus-scoreboard.html'
    content=page.read_text()
    if START in content:
        before,tail=content.split(START,1);_,after=tail.split(END,1)
        content=before+section+after.lstrip('\n')
    else:
        content=content.replace('</body>',section+'</body>')
    page.write_text(content)


if __name__=='__main__':
    render()
