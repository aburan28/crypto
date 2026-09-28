#!/usr/bin/env python3
"""Render the sealed implicit-minor diagnostic without extrapolating a solver."""
import hashlib
import html
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parent.parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dump(path,value):path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def main():
    run=HERE/'run_01';manifest=json.loads((run/'manifest.json').read_text())
    assert set(manifest['files'])=={p.name for p in run.iterdir() if p.is_file() and p.name!='manifest.json'}
    for name,digest in manifest['files'].items():assert sha(run/name)==digest,name
    meta=json.loads((run/'metadata.json').read_text())
    for name,digest in meta['source_hashes'].items():assert sha(HERE/name)==digest,name
    r=json.loads((run/'results.json').read_text())
    assert r['all_complete'] and r['cells']==12 and r['observations']==5280
    blocked=sum(g['decision']=='BLOCKS_MEASURED_FIXED_PREPROCESSING' for g in r['gates'])
    verdict='Exact minors verified; the fixed prepass does not justify a complete candidate'
    assert blocked==6 and not r['complete_candidate_justified'] and not r['complete_candidate_implemented']
    census=[]
    for s in r['summaries']:
        stats=s['statistics']['symbolic_minor4'];aff=s['statistics']['affine_filter']
        census.append({'cell':s['cell'],'n':s['n'],'seed':s['seed'],'family':s['family'],
            'outside_assignments':1<<(s['n']-4),'minor_survivors':stats['accepted_prefixes'],'affine_survivors':aff['accepted_prefixes'],
            'original_solutions':s['original_solutions'],'degrees':stats['degrees'],'supports':stats['supports'],
            'product_toggles':s['symbolic_work']['product_toggles'],'cancelled_toggles':s['symbolic_work']['cancelled_toggles']})
    summary={k:r[k] for k in ['cells','observations','all_complete','gates','n16_ms','campaign_seconds','peak_worker_rss_bytes',
        'complete_candidate_justified','complete_candidate_implemented','full_ic_cost','production_solver_cost','rho_ratio','calibrated_operation_ratio']}
    summary.update(classification='engineering',verdict=verdict,census=census)
    dump(HERE/'SUMMARY.json',summary)
    dispatcher=next(row for row in r['n16_ms'] if row['variant']=='word_dispatch')
    costs=[[row['variant'],'complete solve' if row['kind']=='solver' else 'necessary filter']+
           [f"{row[f]:.6f}" for f in ['planted','cross_planted','unplanted']]+[f"{dispatcher['planted']/row['planted']:.4f}",'PASS'] for row in r['n16_ms']]
    gates=[[str(g['n']),g['family'],f"{g['median']:.6f}",f"[{g['ci95'][0]:.6f}, {g['ci95'][1]:.6f}]",g['decision']] for g in r['gates']]
    structure=[[str(c['n']),str(c['seed']),c['family'],str(c['outside_assignments']),str(c['minor_survivors']),str(c['affine_survivors']),str(c['original_solutions']),str(c['degrees']),str(c['supports'])] for c in census]
    md=lambda rows:'\n'.join('| '+' | '.join(row)+' |' for row in rows)
    note=f'''# {verdict}

The four selected augmented minors are constructed exactly, including Boolean
product cancellations and degree drops. Every minor truth bit matches independent
numeric elimination. The necessary affine masks also match their separate oracle,
and every original solution survives both filters. All **5,280 observations across
12 discovery systems and 55 arms** complete and verify.

This representation does not justify a complete candidate. In all six measured
size/family groups, even granting all subsequent recovery free leaves the optimistic
reference/preprocessing upper confidence bound below **0.082**. The complete
reference already finishes sooner than this measured prepass. No new complete
solver or holdout campaign was launched, and no speedup is claimed.

## Conditional cost boundary

For a full solver retaining this fixed complete-mask prepass,
`T_complete >= T_prepass`. The table divides the fastest retained complete-reference
total by the symbolic arm's preprocessing time, excluding external mask-validation
time from the denominator to favor the candidate. Intervals use the fixed paired
discovery observations. Instrumentation and cap checks remain part of the measured
compiler; the boundary is conditional on keeping this implementation and architecture.
It does not establish a bound for another compiler, representation or early exit.

| Variables | Family | Median reference / prepass | 95% interval | Decision |
|---:|---|---:|---|---|
{md(gates)}

## Complete solver and necessary-filter costs

All values are milliseconds per cold arm plus external validation, pooled over
the two n16 discovery seeds and eight random/reverse repetitions. The kind column
separates complete solves from filter-only work. The final ratio is descriptive:
the retained dispatcher's planted median divided by the arm's planted median.
It is not a complete-solver speedup for a filter row. All 52 retained solvers remain
in the same binary alongside the three new filter arms.

| Arm | Kind | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Dispatcher / arm, planted | Correctness |
|---|---|---:|---:|---:|---:|---|
{md(costs)}

## Exact structural evidence

The fixed minors leave **1,410–1,920** of 4,096 outside assignments at n16, while
exact affine consistency leaves **5–15**. Both predicates remain necessary only;
the final column gives independently enumerated original solution counts. Four
minors drop to degree five in the cross-planted cases; the remaining 44 have
degree six. No chosen minor is zero or duplicated on this discovery grid. Separate
tests retain those degeneracies and the rank-deficient insufficiency counterexample.

| Variables | Seed | Family | Outside assignments | Minor survivors | Affine survivors | Original solutions | Degrees | ANF supports |
|---:|---:|---|---:|---:|---:|---:|---|---|
{md(structure)}

Individual rejection masks, their exact pairwise intersections and their joint
union are retained in the result. Their rates are not added or assumed independent.
The verifier checks **417,792 original assignments** across the complete discovery
grid. A vanishing determinant cannot substitute for an original solution.

The campaign took **{r['campaign_seconds']:.3f} seconds**. Whole-worker peak RSS was
**{r['peak_worker_rss_bytes']:,} bytes**, covering all arms and oracle preparation.
The symbolic DP coefficient payload peaks at 1,024 bytes at n12 and 16,384 bytes
at n16; these are that component's word arrays, not total candidate memory or
allocator peaks. Other candidate-specific memory remains unmeasured.

Validation comprises 76 Rust tests and 15 Python evidence checks, including
independent truth-table replay, explicit permutation parity, cancellation and
degree-drop cases, false models, changed sources/receipts, overlap accounting and
censored-construction rejection. These are producer checks, not external review.

`README.md` states the exact scope and timing exclusions. The frozen protocol,
source, raw samples, receipts and result are sealed in `run_01/manifest.json`;
`RUN_LEDGER.json` binds this report. The earlier projected-fiber result remains
unchanged. `NEXT_EXPERIMENT.md` records a distinct, unimplemented 16-bit syndrome
representation hypothesis, with mandatory original-equation checks and fresh
holdouts for any promising complete candidate.

Full-IC, production, calibrated-operation and rho costs remain **null**. This
construction diagnostic does not meet the active dramatic-gain objective.
'''
    (HERE/'CONCLUSION.md').write_text(note)
    base='https://github.com/aburan28/crypto/blob/main/research/boolean_implicit_minors_20260923/'
    html_rows=lambda rows:'\n'.join('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in row)+'</tr>' for row in rows)
    section=f'''<section class="panel" id="boolean-implicit-minors-20260923">
  <div class="panel-head">
    <h2>{html.escape(verdict)} <span class="chip">engineering</span></h2>
    <p>All 5,280 observations on 12 discovery systems and 55 arms complete and verify.
      All six measured preprocessing opportunity bounds miss two, even granting later recovery free.
      No complete candidate or holdout promotion follows. Full-IC and rho costs remain null.</p>
    <p>Sources: <a href="{base}CONCLUSION.md">costs, overlaps and limitations</a>,
      <a href="{base}run_01/results.json">frozen result</a>,
      <a href="{base}SUMMARY.json">recorded figures</a>,
      <a href="{base}RUN_LEDGER.json">evidence hashes</a>.</p>
  </div>
  <div class="table-wrap"><table><thead><tr><th>Variables</th><th>Family</th><th>Reference / prepass</th><th>95% interval</th><th>Conditional decision</th></tr></thead><tbody>
{html_rows(gates)}
  </tbody></table></div>
  <p>n16 discovery medians, milliseconds including validation. Kind distinguishes complete solves from necessary filters.
    The dispatcher ratio uses pooled planted medians and is descriptive; a filter row does not establish a full-solver speedup.</p>
  <div class="table-wrap"><table><thead><tr><th>Arm</th><th>Kind</th><th>Planted ms</th><th>Cross-planted ms</th><th>Unplanted ms</th><th>Dispatcher / arm, planted</th><th>Correctness</th></tr></thead><tbody>
{html_rows(costs)}
  </tbody></table></div>
  <p>The four minors leave 1,410&ndash;1,920 outside assignments at n16; exact affine consistency leaves 5&ndash;15.
    Every filter bit is checked and all 417,792 original assignments are covered. Product cancellations, degree drops,
    zero/duplicate controls and rank-deficient insufficiency remain explicit. Whole-worker RSS: {r['peak_worker_rss_bytes']:,} bytes.
    Campaign: {r['campaign_seconds']:.3f} seconds. The prior complete projected-fiber comparison remains below.</p>
</section>

'''
    path=REPO/'docs/index-calculus-scoreboard.html';page=path.read_text();marker='<section class="panel" id="boolean-implicit-minors-20260923">'
    if marker in page:
        begin=page.index(marker);end=page.index('</section>',begin)+len('</section>');page=page[:begin]+page[end:].lstrip('\n')
    anchor='<section class="panel" id="boolean-projected-fibers-20260923">';assert anchor in page
    path.write_text(page.replace(anchor,section+anchor,1))
    dump(HERE/'RUN_LEDGER.json',{'schema_version':1,'classification':'engineering','scope':meta['scope'],
        'manifest_sha256':sha(run/'manifest.json'),'results_sha256':sha(run/'results.json'),
        'source_hashes':meta['source_hashes'],'execution_git_base':meta['git_base'],
        'report_sha256':sha(Path(__file__)),'summary_sha256':sha(HERE/'SUMMARY.json'),'conclusion_sha256':sha(HERE/'CONCLUSION.md'),
        'complete_candidate_implemented':False,'full_ic_cost':None,'rho_ratio':None,'calibrated_operation_ratio':None})
    print(json.dumps({k:summary[k] for k in ['verdict','cells','observations','gates','complete_candidate_justified']}))

if __name__=='__main__':main()
