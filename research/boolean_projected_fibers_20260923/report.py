#!/usr/bin/env python3
"""Render a static note and scoreboard from a verified frozen result."""
import hashlib
import html
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parent.parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dump(path,value):path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def main():
    run=HERE/'run_01'
    manifest=json.loads((run/'manifest.json').read_text())
    assert set(manifest['files'])=={p.name for p in run.iterdir() if p.is_file() and p.name!='manifest.json'}
    for name,digest in manifest['files'].items():assert sha(run/name)==digest,name
    meta=json.loads((run/'metadata.json').read_text())
    for name,digest in meta['source_hashes'].items():assert sha(HERE/name)==digest,name
    r=json.loads((run/'results.json').read_text())
    assert r['cells']==216 and r['samples']==89856
    assert r['all_complete']
    rows=r['n24_holdout_ms']
    dispatcher=next(x for x in rows if x['variant']=='word_dispatch')
    table=[]
    for row in rows:
        table.append([row['variant']]+[f"{row[f]:.6f}" for f in ['planted','cross_planted','unplanted']]+[f"{dispatcher['planted']/row['planted']:.3f}",'PASS'])
    decisions=[[c,str(d['dramatic_groups'])+'/18',str(d['incremental_groups'])+'/18',
                'REJECTED' if not d['dramatic_pass'] else 'CONFIRMATION REQUIRED'] for c,d in r['decisions'].items()]
    gate_rows=[[g['candidate'],g['split'],str(g['n']),g['family'],f"{g['median']:.4f}",f"[{g['ci95'][0]:.4f}, {g['ci95'][1]:.4f}]",'PASS' if g['dramatic_pass'] else 'REJECTED'] for g in r['gates']]
    no_promotion=not any(d['dramatic_pass'] for d in r['decisions'].values())
    verdict='Complete recovery is verified; the dramatic-gain gate is rejected' if no_promotion else 'A finite candidate requires fresh unchanged-source confirmation'
    work=[]
    for c in r['decisions']:
        entries=[s['projected_work'][c] for s in r['summaries'] if s['n']==24 and s['split']=='holdout']
        work.append({'variant':c,'cells':len(entries),'prefixes':sum(w['prefixes'] for w in entries),
            'screen_rejected':sum(w['screen_rejected'] for w in entries),'affine_queries':sum(w['affine_queries'] for w in entries),
            'affine_rejected':sum(w['affine_rejected'] for w in entries),'extensions_checked':sum(w['extensions_checked'] for w in entries),
            'original_rejected':sum(w['original_rejected'] for w in entries),
            'annihilated_ranks':sorted({w['annihilated_rank'] for w in entries}),
            'quotient_dimensions':sorted({w['quotient_dimension'] for w in entries})})
    summary={k:r[k] for k in ['all_complete','cells','samples','decisions','campaign_seconds','peak_worker_rss_bytes','production_solver_cost','full_ic_cost','rho_ratio','calibrated_operation_ratio']}
    summary.update(classification='engineering',verdict=verdict,n24_holdout_ms=rows,n24_holdout_work_once_per_fixture=work)
    dump(HERE/'SUMMARY.json',summary)
    md=lambda values:'\n'.join('| '+' | '.join(row)+' |' for row in values)
    note=f'''# {verdict}

This implements the exact projected-fiber contract on bounded generated Boolean
quadratic systems. All **{r['samples']:,} observations on {r['cells']} systems and
52 methods** completed with verified outcomes. Each candidate computes a necessary
affine system, recovers its complete free-variable space, and checks the original
equations before accepting a solution. Degree changes and coefficient cancellations
are exact XOR operations; a rank drop never discards required extensions.

| Candidate | Groups above 2.0 | Groups above 1.0 | Dramatic decision |
|---|---:|---:|---|
{md(decisions)}

The gate compares against the pointwise fastest of all 49 retained references,
with a paired 95% lower confidence bound in every one of 18 regression/holdout
groups. No source was tuned after holdout timing. These are finite generic-solver
measurements, not calibrated operation ratios, asymptotic results or a full
index-calculus result. A positive group does not establish a universal gain.

## Complete cold comparison on n24 fresh holdouts

Times are milliseconds per complete cold solve plus validation, pooled across two
fresh seeds and eight random/reverse repetitions in each family. The last ratio is
the retained dispatcher's pooled planted median divided by the displayed arm's
median; it is descriptive and is not the acceptance statistic. All 52 methods
remain visible. The original fixtures and all 192 predecessor cases are retained.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Dispatcher / arm, planted | Correctness |
|---|---:|---:|---:|---:|---|
{md(table)}

## Mathematical work and limits

Projection removes the quadratic terms internal to the chosen low block, making
the projected subsystem affine. It can also erase constraints. The solver therefore
returns a particular affine answer **and a complete kernel basis**, enumerates the
remaining freedom, and verifies all original equations. The `xy+1` counterexample
is a regression test: a zero projected system must not be accepted at x=y=0.
Caps remain UNKNOWN. The measured suite completed without censored observations.

`SUMMARY.json` records n24 holdout screening, affine elimination and original-check
counts once per fixture, without multiplying by benchmark repetitions. Its ranks
and quotient dimensions are exact structural quantities; they are not predictions
of success or counts of independent conditional constraints.

The campaign took **{r['campaign_seconds']:.3f} seconds**. Whole-worker peak RSS was
**{r['peak_worker_rss_bytes']:,} bytes**, including all methods and common references.
Candidate-specific allocator peaks and calibrated operation costs are unmeasured.
Temporary solver-workspace destruction occurs inside solve time. Fixture generation,
reference preparation, serialization and returned-diagnostic destruction are outside
arm clocks and inside worker process receipts. No phase cost is subtracted.

## Paired gates

| Candidate | Split | Variables | Family | Median reference / candidate | 95% interval | Above 2.0 |
|---|---|---:|---|---:|---|---|
{md(gate_rows)}

The 69 Rust tests include exhaustive complete affine fibers, every pair of
three-variable quadratics, changing coefficients and rank, direct projection
identities, original-equation recovery and caps, together with all retained tests.
The 15 Python evidence checks replay the result, independently regenerate the fixtures,
verify source custody and reject altered models, work, receipts and capped claims.
These are producer checks, not an external review.

The earlier support-envelope construction and construction-cost diagnostic remain
unchanged. This result implements the previously prospective scan; it does not
convert their negative or inconclusive timing cells into gains. See `README.md`
for the exact domain, fixed selection, accounting and replay contract.

`NEXT_EXPERIMENT.md` records a separate unimplemented hypothesis: capped implicit
determinant filters for affine consistency. It requires exact changing-coefficient
products, Boolean cancellations, rank-deficient controls and charged construction
before any complete-solver comparison. No performance conclusion follows from it.

Production solver, full index-calculus, calibrated-operation and rho costs remain
**null**. The original dramatic-gain goal remains open.
'''
    (HERE/'CONCLUSION.md').write_text(note)
    base='https://github.com/aburan28/crypto/blob/main/research/boolean_projected_fibers_20260923/'
    html_rows=lambda values:'\n'.join('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in row)+'</tr>' for row in values)
    section=f'''<section class="panel" id="boolean-projected-fibers-20260923">
  <div class="panel-head">
    <h2>{html.escape(verdict)} <span class="chip">engineering</span></h2>
    <p>Exact equation-space projection, affine kernel recovery and original-equation checks are implemented.
      All {r['samples']:,} observations over 216 generated systems and 52 methods complete and verify.
      The table reports complete generic solver time; full-IC and rho costs remain null.</p>
    <p>Sources: <a href="{base}CONCLUSION.md">results and paired gates</a>,
      <a href="{base}run_01/results.json">frozen measurements</a>,
      <a href="{base}SUMMARY.json">recorded figures</a>,
      <a href="{base}RUN_LEDGER.json">source and evidence hashes</a>.</p>
  </div>
  <div class="table-wrap"><table><thead><tr><th>Candidate</th><th>Groups above 2.0</th><th>Groups above 1.0</th><th>Dramatic decision</th></tr></thead><tbody>
{html_rows(decisions)}
  </tbody></table></div>
  <p>Fresh n24 holdouts; milliseconds per cold solve plus validation. The dispatcher ratio uses pooled planted medians
    and is descriptive. Acceptance requires all 18 paired groups, shown in the linked note. No tuning followed holdout timing.</p>
  <div class="table-wrap"><table><thead><tr><th>Method</th><th>Planted ms</th><th>Cross-planted ms</th><th>Unplanted ms</th><th>Dispatcher / arm, planted</th><th>Correctness</th></tr></thead><tbody>
{html_rows(table)}
  </tbody></table></div>
  <p>Projection is necessary only: every free extension remains subject to the original equations.
    Quotient dimensions do not establish independent conditional constraints or predicted yield.
    Whole-worker peak RSS: {r['peak_worker_rss_bytes']:,} bytes; campaign: {r['campaign_seconds']:.3f} seconds.
    These are finite mathematical solver diagnostics. Earlier negative and inconclusive construction results remain below.</p>
</section>

'''
    page=REPO/'docs/index-calculus-scoreboard.html'
    text=page.read_text()
    marker='<section class="panel" id="boolean-projected-fibers-20260923">'
    if marker in text:
        start=text.index(marker);end=text.index('</section>',start)+len('</section>')
        text=text[:start]+text[end:].lstrip('\n')
    old='<section class="panel" id="boolean-schedule-construction-20260923">'
    assert old in text
    page.write_text(text.replace(old,section+old,1))
    dump(HERE/'RUN_LEDGER.json',{'schema_version':1,'classification':'engineering','scope':meta['scope'],
        'run_manifest_sha256':sha(run/'manifest.json'),'results_sha256':sha(run/'results.json'),
        'source_hashes':meta['source_hashes'],'git_execution_base':meta['git_base'],
        'report_sha256':sha(Path(__file__)),'summary_sha256':sha(HERE/'SUMMARY.json'),'conclusion_sha256':sha(HERE/'CONCLUSION.md'),
        'validation':'Producer correctness and evidence replay; no external review implied.',
        'full_ic_cost':None,'rho_ratio':None,'calibrated_operation_ratio':None})
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
