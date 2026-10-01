#!/usr/bin/env python3
"""Render only hash-bound historical or qualified half-word evidence."""
import hashlib
import html
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parent.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,v):p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def bundle(path):
    manifest=read(path/'manifest.json')
    for name,digest in manifest['files'].items():assert sha(path/name)==digest,name
    return read(path/'results.json')
def main():
    first=bundle(HERE/'probe_01');second=bundle(HERE/'probe_02')
    registry=read(HERE/'QUALIFIED_RUNS.json');qualified={};qualified_history=[]
    for entry in registry.get('historical_qualified_runs',[]):
        path=HERE/entry['path'];assert sha(path/'manifest.json')==entry['manifest_sha256']
        result=bundle(path);assert result['qualified'] and result['all_complete']
        qualified_history.append((entry,result,read(path/'metadata.json')))
    failures=[]
    for entry in registry['failed_attempts']:
        path=HERE/entry['path'];assert sha(path/'manifest.json')==entry['manifest_sha256']
        manifest=read(path/'manifest.json')
        for name,digest in manifest['files'].items():assert sha(path/name)==digest,name
        status=read(path/'EXECUTION_STATUS.json')
        failures.append({'path':entry['path'],'github_run_id':entry['github_run_id'],'github_attempt':entry['github_attempt'],
                         'status':status,'reason':entry['reason']})
    for phase in ['discovery','full']:
        entry=registry[phase]
        if entry is None:continue
        path=HERE/entry['path'];assert sha(path/'manifest.json')==entry['manifest_sha256']
        result=bundle(path);assert result['qualified'] and result['phase']==phase
        qualified[phase]=result
    latest=qualified.get('full',qualified.get('discovery'))
    title=('16-bit predecessor reaches about 1.55x on n24; specialization pending' if qualified_history else '16-bit syndrome implementation verified; qualified comparison pending') if latest is None else '16-bit syndrome qualified comparison recorded; full objective remains open'
    old1={r['variant']:r for r in first['n24_ms']};old2={r['variant']:r for r in second['n24_ms']}
    protocol=read(HERE/'protocol_full.json')
    historical=[]
    for arm in protocol['variants']:
        a,b=old1.get(arm),old2.get(arm)
        historical.append([arm]+([f"{a['planted']:.6f}",f"{b['planted']:.6f}",f"{b['cross_planted']:.6f}",f"{b['unplanted']:.6f}",'historical; unqualified'] if a and b else ['pending']*4+['not measured']))
    gate_rows=[[a,str(d['dramatic_groups'])+'/9',str(d['incremental_groups'])+'/9'] for a,d in second['decisions'].items()]
    md=lambda rows:'\n'.join('| '+' | '.join(row)+' |' for row in rows)
    report=f'''# {title}

The complete 16-bit syndrome implementation preserves the original model and
assignment order and checks every projected zero on the original equations. The
current source additionally provides runtime-gated three-input-XOR kernels for
both 32-bit and 16-bit lanes, with portable fallbacks. No global CPU feature flag
changes the retained control implementations.

The current code passes 79 Rust correctness checks. Resource/evidence tests cover
the exact historical replay, source identity, full-word agreement, corrupted work,
missing/contended isolation records, false resource flags, and A/A work drift.
These are producer checks, not an independent external review or a performance gain.

## Historical discovery, preserved without promotion

Each of the two September discovery probes completed 24 systems and 10,944
observations across 57 methods. Probe 1's 64-point native variant passed six
incremental groups; probe 2 passed five. Neither produced a dramatic group. Probe
2's outlined mutable-result helper regressed. Its source-based explanation and
the subsequent read-only-helper change are in `SOURCE_EVOLUTION.md`.

The repository's current isolation policy requires pinned/reserved CPUs and A/A
calibration. Those older probes lack that evidence. Their raw timings and original
decisions remain immutable historical diagnostics; they do not discharge the
current performance gate and are not pooled with new runs. Their contention is
unknown. `ACCOUNTING_NOTE.md` corrects one prose exclusion: repeat-signature checks
were actually charged inside validation and total time. No numbers were rewritten.

| Probe 2 arm | Historical groups above 2.0 | Historical groups above 1.0 |
|---|---:|---:|
{md(gate_rows)}

## Preserved n24 timing table

Units are milliseconds per complete generic solve plus validation, pooling two
discovery seeds and eight repetitions. These are **unisolated historical values**,
not eligible new performance evidence. The first planted column retains the initial
figure beside the second probe. New hardware variants have no historical value.

| Method | Probe 1 planted ms | Probe 2 planted ms | Probe 2 cross-planted ms | Probe 2 unplanted ms | Qualification |
|---|---:|---:|---:|---:|---|
{md(historical)}

## Current measurement contract

`isolated_run.py` refuses unsupported platforms. It freezes the sources and protocol,
keeps the rebuilt unmodified baseline and actual comparison/test binaries, builds
and tests under the repository's busy lock, and then runs A/A before A/B with pinned
CPU reservations. It records host identity/features, memory, compiler, exact commands,
pressure, context switches, faults and other-process CPU use. Any failed, contended
or missing-isolation stage is retained and stops that campaign; it is not pooled
or silently retried. The A/A symmetric noise spread is an additional gate.

The comparison retains the 52 September 23 solver implementations as a frozen
reference roster. It does not claim to include every later repository optimization
or every known solver. Nine new treatments include matched 32-bit EOR3 controls so
the hardware operation is not attributed solely to narrowing. Runtime capability
records distinguish an actual feature path from a portable fallback.

The intended discovery grid has 61 A/B arms and 24 fixtures: 11,712 comparison
observations plus 384 A/A calibration observations. The full grid has 240 fixtures,
117,120 comparison observations and 3,840 calibration observations. It retains all
216 earlier inputs and adds 24 unused holdouts. The full run must bind unchanged
timed source from qualified discovery. A positive primary comparison still requires
confirmation on further unused holdouts. No source tuning follows holdout timing.

Native timing currently requires Linux affinity and pressure interfaces, so the
local macOS correctness checks do not supply new qualified timings. The Linux ARM64
CI route is described in `RESOURCE_PLAN.md`. VM neighbours and frequency remain
outside the tool's control; the measured A/A spread and hardware class must accompany
any result.

`QUALIFIED_RUNS.json` binds any accepted new run; null means no such measurement
is recorded. `RUN_LEDGER.json` binds this report and the preserved probes. All
full-IC, production, calibrated-operation and rho costs remain **null**. The active
dramatic-gain objective is not achieved by this implementation or its tests.
'''
    summary={'schema_version':1,'classification':'engineering','status':title,'historical_cells':48,'historical_observations':21888,
        'historical_qualification':'unisolated; not eligible under current policy','historical_probe_02_decisions':second['decisions'],
        'qualified_runs':registry,'qualified_decisions':None if latest is None else latest['decisions'],
        'full_ic_cost':None,'production_solver_cost':None,'calibrated_operation_ratio':None,'rho_ratio':None}
    history_html=[];history_summary=[]
    for entry,result,meta in qualified_history:
        reference=next(row for row in result['n24_ms'] if row['variant']=='word_dispatch')
        rows=[[row['variant']]+[f"{row[f]:.6f}" for f in ['planted','cross_planted','unplanted']]+[f"{reference['planted']/row['planted']:.4f}",'PASS'] for row in result['n24_ms']]
        report+='\n## Qualified predecessor: '+entry['source_head']+'\n\n'
        report+='The preserved Linux ARM64 discovery completed 24 fixtures, 11,712 comparison observations and 384 A/A observations. Every resource record passed. The 64-point native policy is about 1.55x the retained dispatcher on pooled n24 medians, but no dramatic group passed. These results belong to the predecessor source, before the complete-budget specialization; they do not measure the current source or discharge the full/holdout gate.\n\n'
        report+='All values below are n24 discovery milliseconds per cold solve plus validation. The ratio is the pooled planted dispatcher median divided by the arm median; the paired gates and A/A noise floors are retained in the result.\n\n| Method | Planted ms | Cross-planted ms | Unplanted ms | Dispatcher / arm, planted | Correctness |\n|---|---:|---:|---:|---:|---|\n'+md(rows)+'\n'
        host={'source_head':entry['source_head'],'phase':result['phase'],'capabilities':meta['capabilities'],'host':meta['host'],'cpus':meta['cpus'],
              'decisions':result['decisions'],'manifest_sha256':entry['manifest_sha256']}
        history_summary.append(host)
        history_html.append('<p>Qualified predecessor '+html.escape(entry['source_head'][:12])+': Linux ARM64, 24 discovery fixtures, 11,712 comparison observations plus A/A. No dramatic group passed. These numbers precede the current specialization.</p><div class="table-wrap"><table><thead><tr><th>Method</th><th>Planted ms</th><th>Cross-planted ms</th><th>Unplanted ms</th><th>Dispatcher / arm, planted</th><th>Correctness</th></tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>')
    summary['historical_qualified_runs']=history_summary
    if latest is not None:
        report+='\n## Qualified run readback\n\n'+json.dumps({'phase':latest['phase'],'cells':latest['cells'],'observations':latest['observations'],'decisions':latest['decisions']},indent=2)+'\n'
    if failures:
        report+='\n## Retained resource failures\n\n'
        for failure in failures:
            report+=f"- GitHub run {failure['github_run_id']}, attempt {failure['github_attempt']}: {failure['reason']} The complete artifact is retained in `{failure['path']}` and contributes no accepted timing samples.\n"
        report+='\n`ISOLATION_ATTEMPTS.md` records the exact failure and any subsequent complete same-source retry.\n'
    summary['retained_resource_failures']=failures
    (HERE/'CONCLUSION.md').write_text(report);dump(HERE/'SUMMARY.json',summary)
    base='https://github.com/aburan28/crypto/blob/main/research/boolean_halfword_20260923/'
    html_rows=lambda rows:'\n'.join('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in row)+'</tr>' for row in rows)
    section=f'''<section class="panel" id="boolean-halfword-20260923">
  <div class="panel-head"><h2>{html.escape(title)} <span class="chip">engineering</span></h2>
    <p>Exact half-word projection, original-equation checks and feature-gated EOR3 kernels are implemented.
      The two preserved probes contain 21,888 observations, but lack current CPU-isolation and A/A evidence.
      No new timing gain or full-IC result is established. Contention in the old runs is unknown.</p>
    <p>Sources: <a href="{base}CONCLUSION.md">scope and results</a>, <a href="{base}QUALIFIED_RUNS.json">qualified-run registry</a>,
      <a href="{base}RESOURCE_PLAN.md">resource and calibration plan</a>, <a href="{base}RUN_LEDGER.json">evidence hashes</a>.</p></div>
  {''.join(history_html)}
  <p>Historical n24 milliseconds, explicitly unqualified under current policy. Prior values remain visible.
    New variants are pending; null or pending is never a zero-cost result.</p>
  <div class="table-wrap"><table><thead><tr><th>Method</th><th>Probe 1 planted ms</th><th>Probe 2 planted ms</th><th>Probe 2 cross-planted ms</th><th>Probe 2 unplanted ms</th><th>Qualification</th></tr></thead><tbody>
{html_rows(historical)}
  </tbody></table></div>
  <p>The 240-fixture primary comparison and further confirmation remain separate gates. Every projected hit is checked
    against the original system. Feature detection has a portable fallback. Full-IC and rho costs remain unmeasured.</p>
</section>

'''
    path=REPO/'docs/index-calculus-scoreboard.html';page=path.read_text();marker='<section class="panel" id="boolean-halfword-20260923">'
    if marker in page:
        begin=page.index(marker);end=page.index('</section>',begin)+len('</section>');page=page[:begin]+page[end:].lstrip('\n')
    anchor='<section class="panel" id="boolean-implicit-minors-20260923">';assert anchor in page
    path.write_text(page.replace(anchor,section+anchor,1))
    dump(HERE/'RUN_LEDGER.json',{'schema_version':1,'classification':'engineering',
        'probe_01_manifest_sha256':sha(HERE/'probe_01/manifest.json'),'probe_02_manifest_sha256':sha(HERE/'probe_02/manifest.json'),
        'qualified_registry_sha256':sha(HERE/'QUALIFIED_RUNS.json'),'summary_sha256':sha(HERE/'SUMMARY.json'),
        'conclusion_sha256':sha(HERE/'CONCLUSION.md'),'report_sha256':sha(Path(__file__)),
        'current_rust_source_hashes':{p.name:sha(p) for p in sorted(HERE.glob('*.rs'))},
        'goal_complete':False,'full_ic_cost':None,'rho_ratio':None,'calibrated_operation_ratio':None})
    print(title)


if __name__=='__main__':main()
