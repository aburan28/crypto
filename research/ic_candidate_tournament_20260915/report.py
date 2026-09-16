#!/usr/bin/env python3
"""Audit a frozen campaign and write its static research table and scoreboard panel."""
import argparse
import hashlib
import html
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys


def read(path):
    return json.loads(path.read_text())


def geometric(values):
    return math.exp(statistics.mean(math.log(x) for x in values))


def frozen_write(path, value):
    data=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if path.exists():
        if path.read_text()!=data:
            raise ValueError('refusing to replace frozen reporting data')
    else:
        with path.open('x') as f:f.write(data)


def stage_table(root,stage):
    records=[read(p) for p in sorted((root/'runs'/stage).glob('*/*/rep-*/receipt.json'))]
    summary=read(root/'summaries'/f'{stage}.json')
    comparisons={c['candidate']:c for c in summary['comparisons']}
    rho=summary.get('rho_over_incumbent')
    ratios={'incumbent':1}
    ratios.update({name:c.get('candidate_over_baseline') for name,c in comparisons.items()})
    ratios['rho']=rho.get('candidate_over_baseline') if rho else None
    rows=[]
    names=['incumbent']+sorted({r['arm'] for r in records}-{'incumbent','rho'})+['rho']
    for name in names:
        rs=[r for r in records if r['arm']==name]
        if not rs:continue
        complete=all(r['status']=='VERIFIED' and r['total_operations'] is not None for r in rs)
        normalized={};floor={}
        if complete:
            groups={}
            for r in rs:groups.setdefault((r['cell'],r['case']),[]).append(r)
            for (cell,_),samples in groups.items():
                normalized.setdefault(cell,[]).append(statistics.median(r['normalized_S'] for r in samples))
                if all(r.get('ratio_to_floor') is not None for r in samples):
                    floor.setdefault(cell,[]).append(statistics.median(r['ratio_to_floor'] for r in samples))
        ratio=ratios.get(name)
        rho_ratio=ratios.get('rho')
        rows.append({'variant':name,'stage':stage,'S_Ir':geometric([geometric(v) for v in normalized.values()]) if normalized else None,
                     'candidate_over_incumbent':ratio,
                     'candidate_over_rho':ratio/rho_ratio if ratio and rho_ratio else None,
                     'candidate_over_floor':geometric([geometric(v) for v in floor.values()]) if floor else None,
                     'verified_runs':sum(r['status']=='VERIFIED' for r in rs),'scheduled_runs':len(rs),
                     'class':'reference' if name in ('incumbent','rho') else 'engineering experiment'})
    return rows


def native_table(root, stage):
    records=[read(p) for p in sorted((root/'runs'/stage).glob('*/*/rep-*/receipt.json'))]
    summary=read(root/'summaries'/f'{stage}.json')
    comps={x['candidate']:x for x in summary['comparisons']}
    comps['rho']=summary.get('rho_over_incumbent') or {}
    rho_ratio=comps['rho'].get('native_wall_candidate_over_baseline')
    rows=[]
    for name in ['incumbent']+sorted({r['arm'] for r in records}-{'incumbent','rho'})+['rho']:
        rs=[r for r in records if r['arm']==name]
        if not rs: continue
        complete=all(r['status']=='VERIFIED' for r in rs)
        groups={}; cells={}
        if complete:
            for r in rs: groups.setdefault((r['cell'],r['case']),[]).append(r['native_process']['process_wall_seconds'])
            for (cell,_),values in groups.items(): cells.setdefault(cell,[]).append(statistics.median(values))
        comparison=comps.get(name,{})
        ratio=1 if name=='incumbent' and complete else comparison.get('native_wall_candidate_over_baseline')
        rows.append({'variant':name,'native_seconds':geometric([geometric(v) for v in cells.values()]) if cells else None,
                     'over_incumbent':ratio,'ci95_over_incumbent':comparison.get('native_wall_ci95'),
                     'over_rho':ratio/rho_ratio if ratio and rho_ratio else None,
                     'verified':sum(r['status']=='VERIFIED' for r in rs),'scheduled':len(rs)})
    return rows


def native_markdown(rows):
    lines=['| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |',
           '|---|---:|---:|---|---:|---:|']
    for r in rows:
        interval=r['ci95_over_incumbent']
        ci='['+', '.join(fmt(x) for x in interval)+']' if interval else 'reference'
        lines.append('| '+ ' | '.join([r['variant'],fmt(r['native_seconds']*1000 if r['native_seconds'] else None),
                     fmt(r['over_incumbent']),ci,fmt(r['over_rho']),f"{r['verified']}/{r['scheduled']}"])+ ' |')
    return '\n'.join(lines)


def matched_base_audit(root):
    count=0
    for stage in ('aa','smoke','development','selection','confirmation','replay'):
        for case in (root/'runs'/stage).iterdir():
            hashes=set()
            for path in case.glob('*/rep-*/profile/stdout.json'):
                report=read(path)
                if report.get('mode')!='ic' or report.get('status')!='complete':
                    continue
                points=sorted(tuple(map(int,p)) for p in report['factor_base'])
                hashes.add(hashlib.sha256(json.dumps(points,separators=(',',':')).encode()).hexdigest())
                count+=1
            if len(hashes)>1:
                raise ValueError('changed paired factor-base support in '+str(case))
    return {'status':'VERIFIED','complete_ic_profiles_checked':count}


def rho_health(root):
    cells={}
    for path in sorted((root/'runs/confirmation').glob('*/rho/rep-*/profile/stdout.json')):
        report=read(path)
        fixture=report['fixture']
        degree=fixture['degree']
        expected=math.sqrt(math.pi*int(fixture['subgroup_order'])/(2*(2*degree)))
        key=f"n{degree}a{fixture['curve_a']}"
        for row in report.get('solutions',[]):
            cells.setdefault(key,[]).append((row['iterations']/expected,
                row['walk_group_additions']/expected,row['restarts']))
    return {key:{'expected_steps_model':'sqrt(pi*r/(2*A)), A=2*n',
                 'median_iterations_over_model':statistics.median(x[0] for x in values),
                 'median_walk_additions_over_model':statistics.median(x[1] for x in values),
                 'max_restarts':max(x[2] for x in values),'verified_runs':len(values)}
            for key,values in cells.items()}


def fmt(value):
    return 'unmeasured' if value is None else f'{value:,.4g}'


def markdown_table(rows):
    out=['| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |',
         '|---|---:|---:|---:|---:|---:|---|']
    for r in rows:
        out.append('| '+ ' | '.join([r['variant'],fmt(r['S_Ir']),fmt(r['candidate_over_incumbent']),
            fmt(r['candidate_over_rho']),fmt(r['candidate_over_floor']),f"{r['verified_runs']}/{r['scheduled_runs']}",r['class']])+' |')
    return '\n'.join(out)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--round',type=Path,required=True)
    p.add_argument('--scoreboard',type=Path)
    args=p.parse_args()
    root=args.round.resolve()
    result=read(root/'decision.json')
    contract=read(root/'contract.json')
    audit_path=root/'audit.json'
    # Always audit before publication, including raw profiles and certificates.
    # A prior receipt hash alone cannot detect a subsequently changed raw file.
    evidence={str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest()
              for x in sorted(root.glob('runs/**/receipt.json'))}
    evidence['decision.json']=hashlib.sha256((root/'decision.json').read_bytes()).hexdigest()
    audit_seal=root/'audit-inputs.json'
    completed=subprocess.run([sys.executable,str(root/'evaluator/tournament.py'),'verify','--round',str(root)],
                             capture_output=True,text=True,check=True)
    audit=json.loads(completed.stdout)
    if audit.get('status')!='VERIFIED':raise ValueError('audit failed')
    audit_path.write_text(json.dumps(audit,indent=2)+'\n')
    audit_seal.write_text(json.dumps(evidence,indent=2,sort_keys=True)+'\n')
    rows={stage:stage_table(root,stage) for stage in ('development','confirmation','replay')}
    data={'schema_version':1,'round':root.name,'unit':contract['unit'],'target_count':contract.get('target_count',1),'decision':result,
          'audit':audit,'matched_base_audit':matched_base_audit(root),'tables':rows,'rho_health':rho_health(root),
          'limits':['Fixed compiler/ISA user-space instructions; kernel/device/external audit work excluded.',
                    'Bounded Koblitz configuration result; no arithmetic exponent or family-wide claim.',
                    'Native wall times are diagnostics; profiled elapsed time is not native performance.']}
    if contract.get('native_timing_protocol'):
        data['native_tables']={stage:native_table(root,stage) for stage in rows}
        data['limits'][-1]='Native cold-process time is secondary, with paired curve/target confidence intervals; profiler time is excluded.'
    frozen_write(root/'measurements.json',data)
    winner=result['winner'] or 'none (inconclusive)'
    text=[f'# IC candidate tournament: {root.name}', '',f"Decision: **{result['status']} — {winner}**.",'',
          'This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).',
          'Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.','']
    text += [f"Each complete cold job recovers **{contract.get('target_count',1)} target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.", '']
    if result.get('instruction_speedup'):
        text += [f"The confirmed instruction speedup is **{result['instruction_speedup']:.3f}x**; "
                 f"candidate/incumbent ratio {result['confirmation']['candidate_over_baseline']:.4f}, "
                 f"paired 95% interval {result['confirmation']['ci95']}. This is an engineering result on the tested workloads.",'']
    text += [f"Independent audit checked **{audit['trial_receipts']} trial receipts**. "
             'Repetitions are grouped within fixtures; intervals resample curve cells and their targets.', '',
             'The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.', '',
             'The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.', '',
             'The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.', '',
             f"Observed rho/winner instruction ratio: **{fmt(result.get('rho_over_winner'))}**. "
             'A value below one means rho costs less. This is not an extrapolated crossover.', '']
    for stage,rs in rows.items():
        text += [f'## {stage.title()}','',markdown_table(rs),'']
    if 'native_tables' in data:
        text += ['## Native runtime and rho parity','',
                 'Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.', '',
                 'Parity verdict: **'+str(result.get('rho_parity',False))+'**. '+result.get('parity_definition',''), '']
        for stage,rs in data['native_tables'].items():
            text += [f'### {stage.title()} native time','',native_markdown(rs),'']
        parity=result.get('winner_over_rho',{})
        for stage,pair in parity.items():
            if pair.get('eligible'):
                text += [f"{stage.title()} winner/rho: instructions {pair['candidate_over_baseline']:.4f}, CI {pair['ci95']}; native time {pair['native_wall_candidate_over_baseline']:.4f}, CI {pair['native_wall_ci95']}.", '']
    text += ['## Rho health','',
             'All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.','',
             '| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |',
             '|---|---:|---:|---:|---:|']
    for cell,health in data['rho_health'].items():
        text.append(f"| {cell} | {health['median_iterations_over_model']:.3f} | {health['median_walk_additions_over_model']:.3f} | {health['max_restarts']} | {health['verified_runs']} |")
    text += ['','## Evidence','',
             '- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).',
             '- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).',
             '- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).',
             '- [Operating commands and skills](' + Path(__import__('os').path.relpath(Path(__file__).resolve().parent/'OPERATIONS.md',root)).as_posix() + ').','']
    if (root.parent/'round-0001').is_dir():
        text += ['The earlier snapshot/build attempt is retained under `../round-0001`.']
    (root/'REPORT.md').write_text('\n'.join(text)+'\n')
    if args.scoreboard:
        board=args.scoreboard.resolve()
        relative=Path(__import__('os').path.relpath(root,board.parent)).as_posix()
        marker='ic-tournament-'+root.name
        parts=[f'<!-- BEGIN {marker} -->',f'<section class="panel" id="{marker}">',
               '<div class="panel-head"><h2>Complete IC candidate tournament</h2>',
               f'<p><strong>{html.escape(result["status"])}: {html.escape(winner)}</strong>. '
               'Cold complete-DLP instruction costs; fixed compiler and amd64 ISA. '
               f"Targets per cold job: {contract.get('target_count',1)}. No arithmetic-complexity or family-wide claim.</p>",
               f'<p>Evidence: <a href="{relative}/REPORT.md">full result</a> · '
               f'<a href="{relative}/measurements.json">frozen reporting data</a> · '
               f'<a href="{relative}/audit.json">independent audit</a>.</p></div>']
        for stage in ('development','confirmation'):
            parts += ['<div class="table-scroll"><table>',f'<caption>{stage.title()} — instruction unit Ir; all variants on matched fixtures.</caption>',
                      '<thead><tr><th>Variant</th><th>S (Ir / √r)</th><th>Cost / incumbent</th><th>Cost / rho</th><th>Cost / floor</th><th>Verified</th><th>Class</th></tr></thead><tbody>']
            for r in rows[stage]:
                values=[r['variant'],fmt(r['S_Ir']),fmt(r['candidate_over_incumbent']),fmt(r['candidate_over_rho']),
                        fmt(r['candidate_over_floor']),f"{r['verified_runs']}/{r['scheduled_runs']}",r['class']]
                parts.append('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in values)+'</tr>')
            parts += ['</tbody></table></div>']
        if 'native_tables' in data:
            parts += ['<h3>Native cold process time — confirmation</h3>',
                      '<p>Paired confidence intervals; fixed CPU/resources. Rho parity: '+str(result.get('rho_parity',False))+'.</p>',
                      '<div class="table-scroll"><table><thead><tr><th>Variant</th><th>ms</th><th>Time / incumbent</th><th>95% CI</th><th>Time / rho</th><th>Verified</th></tr></thead><tbody>']
            for r in data['native_tables']['confirmation']:
                interval=r['ci95_over_incumbent']
                values=[r['variant'],fmt(r['native_seconds']*1000 if r['native_seconds'] else None),fmt(r['over_incumbent']),
                        str(interval) if interval else 'reference',fmt(r['over_rho']),f"{r['verified']}/{r['scheduled']}"]
                parts.append('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in values)+'</tr>')
            parts += ['</tbody></table></div>']
        parts += ['<p>Kernel/device work and external audit are outside Ir. Native timings follow their separate protocol. '

                  'The K-instruction floor applies only to this full-rank collector and cannot establish an algorithmic advance.</p>',
                  '</section>',f'<!-- END {marker} -->']
        chunk='\n'.join(parts)+'\n'
        page=board.read_text()
        if f'<!-- BEGIN {marker} -->' in page:
            start=page.index(f'<!-- BEGIN {marker} -->');end=page.index(f'<!-- END {marker} -->')+len(f'<!-- END {marker} -->')
            page=page[:start]+chunk.rstrip()+page[end:]
        else:
            if '</body>' not in page:raise ValueError('scoreboard has no body closing tag')
            page=page.replace('</body>',chunk+'</body>',1)
        board.write_text(page)
    print(root/'REPORT.md')


if __name__=='__main__':main()
