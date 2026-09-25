"""Render frozen exported rows; no benchmark execution or promotion decisions."""
import argparse
import html
import json
from pathlib import Path


def fmt(value):
    if value is None:
        return 'unknown'
    if isinstance(value, (list,tuple)):
        return '['+', '.join(fmt(v) for v in value)+']'
    return format(value,'.6g') if isinstance(value,float) else str(value)


def comparison(row):
    return row.get('comparison') or {}


def columns(unit):
    if unit=='online':
        return [('Variant',lambda r:r['alias']),('Online ms',lambda r:r['online_ms']),
            ('/ IC reference',lambda r:r['online_over_ic']),
            ('Descriptive 95% interval',lambda r:(comparison(r).get('online') or {}).get('ci95')),
            ('Online rho / variant',lambda r:r['rho_online_over_IC_online']),
            ('Verified / scheduled',lambda r:f"{r['verified']}/{r['scheduled']}")]
    if unit=='instructions':
        return [('Variant',lambda r:r['alias']),('Complete cold Ir',lambda r:r['cold_Ir']),
            ('S = Ir / sqrt(r)',lambda r:r['S_Ir_per_sqrt_r']),
            ('/ IC reference',lambda r:r['cold_Ir_over_ic']),
            ('Descriptive 95% interval',lambda r:comparison(r).get('ci95')),
            ('/ cold rho',lambda r:r['cold_Ir_over_cold_rho']),
            ('/ K floor',lambda r:'not applicable' if r['mode']=='rho' else r['Ir_over_K_floor']),
            ('Verified / scheduled',lambda r:f"{r['verified']}/{r['scheduled']}")]
    return [('Variant',lambda r:r['alias']),('Complete cold ms',lambda r:r['cold_ms']),
            ('/ IC reference',lambda r:r['cold_time_over_ic']),
            ('Descriptive 95% interval',lambda r:comparison(r).get('native_wall_ci95')),
            ('Verified / scheduled',lambda r:f"{r['verified']}/{r['scheduled']}")]


def tables(rows, unit, caption):
    cols=columns(unit) + [('Class',lambda r:'accounting' if r['alias'] in ('incumbent','aa_control','compatibility') or r['mode']=='rho' else 'engineering')]
    heads=[name for name,_ in cols]
    values=[[fmt(fn(row)) for _,fn in cols] for row in rows]
    for row,vals in zip(rows,values):
        if row['alias']=='incumbent':
            if 'Descriptive 95% interval' in heads: vals[heads.index('Descriptive 95% interval')]='reference'
    md=['| '+' | '.join(heads)+' |','|'+'|'.join('---' for _ in heads)+'|']
    md += ['| '+' | '.join(vals)+' |' for vals in values]
    markup=['<div class="table-scroll"><table><caption>'+html.escape(caption)+'</caption><thead><tr>']
    markup += ['<th>'+html.escape(head)+'</th>' for head in heads];markup += ['</tr></thead><tbody>']
    markup += ['<tr>'+''.join('<td>'+('<span class="chip">'+html.escape(value)+'</span>' if head=='Class' else html.escape(value))+'</td>' for head,value in zip(heads,vals))+'</tr>' for vals in values]
    markup += ['</tbody></table></div>']
    return '\n'.join(md),'\n'.join(markup)


def render(data):
    n=data['round'];d=data['decision'];challenger=d['provisional_challenger'];status=d['status']
    url=f'https://github.com/aburan28/crypto/blob/main/research/ic_candidate_tournament_20260915/goal_20260924/improvement/round{n}'
    opening=(f"Round {n}: {status}. Selected challenger: {challenger}; retained/promoted winner: {d['winner']}. "
        f"{data['verified_runs']}/{data['total_runs']} native/profile pairs independently verified by the measured round's frozen checker. "
        'Fresh Linux transport replay is the evidence-PR merge gate. All results are limited to the registered synthetic toy panel.')
    md=[f'# Bounded IC round {n}', '', opening, '',
        'Primary metric: one supplied target, after reusable preparation through scalar replay; fixture generation is outside both timed algorithms. Cold instruction and native process costs are supplementary promotion gates. Each table uses one cost unit. Values are equal-cell geometric means of three-process per-point medians, with no target amortization.', '',
        'The IC reference is the qualified `pairinv` source. `rho` is the separately qualified cold-instruction reference; `rho_online` is the separately qualified online-time reference. Ratios to rho are descriptive. The K-instruction floor applies only to this full-rank collector; it is not a generic IC lower bound.', '',
        'The class column labels engineering experiments and accounting controls. No asymptotic advance is claimed. Variant names are readable aliases; the machine-readable export retains every canonical candidate, workload and run ID. `stop3` uses the legacy adaptive orbit bound, which is not a universal three-column guarantee. Actual admitted bases and columns below are authoritative.', '']
    markup=[f'<!-- BEGIN ic-bounded-round-{n}-20260925 -->',f'<section class="panel" id="ic-bounded-round-{n}-20260925">',
        f'<h2>Bounded IC round {n}: {html.escape(status)}</h2>',f'<p>{html.escape(opening)}</p>',
        '<p>Primary metric: single-target online native time, reusable preparation excluded and scalar replay included. Cold cost remains an additional acceptance gate. Engineering/accounting only; no global-optimum or asymptotic claim. '+
        f'<a href="{url}/README.md">Full report</a> · <a href="{url}/RESULTS.json">Frozen tables and stage diagnostics</a> · <a href="{url}/RUNS.csv">Every run, including failures</a>.</p>']
    for stage in ('confirmation','replay','development','smoke','selection','aa'):
        block=data['stages'][stage]
        md += ['## '+stage.title(),'',f"Verified {block['verified_runs']}/{block['runs']} pairs.",'']
        if stage in ('confirmation','replay','development','smoke'):
            markup += [f'<h3>{stage.title()}</h3>']
        for unit,title in (('online','Single-target online time'),('instructions','Complete cold instructions'),('cold','Complete cold native time')):
            table,html_table=tables(block['table'],unit,stage.title()+': '+title)
            md += ['### '+title,'',table,'']
            if stage in ('confirmation','replay','development','smoke'):markup += [html_table]
        if block.get('retained_portfolio'):
            md += ['Retained portfolio:','']+[f"- `{row['candidate']}`: {row['reason']}." for row in block['retained_portfolio']]+['']
        if block['failures']:
            md += ['Retained failures:','', '```json',json.dumps(block['failures'],indent=2),'```','']
        md += ['### Actual base and matrix sizes', '',
            '| Variant | Curve cell | Usable points B / folded columns K / final rank |',
            '|---|---|---|']
        for row in block['table']:
            if row['mode'] != 'ic':
                continue
            for cell,shapes in row['actual_B_columns_rank'].items():
                md += ['| '+row['alias']+' | '+cell+' | '+ '; '.join(' / '.join(fmt(v) for v in shape) for shape in shapes)+' |']
        md += ['']
    md += ['## Confirmation rule and scope','','The selected challenger must pass both final stages: cold Ir and native ratios <=0.8, online ratio <=1, all three predeclared one-sided familywise upper bounds <1, and each metric\'s largest per-cell ratio <=1.1. The nominal familywise rule allocates alpha across three attempts, two stages and three metrics. Bootstrap coverage is approximate, and replay repeats the same targets in new processes. Ordinary 95% table intervals are descriptive; they do not replace this gate.','']
    for stage in ('confirmation','replay'):
        result=d.get(stage) or {};metrics=(result.get('familywise') or {}).get('metrics',{})
        md += ['### '+stage.title()+' promotion evidence','','| Metric | Candidate / reference | Familywise upper | Largest cell ratio |','|---|---|---|---|']
        for metric in ('online_ns','instructions','cold_ns'):
            row=metrics.get(metric,{});cell=row.get('per_cell') or {}
            md += ['| '+' | '.join([metric,fmt(row.get('ratio')),fmt(row.get('upper')),fmt(max(cell.values()) if cell else None)])+' |']
        md += ['']
    md += ['Decision reasons:','']+(['- '+reason for reason in d['reasons']] or ['- Every declared promotion gate passed.'])+['',
        'Complete run keys, admitted B/column counts, exclusive instruction ledgers, online phase clocks, peak process RSS and collection statistics are retained in `RUNS.csv` and `RESULTS.json`; raw certificates, sources, binaries and profiles are in the archive. Base-only memory was not instrumented and remains unknown. Yield intervals resample whole target/walk clusters, never individual dependent queries; few-target diagnostic intervals have weak coverage.','',
        'F4/F5/SAT and large sparse LA backends were not executed in this round. The tested mechanisms are pair-table PDP policies, Frobenius-orbit base construction and scalar-field Gaussian row kernels. No result here establishes a production-curve crossover or a globally fastest IC method.','']
    markup += ['<p>Complete cold Ir is user-space guest instructions, not curve additions. The K-instruction collector floor is implementation-specific and inapplicable to rho. Engineering rows test implementation changes; accounting rows are references or compatibility controls. Per-target identities, exact admitted bases, phase ledgers, uncertainty, retained failures and the final gate are available in the linked frozen report.</p>', '</section>',f'<!-- END ic-bounded-round-{n}-20260925 -->']
    return '\n'.join(md),'\n'.join(markup)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--markdown',type=Path,required=True)
    p.add_argument('--html',type=Path)
    p.add_argument('--scoreboard',type=Path)
    p.add_argument('--verify',action='store_true')
    a=p.parse_args()
    data=json.loads(a.results.read_text());md,markup=render(data)
    if a.verify:
        if a.markdown.read_text() != md:
            raise ValueError('Markdown tables differ from the frozen export')
        if a.html and a.html.read_text() != markup+'\n':
            raise ValueError('HTML tables differ from the frozen export')
        if a.scoreboard:
            page=a.scoreboard.read_text()
            start=f'<!-- BEGIN ic-bounded-round-{data["round"]}-20260925 -->'
            end=f'<!-- END ic-bounded-round-{data["round"]}-20260925 -->'
            if page.count(start)!=1 or page.count(end)!=1 or page.split(start)[1].split(end)[0] != markup.split(start)[1].split(end)[0]:
                raise ValueError('Scoreboard section differs from the frozen export')
    else:
        if not a.html or a.scoreboard:
            p.error('Rendering requires --html; --scoreboard is verification only')
        a.markdown.write_text(md);a.html.write_text(markup+'\n')

if __name__=='__main__':main()
