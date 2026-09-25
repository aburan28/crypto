#!/usr/bin/env python3
"""Audit the bounded development sequence and freeze its static report/table."""
import collections
import html
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
LAB = HERE.parent
sys.path.insert(0, str(LAB))
from autolab import paired_wall_ratio
from tournament import digest, read, write

PANELS = ('implementation_v2','algebra','factor_base','union_followup','interactions','interaction_replication')


def main():
    panels = []
    for name in PANELS:
        root = HERE/name
        raw = subprocess.check_output([sys.executable,str(root/'evaluator/autolab.py'),
                                       'verify','--round',str(root)],text=True)
        audit = json.loads(raw)
        print(name, raw.strip(), flush=True)
        contract, summary = read(root/'contract.json'), read(root/'summary.json')
        receipts = [read(p) for p in sorted((root/'trials').glob('*/*/rep-*/receipt.json'))]
        table = []
        for row in summary['table']:
            row = dict(row)
            row['paired_over_incumbent'] = paired_wall_ratio(contract,receipts,row['arm'])
            row['paired_over_rho'] = {f'rho_w{w}':paired_wall_ratio(contract,receipts,row['arm'],f'rho_w{w}')
                                      for w in (1,8,32)}
            # Rank/support diagnostics for partial arms remain visible by cell;
            # incomplete workloads still have no aggregate cost or comparison.
            row['verified_rank_by_case'] = {
                c['id']: sorted({r['certificate']['rank'] for r in receipts
                    if r['case']==c['id'] and r['arm']==row['arm'] and r['status']=='VERIFIED'
                    and r['certificate']['rank'] is not None}) for c in contract['cases']}
            table.append(row)
        panels.append({'panel':name,'audit':audit,'seed':contract['seed'],
            'contract_sha256':digest(root/'contract.json'),'binary_sha256':contract['pinned']['worker'],
            'cells':sorted({f"{c['job']['degree']}a{c['job']['curve_a']}" for c in contract['cases']}),
            'cases':len(contract['cases']),'repetitions':contract['repetitions'],
            'aa_paired_ratio':paired_wall_ratio(contract,receipts,'aa_control'),
            'verified':sum(r['status']=='VERIFIED' for r in receipts),'trials':len(receipts),
            'failures':dict(collections.Counter(r.get('reason') for r in receipts if r['status']!='VERIFIED')),
            'table':table})
    failed = HERE/'implementation'
    # Its evaluator's summary has the preserved null-rank bug; verify its raw
    # receipts through the same frozen audit, without inventing a summary.
    code = ('import sys,json; from pathlib import Path; '
            'sys.path.insert(0,sys.argv[1]); import autolab; '
            'c,rows=autolab.audit(Path(sys.argv[2]),complete=False); '
            'print(json.dumps({"status":"REPORTING_FAILURE_RETAINED",'
            '"raw_receipts_verified":len(rows),"verified_jobs":sum(r["status"]=="VERIFIED" for r in rows)}))')
    failure = json.loads(subprocess.check_output([sys.executable,'-c',code,str(failed/'evaluator'),str(failed)],text=True))
    result = {'status':'DEVELOPMENT_ONLY','classification':'accounting and integration engineering',
        'promotion':False,'global_optimum_claim':False,'panels':panels,'retained_first_run':failure,
        'trials':sum(p['trials'] for p in panels),'verified':sum(p['verified'] for p in panels),
        'unit':'complete cold child wall seconds','operation_counts':None,'S':None,
        'ratio_definition':'Equal-cell geometric mean of per-case median ratios; repetitions are not independent targets.',
        'limits':['No pinned affinity or measured memory cap on native path.',
                  'Public development targets only; no final confirmation used.',
                  'These use the general worker, not the archived optimized tiny-IC winner.',
                  'A/A noise prevents interpreting small native differences as improvements.',
                  'One-column IC matrices remain outside this checker protocol, not mathematically disproved.',
                  'Implementation_v2 overlapped a 1.11-second local rho test run; its timing is diagnostic only.',
                  'Other host activity is uncontrolled; no runtime promotion or cross-host claim.'],
        'reporter_sha256':digest(Path(__file__)),'ratio_implementation_sha256':digest(LAB/'autolab.py')}
    write(HERE/'RESULTS.json',result)
    lines=['# IC autolab development checkpoint — 2026-09-24','',
        'The portable loop completed six frozen, audited development screens. '
        f"**{result['verified']}/{result['trials']} jobs verified**; every rejection is retained. "
        'No performance winner is promoted. The first 192-job run completed its solves but failed '
        'during rho null-rank reporting; its original evaluator and receipts remain unchanged.','',
        'All jobs include curve/target generation, base construction, unsuccessful attempts, solving, '
        'relation processing, scalar-field linear algebra, descent, verification and output. '
        'The independent Python audit is outside the child. Operation counts, S and floor ratios '
        'are unmeasured, not zero. Seconds do not substitute for them.','',
        '## Sequence and findings','',
        '| Panel | Cells | Verified / jobs | A/A paired wall ratio |',
        '|---|---|---:|---:|']
    for p in panels:
        lines.append(f"| {p['panel']} | {', '.join(p['cells'])} | {p['verified']} / {p['trials']} | {p['aa_paired_ratio']:.4f} |")
    lines += ['',
        '- The implementation grid tests batch/window/linear-algebra interactions as complete jobs.',
        '- All seven PDP backends, each with dense and sparse relation linear algebra, passed the degree-9 panel.',
        '- Initial two-generator Frobenius unions failed the existing nontrivial-matrix admission rule. '
        'That is a protocol rejection, not evidence that a one-column IC method is impossible.',
        '- Enlarging the union to seed masks `[1,2,4,8]` passed both cells. Other unions remained partly or wholly rejected.',
        '- That union, its batch-1 and dense-LA combinations, and their individual parents passed the degree-13/17 '
        'interaction panel and its separate-seed development replication.',
        '- Rho controls request 1, 8 and 32 walks; the implementation can cap effective widths on tiny groups. '
        'The strongest general reference is not established by this small timing panel.',
        '- The native A/A controls expose uncontrolled timing variation. The local rho test overlap in '
        '`implementation_v2` is recorded explicitly. No native speedup claim follows from these screens.','',
        '## Complete static measurements','',
        'One unit: complete cold child milliseconds. Wall-time columns are diagnostics; ratios pair identical '
        'cases and weight curve cells equally. Repetitions are not new targets. Each rho ratio is reported '
        'separately, avoiding a post-hoc claim about the best reference. Incomplete arms have no aggregate '
        'cost or ratio. Their per-case admitted ranks and failure reasons remain in `RESULTS.json`.','']
    fragment=['<!-- BEGIN ic-autolab-20260924 -->','<section id="ic-autolab-20260924">',
        '<h2>IC autolab — native development checkpoint, 2026-09-24</h2>',
        f'<p>{result["verified"]}/{result["trials"]} verified jobs across six audited development screens; '
        '72 retained rejections and an earlier reporting failure are preserved. No promoted winner, '
        'operation-count result, ECC2K-130 gain or global-optimum claim. '
        'Native wall times have uncontrolled affinity, OS caches and host activity; A/A variation is recorded.</p>',
        '<p>Source: <a href="https://github.com/aburan28/crypto/blob/main/research/ic_candidate_tournament_20260915/autolab_20260924/RESULTS.json">'
        'frozen measurements</a> and <a href="https://github.com/aburan28/crypto/blob/main/research/ic_candidate_tournament_20260915/autolab_20260924/RESULTS.md">'
        'report and limits</a>. S and ratios to an operation floor are unmeasured. '
        'Classification: accounting and integration engineering.</p>']
    def fmt(v):
        return '—' if v is None else f'{v:.4f}'
    for p in panels:
        lines += [f"### {p['panel']}",'',
            '| Variant | Cold ms | Paired / incumbent | / rho 1 | / rho 8 | / rho 32 | Verified |',
            '|---|---:|---:|---:|---:|---:|---:|']
        fragment += [f'<h3>{p["panel"]}</h3>', '<table><thead><tr><th>Variant</th><th>Cold ms</th>'
                     '<th>Paired / incumbent</th><th>/ rho 1</th><th>/ rho 8</th><th>/ rho 32</th><th>Verified</th></tr></thead><tbody>']
        for r in p['table']:
            values=[r['arm'],fmt(r['median_cold_wall_seconds']*1000 if r['median_cold_wall_seconds'] else None),
                fmt(r['paired_over_incumbent'])]+[fmt(r['paired_over_rho'][f'rho_w{w}']) for w in (1,8,32)]+[f"{r['verified']}/{r['scheduled']}"]
            lines.append('| '+' | '.join(values)+' |')
            fragment.append('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in values)+'</tr>')
        lines.append('')
        fragment.append('</tbody></table>')
    lines += ['## Validation and continuation','',
        'Release tests: 62 descent-related tests passed; 36 rho-related tests passed and 3 expensive tests '
        'remained ignored. The Python suites have 50 tournament/autolab tests (2 Linux-only capability skips) '
        'and 14 boundary tests. All four repo skills validate. See `validation/` for exact output.','',
        'Restore raw evidence with `python3 research/ic_candidate_tournament_20260915/evidence/restore.py '
        '--archive autolab-20260924`; audit each completed panel using its frozen `evaluator/autolab.py verify --round PATH`. '
        'The restored artifact checksum does not require executing its archived binary.','',
        'Next: restore the archived optimized IC winner; qualify a strong matched rho reference; run a new '
        'Linux-amd64/Valgrind instruction tournament with fresh, disjoint confirmation points, '
        'the retained portfolio and the combination candidates. For new F4/F5/SAT implementation changes, '
        'also run the applicable frozen algebra regression and fresh algebra holdouts. '
        'Keep the one-column policies in the archive for a separately validated admission protocol.','']
    (HERE/'RESULTS.md').write_text('\n'.join(lines))
    fragment += ['</section>','<!-- END ic-autolab-20260924 -->']
    board=LAB.parents[1]/'docs/index-calculus-scoreboard.html'
    page=board.read_text()
    first,last='<!-- BEGIN ic-autolab-20260924 -->','<!-- END ic-autolab-20260924 -->'
    if first in page:
        begin,end=page.index(first),page.index(last)+len(last)
        page=page[:begin]+'\n'.join(fragment)+page[end:]
    else:
        if '</body>' not in page: raise ValueError('scoreboard lacks body close')
        page=page.replace('</body>','\n'.join(fragment)+'\n</body>',1)
    board.write_text(page)


if __name__=='__main__':
    main()
