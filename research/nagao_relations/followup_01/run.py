"""Run the frozen follow-up panel; preserve completed cases on interruption."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ROOT / 'ecc2k130/codegen'
sys.path.insert(0, str(CODE))
import curves
import field
import nagaocompare
import nagaoannihilator


def main():
    contract = json.loads((HERE / 'contract.json').read_text())
    raw = HERE / 'raw.jsonl'
    output = HERE / 'summary.json'
    if raw.exists() or output.exists():
        raise SystemExit('Evidence exists; use a new run directory for another campaign.')
    sources = [Path(__file__), HERE / 'contract.json'] + [CODE / name for name in
        ('nagaocompare.py', 'nagaodecomp.py', 'nagaoannihilator.py', 'decomp.py',
         'field.py', 'curves.py', 'indexcalc.py', 'ir.py', 'build.py', 'cnf.py')]
    provenance = {
        'command': [sys.executable] + sys.argv,
        'git_commit': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'git_status': subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True),
        'python': platform.python_version(), 'platform': platform.platform(),
        'pycryptosat': importlib.metadata.version('pycryptosat'),
        'source_sha256': {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    }
    rows = []
    with raw.open('x') as handle:
        def record(row):
            rows.append(row)
            handle.write(json.dumps(row)+'\n')
            handle.flush()
        record({'kind':'provenance', **provenance})
        for panel in contract['coefficient_panels']:
            started = time.perf_counter()
            targetCoords = None
            if panel['targets'] != 'all':
                f = field.Onb(panel['n'])
                c = curves.Curve(f)
                targets = random.Random(contract['seed']).sample(nagaocompare.affinePoints(f,c),panel['targets'])
                targetCoords = [[f.toCoords(v) for v in p] for p in targets]
            try:
                result = nagaoannihilator.runExperiment(panel['n'],panel['dimension'],targetCoords)
                record({'kind':'coefficients','panel':panel,'result':result,
                        'diagnostic_seconds':time.perf_counter()-started})
                print('coefficient panel', panel, 'solutions', result['verifiedSolutions'], flush=True)
            except Exception:
                record({'kind':'error','panel':panel,'traceback':traceback.format_exc()})
        for panel in contract['sat_panels']:
            f = field.Onb(panel['n'])
            c = curves.Curve(f)
            started = time.perf_counter()
            table, oracle = nagaocompare.oracleTable(f,c,panel['weight'])
            targets = random.Random(contract['seed']).sample(nagaocompare.affinePoints(f,c),panel['targets'])
            record({'kind':'oracle','panel':panel,'oracle':oracle,
                    'diagnostic_seconds':time.perf_counter()-started})
            for target in targets:
                for variant in contract['sat_variants']:
                    try:
                        result = nagaocompare.compareInstance(f,c,target,panel['weight'],variant,
                            table.get(target,set()),contract['seconds_per_instance'],contract['maximum_projected_models'])
                        record({'kind':'sat','result':result})
                        print('SAT',panel['n'],variant,[f.toCoords(v) for v in target],
                              'complete',result['complete'],'solutions',len(result['projected_solutions']),flush=True)
                    except Exception:
                        record({'kind':'error','panel':panel,'variant':variant,
                                'target':[f.toCoords(v) for v in target], 'traceback':traceback.format_exc()})
    summary = {'provenance':provenance,'coefficient_panels':[],'sat_panels':[],
               'errors':[r for r in rows if r['kind']=='error'],
               'limits':contract['inference']}
    for row in rows:
        if row['kind']=='coefficients':
            r=row['result']
            summary['coefficient_panels'].append({k:r[k] for k in
                ('fieldDegree','subspaceDimension','targetCount','factorBasePoints','verifiedSolutions',
                 'targetHits','allSolutionSetsEqual','nagaoOperations','oracleOperations')})
    for panel in contract['sat_panels']:
        for variant in contract['sat_variants']:
            group=[r['result'] for r in rows if r['kind']=='sat' and
                   r['result']['field_degree']==panel['n'] and r['result']['variant']==variant]
            summary['sat_panels'].append({'n':panel['n'],'weight':panel['weight'],'variant':variant,
                'attempted':len(group),'complete':sum(r['complete'] for r in group),
                'correct_completed':sum(r['correct'] for r in group),
                'errors':sum(bool(r['errors'] or r['extra'] or r['missing']) for r in group),
                'projected_solutions_completed':sum(len(r['projected_solutions']) for r in group if r['complete']),
                'projected_solutions_partial':sum(len(r['projected_solutions']) for r in group if not r['complete']),
                'diagnostic_seconds':sum(sum(r['diagnostic_seconds'].values()) for r in group)})
    with output.open('x') as handle:
        json.dump(summary,handle,indent=2)
        handle.write('\n')
    print(json.dumps(summary['sat_panels']))


if __name__ == '__main__':
    main()
