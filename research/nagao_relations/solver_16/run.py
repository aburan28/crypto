"""solver_16 campaign driver.  Executes contract.json; appends one record per
solver cell to raw.jsonl (resumable: cells already present are skipped) and
writes summary.json at the end.  Nothing here changes the encodings; they are
weil.py, frozen with the contract.

    python3 run.py               # full panel
    python3 run.py --only cms    # one solver family (cms | wdsat | audit)
"""
import argparse
import hashlib
import json
import platform
import random
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import weil      # noqa: E402
import curves    # noqa: E402
import field     # noqa: E402

N = 131
CACHE = Path('/tmp/ec-baseline-cache')
WDSAT_SRC = CACHE / 'vendor/WDSat/src'
WORK = Path('/tmp/solver16-campaign')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def buildWdsat(configText, dest):
    dest.mkdir(parents=True, exist_ok=True)
    src = dest / 'src'
    src.mkdir(exist_ok=True)
    for p in WDSAT_SRC.iterdir():
        if p.suffix in ('.c', '.h'):
            (src / p.name).write_bytes(p.read_bytes())
    (src / 'config.h').write_text(configText)
    exe = dest / 'wdsat_solver'
    if not exe.exists():
        cmd = ['gcc', '-O3', '-w'] + sorted(str(p) for p in src.glob('*.c')) + ['-lm', '-o', str(exe)]
        subprocess.run(cmd, check=True, capture_output=True)
    return exe


def runWdsat(exe, inp, d, args, budget, nvars):
    t0 = time.time()
    try:
        p = subprocess.run([str(exe), '-i', str(inp), '-n', str(N), '-l', str(d), '-m', '3'] + args,
                           capture_output=True, text=True, timeout=budget)
    except subprocess.TimeoutExpired:
        return {'status': 'TIMEOUT', 'wall': time.time() - t0, 'decisions': None, 'witnesses': [], 'stdout_tail': None}
    wall = time.time() - t0
    lines = [s for s in p.stdout.splitlines() if s.strip()]
    witnesses = [s for s in lines if len(s) == nvars and set(s) <= {'0', '1'}]
    partial = [s for s in lines if len(s) == nvars and set(s) <= {'0', '1', '2'} and '2' in s]
    decisions = None
    for s in reversed(lines):
        t = s.strip()
        if t.startswith('conf:'):
            t = t[5:]
        if t.isdigit():
            decisions = int(t)
            break
    if any('UNSAT on XORGAUSS init' in s for s in lines):
        status = 'UNSAT_LINEAR'
        decisions = 0
    elif witnesses:
        status = 'SAT'
    elif partial:
        status = 'FALSE_SAT_PARTIAL_ASSIGNMENT'
    elif any(s.strip() == 'UNSAT' for s in lines):
        status = 'UNSAT'
    else:
        status = 'ERROR'
    return {'status': status, 'wall': wall, 'decisions': decisions, 'witnesses': witnesses,
            'returncode': p.returncode, 'stderr': p.stderr[-500:] if p.stderr else '',
            'stdout_tail': '\n'.join(lines[-3:])[:300]}


def verifyAll(pb, R, d, bitstrings, planted):
    out = []
    plantedFound = False
    for s in bitstrings:
        xs = weil.decodeBlocks(s, d)
        v = weil.verifyProjected(pb, R, xs)
        v['xs'] = xs
        v['genuine'] = bool(v['relation'] and not v['degenerate'])
        if planted is not None and sorted(xs) == sorted(planted):
            plantedFound = True
        out.append(v)
    return out, plantedFound


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', choices=['cms', 'wdsat', 'audit'])
    args = ap.parse_args()
    contract = json.loads((HERE / 'contract.json').read_text())
    panel = contract['panel']
    budget = panel['budget_seconds_per_solve']
    poly, terms = curves.findIrreduciblePoly(N)
    pb = field.Pb(N, poly)
    assert pb.isIrreducible() and curves.curveOrder(N) == 4 * weil.CHALLENGE_R
    WORK.mkdir(exist_ok=True)
    import pycryptosat
    provenance = {
        'contract_sha256': sha256(HERE / 'contract.json'), 'weil_sha256': sha256(HERE / 'weil.py'),
        'run_sha256': sha256(__file__), 'wdsat_source_sha256': {p.name: sha256(p) for p in sorted(WDSAT_SRC.iterdir()) if p.suffix in ('.c', '.h')},
        'pycryptosat': pycryptosat.__version__, 'python': platform.python_version(), 'platform': platform.platform(),
        'field_poly_terms': list(terms), 'git_head': subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=HERE, capture_output=True, text=True).stdout.strip(),
    }
    raw = HERE / 'raw.jsonl'
    done = set()
    if raw.exists():
        for line in raw.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                done.add(r['cell'])
    stream = raw.open('a')

    def emit(rec):
        rec['recorded_at'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        stream.write(json.dumps(rec, sort_keys=True) + '\n')
        stream.flush()
        done.add(rec['cell'])
        print(json.dumps({k: rec.get(k) for k in ('cell', 'status', 'wall', 'conflicts', 'decisions', 'genuine_relations', 'planted_found')}), flush=True)

    allD = sorted(set(panel['d_values_cms_s4']) | set(panel['d_values_cms_rr']) | set(panel['d_values_wdsat_plain_s4']) | set(panel['d_values_wdsat_xg_audit_s4']))
    rng = random.Random(panel['seed'])
    targets = {}
    for d in allD:
        A = len(weil.admissibleAbscissae(pb, d))
        rows = []
        for i in range(panel['targets_per_d']['uniform']):
            R = weil.uniformTarget(pb, rng)
            rows.append({'stratum': 'uniform', 'index': i, 'R': [str(R[0]), str(R[1])], 'planted': None})
        for i in range(panel['targets_per_d']['planted']):
            R, xs, _ = weil.plantedTarget(pb, d, rng)
            rows.append({'stratum': 'planted', 'index': i, 'R': [str(R[0]), str(R[1])], 'planted': xs})
        targets[d] = {'admissible': A, 'pairs': weil.pairEnumerationCount(A), 'rows': rows}
    (HERE / 'targets.json').write_text(json.dumps({'seed': panel['seed'], 'targets': targets, 'provenance': provenance}, indent=1) + '\n')

    systems = {}

    def system(enc, d, row):
        key = (enc, d, row['stratum'], row['index'])
        if key not in systems:
            R = (int(row['R'][0]), int(row['R'][1]))
            t0 = time.time()
            s = weil.buildS4(pb, d, R[0]) if enc == 's4' else weil.buildRR(pb, d, R[0], R[1])
            s.genSeconds = time.time() - t0
            systems[key] = s
        return systems[key]

    base = {'provenance': provenance}

    # ---------------- CryptoMiniSat, both encodings
    if args.only in (None, 'cms'):
        for enc, ds in (('s4', panel['d_values_cms_s4']), ('rr', panel['d_values_cms_rr'])):
            for d in ds:
                for row in targets[d]['rows']:
                    cell = 'cms/%s/d%d/%s%d' % (enc, d, row['stratum'], row['index'])
                    if cell in done:
                        continue
                    s = system(enc, d, row)
                    R = (int(row['R'][0]), int(row['R'][1]))
                    st = s.stats()
                    res = weil.cmsSolve(s, d, budget, enumerate=(row['stratum'] == 'planted'))
                    ver, plantedFound = verifyAll(pb, R, d, res['witnesses'], row['planted'])
                    bracket = None
                    if res['status'] != 'TIMEOUT' and res['wall'] < 30:
                        # bracket the FIRST decision (single solve), not the enumeration
                        bracket = weil.cmsConflictBracket(s, d, budget)
                    rec = dict(base, cell=cell, solver='cms', encoding=enc, d=d, stratum=row['stratum'], target_index=row['index'],
                               R=row['R'], planted=row['planted'], admissible=targets[d]['admissible'], pairs=targets[d]['pairs'],
                               system_stats=st, gen_seconds=round(s.genSeconds, 2), status=res['status'], wall=round(res['wall'], 3),
                               witnesses=res['witnesses'], verification=ver, genuine_relations=sum(1 for v in ver if v['genuine']),
                               planted_found=plantedFound if row['stratum'] == 'planted' else None,
                               conflicts=(bracket['decided_at'] if bracket else None), conflict_bracket=bracket,
                               cnf={k: res[k] for k in ('cnf_vars', 'clauses', 'xors')}, budget_seconds=budget)
                    emit(rec)

    # ---------------- WDSat plain (-b) on s4: brute-force tree control
    if args.only in (None, 'wdsat'):
        for d in panel['d_values_wdsat_plain_s4']:
            for mode in ('first', 'enum'):
                for row in targets[d]['rows']:
                    if mode == 'enum' and row['stratum'] == 'uniform':
                        continue   # a single decision already exhausts an UNSAT instance
                    cell = 'wdsat_plain/s4/d%d/%s%d/%s' % (d, row['stratum'], row['index'], mode)
                    if cell in done:
                        continue
                    s = system('s4', d, row)
                    R = (int(row['R'][0]), int(row['R'][1]))
                    st = s.stats()
                    cfg, cfgVals = weil.wdsatConfig(st, findAll=(mode == 'enum'))
                    exe = buildWdsat(cfg, WORK / ('s4-d%d-%s' % (d, mode)))
                    inp = WORK / ('s4-d%d-%s%d.anf' % (d, row['stratum'], row['index']))
                    if not inp.exists():
                        inp.write_text(s.anfText())
                    res = runWdsat(exe, inp, d, ['-b'], budget, s.nvars)
                    ver, plantedFound = verifyAll(pb, R, d, res['witnesses'], row['planted'])
                    rec = dict(base, cell=cell, solver='wdsat_plain', args=['-b'], mode=mode, encoding='s4', d=d, stratum=row['stratum'],
                               target_index=row['index'], R=row['R'], planted=row['planted'], admissible=targets[d]['admissible'],
                               pairs=targets[d]['pairs'], system_stats=st, wdsat_config=cfgVals, wdsat_binary_sha256=sha256(exe),
                               status=res['status'], wall=round(res['wall'], 3), decisions=res['decisions'], witnesses=res['witnesses'],
                               verification=ver, genuine_relations=sum(1 for v in ver if v['genuine']),
                               planted_found=plantedFound if row['stratum'] == 'planted' else None,
                               full_tree_over_sym='%d' % ((1 << (3 * d)) // 6), stdout_tail=res.get('stdout_tail'), budget_seconds=budget)
                    emit(rec)

    # ---------------- WDSat -x completeness audit on planted s4 (never a cost)
    if args.only in (None, 'audit'):
        for d in panel['d_values_wdsat_xg_audit_s4']:
            for row in targets[d]['rows']:
                if row['stratum'] != 'planted':
                    continue
                for xargs in (['-x', '-b'], ['-x']):
                    cell = 'wdsat_xg_audit/s4/d%d/planted%d/%s' % (d, row['index'], ''.join(a.strip('-') for a in xargs))
                    if cell in done:
                        continue
                    s = system('s4', d, row)
                    R = (int(row['R'][0]), int(row['R'][1]))
                    cfg, cfgVals = weil.wdsatConfig(s.stats(), findAll=False)
                    exe = buildWdsat(cfg, WORK / ('s4-d%d-first' % d))
                    inp = WORK / ('s4-d%d-planted%d.anf' % (d, row['index']))
                    if not inp.exists():
                        inp.write_text(s.anfText())
                    res = runWdsat(exe, inp, d, xargs, budget, s.nvars)
                    ver, plantedFound = verifyAll(pb, R, d, res['witnesses'], row['planted'])
                    found = res['status'] == 'SAT' and any(v['relation'] for v in ver)
                    rec = dict(base, cell=cell, solver='wdsat_xg_audit', args=xargs, encoding='s4', d=d, stratum='planted',
                               target_index=row['index'], R=row['R'], planted=row['planted'], status=res['status'], wall=round(res['wall'], 3),
                               decisions_not_a_cost=res['decisions'], witnesses=res['witnesses'], verification=ver,
                               genuine_relations=sum(1 for v in ver if v['genuine']), planted_found=plantedFound,
                               completeness=('ok' if found else 'FAILED: decomposable target not found'), budget_seconds=budget)
                    emit(rec)
    stream.close()
    print('done', flush=True)


if __name__ == '__main__':
    main()
