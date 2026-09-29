"""Sizing pilot (DEVELOPMENT DATA, not campaign evidence): one uniform and one
planted target per d on the real ECC2K-130 curve, both encodings, to fix the
d range, budgets and WDSat static sizes before the contract is frozen."""
import json
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
from test_weil import buildWdsat  # noqa: E402

N = 131


def theField():
    poly, terms = curves.findIrreduciblePoly(N)
    pb = field.Pb(N, poly)
    assert pb.isIrreducible()
    order = curves.curveOrder(N)
    assert order == 4 * weil.CHALLENGE_R
    return pb, terms


def main():
    pb, terms = theField()
    print('field poly terms', terms, 'order = 4r verified', flush=True)
    rng = random.Random(20260920)
    out = {'field_poly_terms': terms, 'rows': []}
    work = Path('/tmp/solver16-pilot')
    work.mkdir(exist_ok=True)
    ds = [int(a) for a in sys.argv[1:]] or [4, 5, 6, 7]
    for d in ds:
        adm = None
        if d <= 20:
            adm = len(weil.admissibleAbscissae(pb, d))
        Ru = weil.uniformTarget(pb, rng)
        Rp, planted, _ = weil.plantedTarget(pb, d, rng)
        for stratum, R in (('uniform', Ru), ('planted', Rp)):
            for enc in ('s4', 'rr'):
                t0 = time.time()
                sysm = weil.buildS4(pb, d, R[0]) if enc == 's4' else weil.buildRR(pb, d, R[0], R[1])
                gen = time.time() - t0
                st = sysm.stats()
                cfg, cfgVals = weil.wdsatConfig(st, findAll=True)
                mem = weil.historyBytes(cfgVals)
                exe = buildWdsat(cfg, work / ('%s-d%d' % (enc, d)))
                inp = work / ('%s-d%d-%s.anf' % (enc, d, stratum))
                inp.write_text(sysm.anfText())
                for args in (['-x', '-b'], ['-x']):
                    t0 = time.time()
                    try:
                        p = subprocess.run([str(exe), '-i', str(inp), '-n', str(N), '-l', str(d), '-m', '3'] + args,
                                           capture_output=True, text=True, timeout=600)
                        wall = time.time() - t0
                        lines = [s for s in p.stdout.splitlines() if s.strip()]
                        status = 'UNSAT' if 'UNSAT' in p.stdout else 'SAT?'
                        sols = [s for s in lines if len(s) == sysm.nvars and set(s) <= {'0', '1'}]
                        conf = int(lines[-1]) if lines and lines[-1].strip().lstrip('-').isdigit() else None
                        warn = [s for s in lines if '!!!' in s]
                        verified = None
                        if sols:
                            verified = [weil.verifyProjected(pb, R, weil.decodeBlocks(s, d)) for s in sols]
                    except subprocess.TimeoutExpired:
                        wall, status, sols, conf, warn, verified = 600.0, 'TIMEOUT', [], None, [], None
                    row = {'d': d, 'stratum': stratum, 'encoding': enc, 'args': args, 'admissible': adm,
                           'pairs': weil.pairEnumerationCount(adm) if adm else None,
                           'gen_seconds': round(gen, 2), 'stats': st, 'history_mb': round(mem['total'] / 1e6, 1),
                           'status': status, 'solutions': len(sols), 'conflicts': conf, 'wall': round(wall, 3),
                           'warnings': warn, 'verified': verified, 'planted': planted if stratum == 'planted' else None}
                    out['rows'].append(row)
                    print(json.dumps({k: row[k] for k in ('d', 'stratum', 'encoding', 'args', 'pairs', 'gen_seconds',
                                                            'history_mb', 'status', 'solutions', 'conflicts', 'wall', 'warnings')}),
                          flush=True)
                    if verified:
                        print('   verified:', [(v['relation'], v['degenerate']) for v in verified], flush=True)
        (HERE / 'pilot.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
