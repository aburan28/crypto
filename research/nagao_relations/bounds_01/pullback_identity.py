"""Check a general-coefficient simplification for the next block experiment.

This is an algebraic identity test, not a new decomposition solver benchmark.
The deterministic test corpus is fixed in this source before execution.
"""
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'structured_01'))
import run as prior


def main():
    rng = random.Random(2026091431)
    checks = targetChecks = 0
    cases = []
    for n in (6, 18, 30):
        f = prior.field.Onb(n)
        alpha = prior.algebra.subfieldGenerator(f, 2)
        coefficients = [(a, b) for a in (0, f.one(), alpha, f.add(alpha, f.one()))
                        for b in (f.one(), alpha, f.add(alpha, f.one()))] if n == 6 else [(alpha, alpha)]
        vspace = prior.space.Space(f, [1 << i for i in range(4)])
        for aa, bb in coefficients:
            c = prior.algebra.Curve(f, aa, bb)
            targets = []
            if n == 6:
                for x in range(1, 64):
                    p = c.pointFromX(f.fromCoords(x))
                    if p:
                        targets += [p, c.neg(p)]
                for r, s in targets:
                    invR = f.inv(r)
                    dd = f.mul(f.add(r, s), invR)
                    if f.add(f.add(f.add(f.sqr(dd), dd), r), aa) != f.mul(bb, f.sqr(invR)):
                        raise ArithmeticError('curve-equation cancellation failed')
                    targetChecks += 1
            accepted = 0
            repetitions = 128 if n == 6 else 256
            for _ in range(repetitions):
                while True:
                    target = rng.choice(targets) if n == 6 else prior.previous.uniform(f, c, rng)
                    r, s = target
                    if not r:
                        continue
                    h2 = rng.choice(vspace.values)
                    roots = [b for b in c.asRoots(f.add(h2, r)) if b]
                    if roots:
                        break
                b = rng.choice(roots)
                z = rng.choice([x for x in vspace.values if x not in (0, r, h2)])
                u, t = f.add(h2, z), f.add(r, z)
                w = rng.choice(vspace.values)
                value = f.add(f.sqr(w), f.mul(u, w))
                invR, invB = f.inv(r), f.inv(b)
                gamma = f.mul(f.mul(t, invR), invB)
                eta = f.mul(f.mul(z, u), invB)
                dd = f.mul(f.add(r, s), invR)
                delta = f.add(f.add(r, f.mul(b, dd)), eta)
                a = f.add(f.mul(gamma, value), delta)
                lhs = f.add(f.mul(f.sqr(gamma), f.sqr(value)), f.mul(f.mul(z, invR), value))
                rhs = f.add(f.sqr(eta), f.mul(f.sqr(b), f.mul(bb, f.sqr(invR))))
                error = f.add(lhs, rhs)
                h, _ = prior.algebra.residualNorm(f, c, target, a, b)
                actualV = f.add(h[1], f.mul(z, u))
                if f.add(actualV, value) != error:
                    raise ArithmeticError('normalized support residual differs')
                if prior.solvers.scalar.polyEval(f, h, z) != f.mul(t, error):
                    raise ArithmeticError('conditioned-root residual differs')
                accepted += not error
                checks += 1
            cases.append({'n': n, 'a2': f.toCoords(aa), 'a6': f.toCoords(bb),
                          'identity_checks': repetitions, 'zero_residual_cases': accepted})
    result = {'scope': 'algebraic residual identities, not a solver speedup',
              'seed': 2026091431, 'identity_checks': checks,
              'exhaustive_nonzero_GF64_target_checks': targetChecks, 'cases': cases,
              'failures': 0, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()}
    with (HERE / 'pullback_identity_results.json').open('x') as out:
        out.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'cases'}))


if __name__ == '__main__':
    main()
