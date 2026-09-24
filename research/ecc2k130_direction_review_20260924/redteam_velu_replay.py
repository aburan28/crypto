#!/usr/bin/env python3
"""Independent, bounded full-point replay of saved twist-kernel representatives.

This imports no production arithmetic. It uses bit-polynomial multiplication,
Fermat inversion, explicit arithmetic in F_q[s]/(s^2+s+1), and direct paired
Vélu sums. It is an audit prototype, not a production transport interface.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


class Field:
    def __init__(self, n, irr):
        self.n, self.irr = n, irr
        self.counts = dict(mul=0, sqr=0, inv=0)

    def mul(self, x, y):
        self.counts['mul'] += 1
        z = 0
        while y:
            if y & 1:
                z ^= x
            y >>= 1
            x <<= 1
            if x >> self.n:
                x ^= self.irr
        return z

    def sqr(self, x):
        self.counts['sqr'] += 1
        return self.mul(x, x)

    def inv(self, x):
        self.counts['inv'] += 1
        assert x
        r, y, e = 1, x, (1 << self.n)-2
        while e:
            if e & 1:
                r = self.mul(r, y)
            y = self.sqr(y)
            e >>= 1
        assert self.mul(r, x) == 1
        return r

    def batch(self, xs):
        accum = [1]
        for x in xs:
            assert x
            accum.append(self.mul(accum[-1], x))
        y, out = self.inv(accum[-1]), [0]*len(xs)
        for i in range(len(xs)-1, -1, -1):
            out[i] = self.mul(y, accum[i])
            y = self.mul(y, xs[i])
        return out


class Replay:
    def __init__(self, field):
        self.F = field
        self.add_calls = 0

    def add(self, P, Q, a=0):
        self.add_calls += 1
        F = self.F
        if P is None:
            return Q
        if Q is None:
            return P
        x, y = P
        u, v = Q
        if x == u:
            if y ^ v == x or x == 0:
                return None
            lam = x ^ F.mul(y, F.inv(x))
            z = F.sqr(lam) ^ lam ^ a
            return z, F.sqr(x) ^ F.mul(lam ^ 1, z)
        lam = F.mul(y ^ v, F.inv(x ^ u))
        z = F.sqr(lam) ^ lam ^ x ^ u ^ a
        return z, F.mul(lam, x ^ z) ^ z ^ y

    def mul(self, P, k, a=0):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P, a)
            P = self.add(P, P, a)
            k >>= 1
        return R

    def on_curve(self, P, a=0, b=1):
        if P is None:
            return True
        F, (x, y) = self.F, P
        return F.sqr(y) ^ F.mul(x, y) == F.mul(x, F.sqr(x)) ^ F.mul(a, F.sqr(x)) ^ b

    @staticmethod
    def extension_add(x, y):
        return x[0] ^ y[0], x[1] ^ y[1]

    def extension_mul(self, x, y):
        # s²=s+1, so (a+b*s)(c+d*s)=(ac+bd)+(ad+bc+bd)*s.
        F = self.F
        ac, bd = F.mul(x[0], y[0]), F.mul(x[1], y[1])
        middle = F.mul(x[0] ^ x[1], y[0] ^ y[1]) ^ ac
        return ac ^ bd, middle

    def extension_sqr(self, x):
        a2, b2 = self.F.sqr(x[0]), self.F.sqr(x[1])
        return a2 ^ b2, b2

    def direct_velu(self, P, kernel):
        """Twist points (u,v) become original points (u,v+s*u)."""
        F, (x, y) = self.F, P
        Eadd, Emul, Esqr = self.extension_add, self.extension_mul, self.extension_sqr
        t, X, Y = 0, (x, 0), (y, 0)
        invs = F.batch([x ^ u for u, v in kernel])
        for (u, v), inverse in zip(kernel, invs):
            t ^= u
            pair = []
            for w in (v, v ^ u):
                lam = F.mul(y ^ w, inverse), F.mul(u, inverse)
                x3 = Eadd(Eadd(Esqr(lam), lam), (x ^ u, 0))
                y3 = Eadd(Eadd(Emul(lam, Eadd((x, 0), x3)), x3), (y, 0))
                pair.append((x3, y3))
            # Q and -Q have equal x; their y-coordinate sum is u.
            X = Eadd(X, Eadd(pair[0][0], pair[1][0]))
            Y = Eadd(Y, Eadd(Eadd(pair[0][1], pair[1][1]), (u, 0)))
        assert X[1] == Y[1] == 0, 'nonrational normalized image'
        return X[0], Y[0] ^ t


def point_json(P):
    return None if P is None else [hex(P[0]), hex(P[1])]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=root/'twist_torsion_results.json')
    parser.add_argument('--out', type=Path, default=root/'redteam_velu_replay.json')
    args = parser.parse_args()
    assert not args.out.exists(), 'preserve prior evidence; use a fresh --out'
    begin = time.perf_counter()
    source = json.loads(args.input.read_text())
    P = (int('051C99BFA6F18DE467C80C23B98C7994AA', 16), int('042EA2D112ECEC71FCF7E000D7EFC978BD', 16))
    Q = (int('06C997F3E7F2C66A4A5D2FDA13756A37B1', 16), int('04A38D11829D32D347BD0C0F584D546E9A', 16))
    ell, results = 263, []
    trace0, trace1 = 2, -1
    for _ in range(130):
        trace0, trace1 = trace1, -trace1 - 2*trace0
    for run in source['runs']:
        assert time.perf_counter()-begin < 60, 'review time budget'
        F = Field(131, int(run['field_modulus_hex'], 16))
        R = Replay(F)
        k, q = run['integer_parameters'], 1 << 131
        assert int(k['trace']) == trace1
        assert int(k['order']) == q+1-trace1
        assert int(k['twist_order']) == q+1+trace1
        assert int(k['tau131_A']) % ell == ell-1 and int(k['tau131_B']) % ell == 0
        assert (q+1+trace1) % ell**2 == 0 and (q+1+trace1) % ell**3 != 0
        assert R.on_curve(P) and R.on_curve(Q) and not R.on_curve(P, 1)
        representatives = []
        for saved in run['representative_x_maps']:
            assert time.perf_counter()-begin < 60, 'review time budget'
            G = tuple(int(z, 16) for z in saved['kernel_generator_on_twist'])
            assert G is not None and R.on_curve(G, 1) and R.mul(G, ell, 1) is None
            T, kernel, t = G, [], 0
            for _ in range(131):
                assert R.on_curve(T, 1)
                kernel.append(T)
                t ^= T[0]
                T = R.add(T, G, 1)
            assert {u for u, v in kernel} == {int(z, 16) for z in saved['kernel_abscissae']}
            destb = int(saved['codomain_b'], 16)
            assert destb == 1 ^ t ^ F.sqr(t)
            targets = {'P': P, 'Q': Q, 'P+Q': R.add(P, Q), '2P': R.add(P, P)}
            images = {name: R.direct_velu(point, kernel) for name, point in targets.items()}
            assert all(R.on_curve(point, 0, destb) for point in images.values())
            assert images['P'][0] == int(saved['image_P_x'], 16)
            assert images['Q'][0] == int(saved['image_Q_x'], 16)
            assert images['P+Q'] == R.add(images['P'], images['Q'])
            assert images['2P'] == R.add(images['P'], images['P'])
            minus_q = images['Q'][0], images['Q'][1] ^ images['Q'][0]
            assert images['P+Q'] != R.add(images['P'], minus_q), 'negative sign control'
            representatives.append({
                'kernel_line': saved['line'], 'orbit_length': saved['orbit_length'],
                'kernel_generator_on_twist': point_json(G), 'half_kernel_size': len(kernel),
                'half_kernel_sum_t': hex(t), 'codomain_b': hex(destb),
                'full_images': {name: point_json(point) for name, point in images.items()},
                'codomain_checks': 4, 'saved_x_checks': 2, 'full_additivity_checks': 2,
                'negative_sign_control': 'PASS', 'status': 'PASS',
            })
        results.append({'seed': run['seed'], 'representatives': representatives,
                        'integer_checks': 'PASS', 'native_field_counts': F.counts,
                        'native_group_add_calls': R.add_calls})
    count = sum(len(run['representatives']) for run in results)
    output = {
        'schema': 'ecc2k130-redteam-full-velu-replay-v1', 'status': 'PASS',
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'command_argv': [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        'host': platform.platform(), 'python': sys.version,
        'script_sha256': sha(__file__), 'input_sha256': sha(args.input),
        'source_preflight_sha256': sha(root/'twist_torsion_preflight.py'),
        'contract_sha256': sha(root/'twist_torsion_contract.md'),
        'arithmetic_dependencies': 'Python standard library only; no production arithmetic imported',
        'runs': results, 'elapsed_seconds': time.perf_counter()-begin,
        'checks': {'representative_kernels': count, 'kernel_abscissae': 131*count,
                   'full_point_images': 4*count, 'codomain_equations': 4*count,
                   'saved_x_matches': 2*count, 'full_additivity_equalities': 2*count,
                   'negative_sign_controls': count, 'disagreements': 0},
        'counter_scope': 'mul includes sqr and inv internals; counters are diagnostic native units, not calibrated work',
        'limitations': ['four representative lines per seed; not all 264 full-point maps',
                       'finite checks do not alone prove an arbitrary-input map theorem',
                       'kernel and infinity inputs are not implemented in direct_velu',
                       'no production transport integration, PDP advantage, relation matrix or DLP result',
                       'no comparison against rho or production arithmetic performance'],
    }
    args.out.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'status': output['status'], 'checks': output['checks'],
                      'elapsed_seconds': output['elapsed_seconds'], 'output': str(args.out)}, indent=2))


if __name__ == '__main__':
    main()
