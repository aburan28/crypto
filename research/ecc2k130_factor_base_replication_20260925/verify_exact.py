#!/usr/bin/env python3
"""Independent exact-coordinate replay of the bounded n=131 smoke."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/ecc2k130_relations"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_oriented_transport_20260924"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_dual_transport_20260925"))
import verify as toy_verify
from fastfield import FastGF2m, IRR131
from relations import (Koblitz, CHALLENGE_PX, CHALLENGE_PY, CHALLENGE_QX,
                       CHALLENGE_QY, CHALLENGE_ELL)
from oriented_velu import BinaryVeluMap
from dual_transport import DualTransport


def from_hex(pair):
    return tuple(int(z, 16) for z in pair)


def sha_int(label):
    return int.from_bytes(hashlib.sha256(label.encode()).digest(), 'big')


def prefix(E):
    base, seen = [], set()
    for x in range(4096):
        for P in E.points_over(x):
            Q = E.mul(P, 4)
            if Q is not None and Q not in seen:
                seen.add(Q)
                base.append(Q)
                if len(base) == 16:
                    return base
    raise AssertionError('incomplete exact base')


def pairs(E, base):
    table = defaultdict(list)
    for i in range(16):
        for j in range(i, 16):
            table[E.add(base[i], base[j])].append([i, j])
    return table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    assert not args.out.exists()
    raw = json.loads(args.input.read_text())
    assert raw['status'] == 'PASS'
    preflight_path = ROOT / 'research/ecc2k130_direction_review_20260924/twist_torsion_results.json'
    assert hashlib.sha256(preflight_path.read_bytes()).hexdigest() == raw['preflight_sha256']
    preflight = json.loads(preflight_path.read_text())
    frozen = next(run for run in preflight['runs'] if run['seed'] == 20260924)
    listed = {tuple(x['line']): x for x in frozen['representative_x_maps']}
    V = from_hex(frozen['basis_on_twist'][1])
    F = FastGF2m(131, IRR131)
    E0, Twist = Koblitz(F, 0, 1), Koblitz(F, 1, 1)
    P, Q = (CHALLENGE_PX, CHALLENGE_PY), (CHALLENGE_QX, CHALLENGE_QY)
    assert E0.on_curve(P) and E0.on_curve(Q)
    assert E0.mul(P, CHALLENGE_ELL) is None
    assert E0.mul(Q, CHALLENGE_ELL) is None
    public = raw['public_challenge_input']
    assert public['P'] == [hex(P[0]), hex(P[1])]
    assert public['Q'] == [hex(Q[0]), hex(Q[1])]
    assert public['subgroup_order'] == str(CHALLENGE_ELL)
    assert public['repository_sha256'] == hashlib.sha256((ROOT / public['repository_path']).read_bytes()).hexdigest()
    checked = []
    for line in raw['lines']:
        key = tuple(line['line'])
        assert key in ((1, 0), (1, 4))
        item = listed[key]
        H = from_hex(item['kernel_generator_on_twist'])
        phi = BinaryVeluMap.from_generator(E0, Twist, H, 263)
        D = DualTransport(E0, Twist, H, V, 263, P)
        E1 = Koblitz(F, phi.codomain.a, phi.codomain.b)
        assert hex(E1.b) == line['codomain_b'] == item['codomain_b']
        assert D.compose(P) == E0.mul(P, 263)
        assert D.compose(Q) == E0.mul(Q, 263)
        B0, B1 = prefix(E0), prefix(E1)
        Bt = [phi(S) for S in B0]
        Bb = [E0.mul(D.dual(S), pow(263, -1, CHALLENGE_ELL)) for S in B1]
        bases = {'original': B0, 'transported': Bt,
                 'descendant_native': B1, 'pullback': Bb}
        assert {k: [from_hex(p) for p in line['bases'][k]] for k in bases} == bases
        assert all(phi(a) == b for a, b in zip(Bb, B1))
        assert all(E1.on_curve(phi(E0.frobenius(a))) for a in Bb)
        P1, Q1 = phi(P), phi(Q)
        streams = {'source': [], 'leaf': []}
        tag = f'{key[0]},{key[1]}'
        for i, pair in enumerate(line['coefficient_pairs']):
            u = sha_int(f'degree263-smoke-v1|{tag}|{i}|u') % CHALLENGE_ELL
            v = 1 + sha_int(f'degree263-smoke-v1|{tag}|{i}|v') % (CHALLENGE_ELL - 1)
            assert pair == [u, v]
            T0 = E0.add(E0.mul(P, u), E0.mul(Q, v))
            T1 = E1.add(E1.mul(P1, u), E1.mul(Q1, v))
            assert phi(T0) == T1
            streams['source'].append((u, v, T0))
            streams['leaf'].append((u, v, T1))
        for name, base in bases.items():
            leaf = name in ('transported', 'descendant_native')
            E = E1 if leaf else E0
            table = pairs(E, base)
            rows = []
            cases = line['results'][name]['cases']
            assert len(cases) == 32
            for (u, v, T), case in zip(streams['leaf' if leaf else 'source'], cases):
                witness = table.get(T, [])
                assert case['witness_count'] == len(witness)
                assert case['first_witness'] == (witness[0] if witness else None)
                old_rank = toy_verify.independent_rank(rows, CHALLENGE_ELL)[0]
                if witness:
                    a, b = witness[0]
                    row = [0] * 17
                    row[a] += 1
                    row[b] += 1
                    row[-1] = -v % CHALLENGE_ELL
                    rows.append(row)
                new_rank = toy_verify.independent_rank(rows, CHALLENGE_ELL)[0]
                assert case['independent'] == (new_rank > old_rank)
            assert line['results'][name]['independent_rank'] == toy_verify.independent_rank(rows, CHALLENGE_ELL)[0]
            assert line['results'][name]['hits'] == sum(c['witness_count'] > 0 for c in cases)
            assert line['results'][name]['misses'] == 32 - line['results'][name]['hits']
            assert line['results'][name]['pair_entries'] == sum(map(len, table.values())) == 136
        checked.append({'line': line['line'], 'targets': 32, 'arms': 4,
                        'hits': {name: line['results'][name]['hits'] for name in bases},
                        'ranks': {name: line['results'][name]['independent_rank'] for name in bases}})
    assert len(checked) == 2
    receipt = {'schema': 'ecc2k130-degree263-bounded-smoke-independent-replay-v1',
               'status': 'PASS', 'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
               'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'lines': checked, 'total_exact_target_replays': 64,
               'total_pair_case_replays': 256}
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': receipt['status'],
                      'target_replays': receipt['total_exact_target_replays'],
                      'pair_case_replays': receipt['total_pair_case_replays'],
                      'output': str(args.out)}, indent=2))


if __name__ == '__main__':
    main()
