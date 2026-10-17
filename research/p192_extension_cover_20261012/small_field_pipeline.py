"""End-to-end genus-2 Jacobian relation control, separate from P-192.

Run with /Volumes/SSD990/cryptanalysis/sage -python.  This deliberately uses
the same-field cover H -> E over F_101, so its measurements validate pipeline
arithmetic and accounting, not a descended P-192 cover.
"""
import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, HyperellipticCurve, Integer, PolynomialRing, matrix, vector

HERE = Path(__file__).resolve().parent
SOURCE_SHA = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
SEED = 0x19220261014
rng = random.Random(SEED)
timing = {}


def timed(name, fn):
    start = time.monotonic_ns()
    result = fn()
    timing[name + '_ns'] = time.monotonic_ns() - start
    return result


def setup():
    F = GF(101)
    R = PolynomialRing(F, 'u')
    u = R.gen()
    E = EllipticCurve(F, [1, 46])
    assert E.cardinality() == 109 and Integer(109).is_prime()
    assert (u**3 + u + 46).is_irreducible()
    H = HyperellipticCurve(u**6 + u**2 + 46)
    assert H.genus() == 2 and (u**6 + u**2 + 46).is_squarefree()
    J = H.jacobian()
    assert J.cardinality() == 11663 == 109 * 107
    G = E(F(2), F(37))
    assert G.order() == 109
    D = J([u**2-G[0], R(G[1])])
    assert D != J(0) and 109*D == J(0)
    return F, R, u, E, H, J, G, D


F, R, u, E, H, J, G, D = timed('geometry_and_subgroup', setup)
FR = GF(109)


def build_base():
    points = []
    for x in F:
        y2 = x**6 + x**2 + 46
        if not y2.is_square():
            continue
        y = y2.sqrt()
        for v in sorted({y, -y}, key=int):
            P = J([u-x, R(v)])
            Q = 107*P
            if Q != J(0):
                assert 109*Q == J(0)
                points.append((P, Q, int(x), int(v)))
    rng.shuffle(points)
    distinct = {}
    for entry in points:
        distinct.setdefault(str(entry[1]), entry)
    base = list(distinct.values())[:32]
    assert len(base) == 32
    return base, len(points), len(distinct)


base, affine_points_nonzero_projection, distinct_projected = timed('factor_base', build_base)
B = len(base)


def build_pair_index():
    index = {}
    for i in range(B):
        for j in range(i, B):
            key = str(base[i][0] + base[j][0])
            index.setdefault(key, []).append((i, j))
    return index


pair_index = timed('pair_index', build_pair_index)


def find_three_sum(T):
    indices = list(range(B))
    rng.shuffle(indices)
    for k in indices:
        for i, j in pair_index.get(str(T-base[k][0]), []):
            if base[i][0] + base[j][0] + base[k][0] == T:
                return i, j, k
    return None


def collect_relations():
    rows, rhs, witnesses = [], [], []
    attempts = 0
    verified = 0
    rank = 0
    while attempts < 500 and rank < B:
        a = rng.randrange(1, 109)
        attempts += 1
        T = a*D
        decomp = find_three_sum(T)
        if decomp is None:
            continue
        i, j, k = decomp
        assert base[i][0] + base[j][0] + base[k][0] == T
        verified += 1
        row = [0]*B
        for t in (i, j, k):
            row[t] = (row[t]+1) % 109
        candidate = matrix(FR, rows + [row])
        if candidate.rank() > rank:
            rows.append(row)
            rhs.append((107*a) % 109)
            witnesses.append({'a': a, 'indices': [i, j, k]})
            rank += 1
    return rows, rhs, witnesses, attempts, verified, rank


rows, rhs, witnesses, attempts, verified, rank = timed('ordinary_relations_and_rank', collect_relations)
assert rank == B, f'rank={rank}, base={B}, attempts={attempts}'


def solve_logs():
    A = matrix(FR, rows)
    answer = A.solve_right(vector(FR, rhs))
    assert A*answer == vector(FR, rhs)
    for i in range(B):
        assert int(answer[i])*D == base[i][1]
    return [int(v) for v in answer]


logs = timed('final_linear_algebra_and_base_replay', solve_logs)

# The previously unseen target point is constructed outside either online
# interval.  The scalar is used only for the final independent fixture check.
target_scalar = rng.randrange(1, 109)
target = target_scalar*D
assert target != J(0)
target_digest = hashlib.sha256(str(target).encode()).hexdigest()


def target_online():
    for offset in range(109):
        decomp = find_three_sum(target + offset*D)
        if decomp is not None:
            result = (pow(107, -1, 109)*sum(logs[i] for i in decomp)-offset) % 109
            assert result*D == target
            return result, offset, list(decomp)
    raise AssertionError('target descent cap exhausted')


recovered, target_offset, target_indices = timed('target_online', target_online)
assert recovered == target_scalar


def rho_online():
    a, b = 1, 1
    P = D + target
    seen = {}
    for steps in range(1, 10001):
        key = str(P)
        if key in seen:
            old_a, old_b = seen[key]
            denominator = (b-old_b) % 109
            if denominator:
                answer = ((old_a-a)*pow(denominator, -1, 109)) % 109
                if answer*D == target:
                    return answer, steps
        else:
            seen[key] = (a, b)
        branch = hashlib.sha256(key.encode()).digest()[0] % 3
        if branch == 0:
            P += D
            a = (a+1) % 109
        elif branch == 1:
            P += target
            b = (b+1) % 109
        else:
            P += P
            a = (2*a) % 109
            b = (2*b) % 109
    raise AssertionError('rho cap exhausted')


rho_answer, rho_steps = timed('rho_online', rho_online)
assert rho_answer == recovered

out = {
    'status': 'verified_same_field_genus2_jacobian_pipeline_control',
    'source_sha256': SOURCE_SHA,
    'seed': SEED,
    'field': 'F_101', 'curve': 'y^2=x^3+x+46', 'curve_order': 109,
    'cover': 'v^2=u^6+u^2+46', 'cover_genus': 2,
    'cover_map': '(u,v)->(u^2,v)', 'jacobian_order': 11663,
    'jacobian_prime_subgroup_cofactor': 107,
    'factor_base_law': 'uniformly shuffled affine rational H points, one representative per distinct 107-projected class, first 32',
    'affine_points_nonzero_projection': affine_points_nonzero_projection,
    'distinct_projected_classes': distinct_projected,
    'actual_factor_base_size': B,
    'factor_base_sha256': hashlib.sha256(json.dumps([[v[2],v[3]] for v in base], separators=(',', ':')).encode()).hexdigest(),
    'pair_index_entries': sum(len(v) for v in pair_index.values()),
    'ordinary_query_law': 'MT19937 uniform a in [1,108], T=aD; search verified three-summand H-point decompositions',
    'ordinary_attempts': attempts, 'verified_relations': verified,
    'independent_rows': len(rows), 'final_rank': rank,
    'relation_witnesses': witnesses,
    'target_point_sha256': target_digest,
    'target_input_law': 'one fresh MT19937 scalar after precomputation; fixture scalar excluded from online algorithm',
    'target_offset_queries': target_offset+1,
    'target_decomposition_indices': target_indices,
    'target_scalar_recovered': recovered,
    'scalar_replay': bool(recovered*D == target),
    'rho_same_point_scalar_recovered': rho_answer,
    'rho_steps': rho_steps,
    'timings': timing,
    'timing_status': 'shared-host diagnostic, no isolated CPU speedup claim',
}
(HERE/'small_field_pipeline.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
print(json.dumps(out, indent=2, sort_keys=True))
