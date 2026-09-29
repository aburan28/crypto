#!/usr/bin/env python3
"""Exact toy S3 morphology audit. This does not benchmark full ECDLP."""
import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import subprocess
import time


def echelon(values):
    pivots = {}
    for value in values:
        while value:
            p = value.bit_length() - 1
            if p not in pivots:
                pivots[p] = value
                break
            value ^= pivots[p]
    return pivots


def contains(pivots, value):
    while value:
        p = value.bit_length() - 1
        if p not in pivots:
            return False
        value ^= pivots[p]
    return True


def span(basis):
    out = [0]
    for b in basis:
        out += [x ^ b for x in out]
    return out


def profile(F, basis, t):
    """Input is field, public subspace and target x only; no witness/scalar."""
    l, n = len(basis), F.deg
    products = [F.mul(a, b) for a in basis for b in basis]
    W = echelon(products)
    cross = [F.sqr(z) ^ F.mul(t, z) for z in products]
    linear = [F.mul(F.sqr(t), F.sqr(u)) for u in basis] * 2
    # [constant, linear columns, cross quadratic columns]. High-pivot
    # elimination eliminates all quadratic columns before affine columns.
    rows = []
    offset = 1 + 2*l
    for bit in range(n):
        row = int(bit == 0)  # a6=1
        for j, value in enumerate(linear):
            row |= ((value >> bit) & 1) << (j+1)
        for j, value in enumerate(cross):
            row |= ((value >> bit) & 1) << (offset+j)
        rows.append(row)
    full = echelon(rows)
    affine = [row for row in full.values() if row < (1 << offset)]
    quadratic_rank = len(echelon(cross))
    assert quadratic_rank == len(W) - int(contains(W, t))
    return {
        'product_span_dim': len(W),
        'square_closure_defect': len(echelon(basis+[F.sqr(x) for x in basis]))-l,
        'target_in_product_span': contains(W, t),
        'quadratic_rank': quadratic_rank,
        'affine_rank': len(echelon(row >> 1 for row in affine)),
        'affine_inconsistent': 1 in full.values(),
    }, affine


def linear_models(columns, rhs, n):
    """Solve columns*y=rhs; enumerate all models using RREF, count row XORs."""
    l = len(columns)
    rows = [sum(((v >> bit) & 1) << j for j, v in enumerate(columns)) |
            (((rhs >> bit) & 1) << l) for bit in range(n)]
    pivot_columns, pos, xors = [], 0, 0
    for col in range(l):
        pivot = next((i for i in range(pos, n) if (rows[i] >> col) & 1), None)
        if pivot is None:
            continue
        rows[pos], rows[pivot] = rows[pivot], rows[pos]
        for i in range(n):
            if i != pos and ((rows[i] >> col) & 1):
                rows[i] ^= rows[pos]
                xors += 1
        pivot_columns.append(col)
        pos += 1
    if any(row == (1 << l) for row in rows):
        return [], xors
    free = [j for j in range(l) if j not in pivot_columns]
    models = []
    for mask in range(1 << len(free)):
        model = sum(((mask >> i) & 1) << j for i, j in enumerate(free))
        for i, col in enumerate(pivot_columns):
            value = ((rows[i] >> l) & 1) ^ ((rows[i] & model).bit_count() & 1)
            model |= value << col
        models.append(model)
    return models, xors


def roots_fibers(F, basis, t):
    roots, xors = set(), 0
    xs = span(basis)
    for x in xs:
        columns = [F.mul(F.sqr(x ^ t), F.sqr(u)) ^ F.mul(F.mul(x, t), u)
                   for u in basis]
        rhs = F.sqr(F.mul(x, t)) ^ 1
        models, count = linear_models(columns, rhs, F.deg)
        xors += count
        roots.update((x, xs[mask]) for mask in models)
    return roots, xors


def roots_reference(F, basis, t):
    roots = set()
    xs = span(basis)
    for x in xs:
        for y in xs:
            z = F.mul(x, y)
            value = F.sqr(z) ^ F.mul(z, t) ^ F.mul(F.sqr(x ^ y), F.sqr(t)) ^ 1
            if value == 0:
                roots.add((x, y))
    return roots


def group_roots(curve, points, target):
    result = set()
    for x, left in points.items():
        for y, right in points.items():
            if any(curve.add(P, Q) == target for P in left for Q in right):
                result.add((x, y))
    return result


def bases_for(F, l):
    poly = [1 << i for i in range(l)]
    yield 'polynomial', 0, poly
    for seed in (17, 937):
        rng = random.Random(seed)
        scale = rng.randrange(1, 1 << F.deg)
        yield 'scaled_polynomial', seed, [F.mul(scale, b) for b in poly]
        basis = []
        while len(basis) < l:
            x = rng.randrange(1, 1 << F.deg)
            if not contains(echelon(basis), x):
                basis.append(x)
        yield 'random', seed, basis


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-repo', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    source = args.source_repo / 'scripts/ecc2k130_point_decomposition.py'
    spec = importlib.util.spec_from_file_location('pdp_source', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    started, cpu = time.perf_counter(), time.process_time()
    records, base_records = [], []
    for n, l in ((7,3),(9,4),(13,5),(17,6)):
        F = module.GF2m(n, module.find_irreducible(n), tables=True)
        curve = module.Koblitz(F)
        order = module.curve_order(n)
        r = max(module.prime_factors(order))
        assert module.is_prime(r)
        cofactor = order // r
        rng = random.Random(20260924+n)
        G = None
        while G is None:
            points = curve.points_over(rng.randrange(1, 1 << n))
            if points:
                G = curve.mul(points[0], cofactor)
        assert curve.on_curve(G) and curve.mul(G, r) is None
        natural = [curve.mul(G, rng.randrange(1, r)) for _ in range(32)]
        for family, seed, basis in bases_for(F, l):
            xs = span(basis)
            points = {x: curve.points_over(x) for x in xs}
            flat = [P for pp in points.values() for P in pp]
            planted = []
            attempts = 0
            prng = random.Random(3000+n+seed)
            while len(planted) < 4 and attempts < 512:
                attempts += 1
                T = curve.add(prng.choice(flat), prng.choice(flat))
                if T is not None and T[0] != 0 and curve.mul(T, r) is None:
                    planted.append(T)
            control_support_empty = False
            if len(planted) < 4:
                support = sorted({T for P in flat for Q in flat
                                  if (T := curve.add(P,Q)) is not None and T[0] != 0
                                  and curve.mul(T,r) is None})
                if support:
                    planted = [prng.choice(support) for _ in range(4)]
                else:
                    planted = []
                    control_support_empty = True
            base_id = f'n{n}-l{l}-{family}-{seed}'
            base_records.append({'base_id':base_id,'basis':basis,'curve_order':order,
                                 'prime_order':r,'cofactor':cofactor,'generator':G,
                                 'field_polynomial':F.irr,'factor_base_points':len(flat),
                                 'control_generation_attempts':attempts,'control_support_empty':control_support_empty})
            for kind, targets in (('natural',natural),('planted',planted)):
                for index, T in enumerate(targets):
                    assert T is not None and T[0] != 0 and curve.on_curve(T)
                    t0 = time.perf_counter()
                    features, affine = profile(F,basis,T[0])
                    profile_s = time.perf_counter()-t0
                    t0 = time.perf_counter()
                    reference = roots_reference(F,basis,T[0])
                    reference_s = time.perf_counter()-t0
                    t0 = time.perf_counter()
                    candidate, xors = roots_fibers(F,basis,T[0])
                    fiber_s = time.perf_counter()-t0
                    assert candidate == reference, (base_id,T,'fiber disagreement')
                    if features['affine_inconsistent']:
                        assert not reference
                    coordinates = {x: i for i,x in enumerate(xs)}
                    for x,y in reference:
                        assignment = 1 | (coordinates[x] << 1) | (coordinates[y] << (1+l))
                        assert all((row & assignment).bit_count()%2 == 0 for row in affine)
                    lifted = {(x,y) for x,y in reference if
                              any(curve.add(P,Q)==T for P in points[x] for Q in points[y])}
                    independent = group_roots(curve,points,T)
                    assert lifted == independent, (base_id,T,'group disagreement')
                    if kind == 'planted':
                        assert lifted
                    # Deterministic invertible change within the same V.
                    changed = basis[:]
                    changed[0] ^= changed[1]
                    changed.reverse()
                    changed_profile,_ = profile(F,changed,T[0])
                    assert changed_profile == features
                    changed_roots,_ = roots_fibers(F,changed,T[0])
                    assert changed_roots == reference
                    records.append({'base_id':base_id,'n':n,'l':l,'family':family,'seed':seed,
                                    'kind':kind,'index':index,'target':T,**features,
                                    'algebraic_roots':len(reference),'group_valid_x_pairs':len(lifted),
                                    'exhaustive_assignments':1<<(2*l),'fiber_systems':1<<l,
                                    'fiber_row_xors':xors,'profile_wall_s':profile_s,
                                    'reference_wall_s':reference_s,'fiber_wall_s':fiber_s,
                                    'all_gates_pass':True})
            print(base_id, 'complete', flush=True)
    groups = defaultdict(list)
    for row in records:
        if row['kind']=='natural':
            groups[(row['n'],row['family'])].append(row)
    summary=[]
    for (n,family),rows in groups.items():
        summary.append({'n':n,'family':family,'draws':len(rows),
                        'distinct_target_x':len({r['target'][0] for r in rows}),
                        'product_span_dims':sorted({r['product_span_dim'] for r in rows}),
                        'quadratic_ranks':sorted({r['quadratic_rank'] for r in rows}),
                        'affine_ranks':sorted({r['affine_rank'] for r in rows}),
                        'affine_refuted':sum(r['affine_inconsistent'] for r in rows),
                        'natural_group_hits':sum(bool(r['group_valid_x_pairs']) for r in rows),
                        'algebraic_only_hits':sum(bool(r['algebraic_roots']) and not r['group_valid_x_pairs'] for r in rows)})
    output = {'schema':'ecc2k130-pdp-morphology-v1','status':'TOY_EVIDENCE_PENDING_INDEPENDENT_REVIEW',
              'claim':'S3 structural diagnostic and exact-root verification; no full-DLP speed claim',
              'source_repo':str(args.source_repo.resolve()),'source_sha256':sha(source),
              'script_sha256':sha(Path(__file__)),'contract_sha256':sha(Path(__file__).with_name('pdp_profile_contract.md')),
              'source_head':subprocess.check_output(['git','-C',str(args.source_repo),'rev-parse','HEAD'],text=True).strip(),
              'host':platform.platform(),'records':records,'bases':base_records,'summary':summary,
              'checks':{'cases':len(records),'planted':sum(r['kind']=='planted' for r in records),
                        'root_set_disagreements':0,'group_disagreements':0,'rank_identity_disagreements':0,
                        'basis_change_disagreements':0},
              'timing':{'wall_s':time.perf_counter()-started,'cpu_s':time.process_time()-cpu},
              'full_dlp_total_operations':None,'rho_ratio':None,'solving_degree':None,
              'limitations':['m=2 only','complete roots workload, not first witness','natural draws repeated at small n',
                             'fixed dimension, unequal point cardinality','no learned policy or inferential success test',
                             'no n131 solver or representation transfer established']}
    (args.out/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'checks':output['checks'],'summary':summary,'timing':output['timing']},indent=2))

if __name__=='__main__':
    main()
