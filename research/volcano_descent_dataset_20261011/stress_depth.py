import sys, time, random, json
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from volcano_descent import *
rng = random.Random(99)
res = []
for (ell, D, k, bits) in [(2, -7, 20, 56), (3, -8, 12, 64)]:
    inst = find_instance(ell, D, 1, k, 0, bits, rng, max_tries=100000)
    t0 = time.time()
    V = explore_from_crater(inst, max_vertices=5_000_000)
    dt = time.time() - t0
    ok = level_counts(V) == expected_level_counts(inst)
    jc = [j for j, l in V.level.items() if l == 0][0]
    t1 = time.time(); r = descend_to_floor(jc, inst.p, inst.t, ell, rng); d1 = time.time() - t1
    row = dict(ell=ell, D=D, depth=k, p_bits=inst.p.bit_length(), vertices=len(V.level), level_counts_match_class_numbers=ok,
               bfs_seconds=round(dt, 1), blind_descent_from_crater_seconds=round(d1, 3), steps=len(r['path']) - 1,
               floor_level=V.level[r['floor_j']])
    print(json.dumps(row), flush=True)
    res.append(row)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'stress_depth.json'), 'w'), indent=1)
