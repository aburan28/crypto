"""Frozen support-versus-uniform scalar character-cell diagnostic on P-192."""

import hashlib
import itertools
import json
import math
import random
import time
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, Integer, prod

import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
base_path = HERE/"type_i_line_guided_collector.json"
support_path = HERE/"type_i_line_three_sum_coverage.json"
base_receipt = json.loads(base_path.read_text())
support_receipt = json.loads(support_path.read_text())
assert support_receipt["distinct_three_sum_targets"] == math.comb(143,3) == 477191
p,r,L = dt.p,dt.r,dt.L
E = EllipticCurve(L,[L(dt.gm.a),L(dt.gm.b)])
G = E(L(dt.gm.gx),L(dt.gm.gy))
O = E(0)
base = [E(L(Integer(x)),L(Integer(y))) for x,y in base_receipt["base_points_fp"]]
assert len(base)==141 and r*G==O and G!=O
inv2 = Integer(pow(2,-1,int(r)))

def character(P):
    if P==O or P[0]==dt.c0:
        return "exceptional"
    X = dt.inf + dt.d0/(P[0]-dt.c0)
    rhs = prod(X-root**p for root in dt.quartic_roots)
    assert rhs**(p**3)==rhs
    if rhs==0:
        return "exceptional"
    return "+" if rhs**((p**3-1)//2)==1 else "-"

def cell(T):
    H=-inv2*T
    Q=2*H
    a,b=character(Q),character(H)
    return a+b if a!="exceptional" and b!="exceptional" else "exceptional"

seed=0x19220261016
rng=random.Random(seed)
n=512
indices=sorted(rng.sample(range(477191),n))
scalar_values=[Integer(rng.randrange(1,int(r))) for _ in range(n)]
start=time.monotonic_ns()
support_rows=[]
wanted=set(indices)
for index,(i,j,k) in enumerate(itertools.combinations_with_replacement(range(141),3)):
    if index in wanted:
        T=base[i]+base[j]+base[k]
        assert T!=O
        support_rows.append({"triple_index":index,"indices":[i,j,k],"cell":cell(T)})
        if len(support_rows)%64==0:
            print(f"LINE_CELL_SUPPORT_PROGRESS {len(support_rows)}/{n}",flush=True)
assert len(support_rows)==n
uniform_rows=[]
for index,a in enumerate(scalar_values):
    T=-2*a*G
    assert T!=O and r*T==O
    uniform_rows.append({"trial_index":index,"scalar":str(a),"cell":cell(T)})
    if (index+1)%64==0:
        print(f"LINE_CELL_UNIFORM_PROGRESS {index+1}/{n}",flush=True)
support_cells=Counter(row["cell"] for row in support_rows)
uniform_cells=Counter(row["cell"] for row in uniform_rows)
out={
    "status":"bounded_p192_line_support_vs_uniform_character_cell_sample",
    "candidate_id":None,
    "field_prime":str(p),"subgroup_order":str(r),
    "base_classes":141,"support_population":477191,
    "support_sample_without_replacement":n,
    "uniform_scalar_sample_with_replacement":n,
    "seed":seed,
    "support_cells":{c:support_cells[c] for c in ("++","+-","-+","--","exceptional")},
    "uniform_cells":{c:uniform_cells[c] for c in ("++","+-","-+","--","exceptional")},
    "support_rows":support_rows,
    "uniform_rows":uniform_rows,
    "elapsed_ns_unisolated":time.monotonic_ns()-start,
    "source_sha256":sha(__file__),
    "base_receipt_sha256":sha(base_path),
    "support_receipt_sha256":sha(support_path),
    "runtime_sha256":sha(HERE/"type_i_line_cell_sample_runtime_info.json"),
    "boundary":"Frozen exploratory 512+512 sample; support draw is uniform on all 477191 distinct unordered triples and ambient draw uses uniform nonzero scalars. Cell labels use exact F_(p^3) quadratic character. Ratios are descriptive with sampling uncertainty, not full-cell counts.",
}
(HERE/"type_i_line_cell_enrichment_sample.json").write_text(
    json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LINE_CELL_SAMPLE_PASS",json.dumps({
    "support_cells":out["support_cells"],"uniform_cells":out["uniform_cells"]},
    sort_keys=True),flush=True)
