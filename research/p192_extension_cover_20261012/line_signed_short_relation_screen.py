"""Exact signed mixed-length additive-relation screen on 141 P-192 classes."""

import hashlib
import json
import math
import time
from pathlib import Path

HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source_path=HERE/"type_i_line_guided_collector.json"
prior_path=HERE/"type_i_line_three_sum_coverage.json"
source=json.loads(source_path.read_text())
prior=json.loads(prior_path.read_text())
p=int(source["field_prime"])
r=int(source["subgroup_order"])
a=p-3
b=int("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1",16)
assert p==int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF",16)
base=[tuple(map(int,xy)) for xy in source["base_points_fp"]]
assert len(base)==len(set(base))==141

def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x,y=P
    u,v=Q
    if x==u and (y+v)%p==0: return None
    if P==Q:
        slope=(3*x*x+a)*pow(2*y,-1,p)%p
    else:
        slope=(v-y)*pow((u-x)%p,-1,p)%p
    z=(slope*slope-x-u)%p
    return z,(slope*(x-z)-y)%p

def neg(P):
    return None if P is None else (P[0],(-P[1])%p)

def mul(n,P):
    R=None
    while n:
        if n&1: R=add(R,P)
        P=add(P,P)
        n>>=1
    return R

assert all((y*y-x*x*x-a*x-b)%p==0 and mul(r,(x,y)) is None
           for x,y in base)
assert r%2==1
start=time.monotonic_ns()
support={1:{P:(i,) for i,P in enumerate(base)},2:{},3:{}}
for i in range(141):
    for j in range(i,141):
        P=add(base[i],base[j])
        assert P not in support[2]
        support[2][P]=(i,j)
        for k in range(j,141):
            T=add(P,base[k])
            assert T not in support[3]
            support[3][T]=(i,j,k)
    if (i+1)%32==0:
        print(f"LINE_SIGNED_SHORT_PROGRESS {i+1}/141",flush=True)
assert [len(support[n]) for n in (1,2,3)]==[141,10011,477191]
assert len(support[3])==prior["distinct_three_sum_targets"]
assert None not in support[3]
support_hash=hashlib.sha256()
for x,y in sorted(support[3]):
    support_hash.update(x.to_bytes(24,"big")+y.to_bytes(24,"big"))
assert support_hash.hexdigest()==prior["support_sha256_sorted_xy_be24"]

witnesses=[]
counts={}
for n in (1,2,3):
    for m in range(n,4):
        if n<m:
            key=f"{n}_vs_{m}_same_sign"
            matches=[]
            for P,left in support[n].items():
                if P in support[m]:
                    matches.append((left,support[m][P],P))
            counts[key]=len(matches)
            witnesses += [{"kind":key,"left":list(L),"right":list(R)}
                          for L,R,_ in matches]
        key=f"{n}_vs_{m}_opposite_sign"
        matches=[]
        for P,left in support[n].items():
            other=support[m].get(neg(P))
            if other is not None:
                if n==m and left>other:
                    continue
                assert neg(P) in support[m]
                matches.append((left,other,P))
        counts[key]=len(matches)
        witnesses += [{"kind":key,"left":list(L),"right":list(R)}
                      for L,R,_ in matches]
    counts[f"{n}_term_zero"]=int(None in support[n])

def rank_mod(rows,modulus):
    pivots={}
    for row in rows:
        v=[x%modulus for x in row]
        for j in sorted(pivots):
            if v[j]:
                q=v[j]
                for k in range(j,len(v)):
                    v[k]=(v[k]-q*pivots[j][k])%modulus
        lead=next((j for j,x in enumerate(v) if x),None)
        if lead is not None:
            inv=pow(v[lead],-1,modulus)
            pivots[lead]=[(x*inv)%modulus for x in v]
    return len(pivots)

rows=[]
for w in witnesses:
    row=[0]*141
    for i in w["left"]: row[i]+=1
    sign=-1 if "opposite_sign" in w["kind"] else 1
    for i in w["right"]: row[i]-=sign
    assert sum(row)==len(w["left"])-sign*len(w["right"])
    rows.append(row)
assert all(count==0 for count in counts.values())
assert not witnesses
out={
 "status":"verified_p192_line_base_signed_short_relation_screen",
 "candidate_id":None,
 "field_prime":str(p),"subgroup_order":str(r),"base_classes":141,
 "support_counts":{"1":141,"2":10011,"3":477191},
 "support_sha256_sorted_xy_be24":support_hash.hexdigest(),
 "overlap_counts":counts,
 "witnesses":witnesses,
 "collision_relation_rank":rank_mod(rows,r),
 "max_covered_terms":6,
 "scope":"same- and opposite-sign overlaps of all 1-, 2-, and 3-term unordered supports, including cross-length equalities and identity targets; no longer support tested",
 "elapsed_ns_unisolated":time.monotonic_ns()-start,
 "source_sha256":sha(__file__),
 "input_sha256":sha(source_path),
 "prior_support_sha256":sha(prior_path),
 "boundary":"Exact 141-class group-equation screen. Zero short-support overlaps do not bound longer relations or geometry-conditioned known-scalar target yield.",
}
(HERE/"type_i_line_signed_short_relation_screen.json").write_text(
    json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LINE_SIGNED_SHORT_PASS",json.dumps({
    "overlap_counts":counts,"rank":out["collision_relation_rank"]},sort_keys=True),flush=True)
