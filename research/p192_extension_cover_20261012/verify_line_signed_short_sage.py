"""Independent Sage replay of all signed 1/2/3-support intersections."""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer

HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
input_path=HERE/"type_i_line_guided_collector.json"
screen_path=HERE/"type_i_line_signed_short_relation_screen.json"
input=json.loads(input_path.read_text())
screen=json.loads(screen_path.read_text())
p=Integer(input["field_prime"])
r=Integer(input["subgroup_order"])
assert p==Integer("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF",16)
F=GF(p)
E=EllipticCurve(F,[F(p-3),F(Integer("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1",16))])
O=E(0)
base=[E(F(x),F(y)) for x,y in input["base_points_fp"]]
assert len(base)==len(set(base))==141
assert all(Q!=O and r*Q==O for Q in base)

def key(Q):
    return None if Q==O else (int(Q[0]),int(Q[1]))

def neg(P):
    return None if P is None else (P[0],int((-F(P[1]))))

supports={1:{key(Q) for Q in base},2:set(),3:set()}
for i in range(141):
    for j in range(i,141):
        S=base[i]+base[j]
        supports[2].add(key(S))
        for k in range(j,141):
            supports[3].add(key(S+base[k]))
    if (i+1)%32==0:
        print(f"LINE_SIGNED_SAGE_PROGRESS {i+1}/141",flush=True)
assert [len(supports[n]) for n in (1,2,3)]==[141,10011,477191]
assert None not in supports[1] and None not in supports[2] and None not in supports[3]
h=hashlib.sha256()
for x,y in sorted(supports[3]):
    h.update(x.to_bytes(24,"big")+y.to_bytes(24,"big"))
assert h.hexdigest()==screen["support_sha256_sorted_xy_be24"]
counts={}
for n in (1,2,3):
    for m in range(n,4):
        if n<m:
            key_name=f"{n}_vs_{m}_same_sign"
            counts[key_name]=len(supports[n]&supports[m])
        key_name=f"{n}_vs_{m}_opposite_sign"
        counts[key_name]=len({neg(P) for P in supports[n]}&supports[m])
    counts[f"{n}_term_zero"]=int(None in supports[n])
assert counts==screen["overlap_counts"]
assert all(v==0 for v in counts.values())
out={
 "status":"LINE_SIGNED_SHORT_SAGE_INDEPENDENT_PASS",
 "base_classes":141,
 "support_counts":{"1":141,"2":10011,"3":477191},
 "overlap_counts":counts,
 "support_sha256_sorted_xy_be24":h.hexdigest(),
 "source_sha256":sha(__file__),
 "input_sha256":sha(input_path),
 "screen_sha256":sha(screen_path),
 "runtime_sha256":sha(HERE/"type_i_line_signed_sage_runtime_info.json"),
}
(HERE/"type_i_line_signed_short_sage_verification.json").write_text(
    json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LINE_SIGNED_SHORT_SAGE_INDEPENDENT_PASS",json.dumps({
 "support_counts":out["support_counts"],"overlap_count":sum(counts.values())},
 sort_keys=True),flush=True)
