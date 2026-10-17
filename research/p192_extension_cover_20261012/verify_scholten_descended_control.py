"""Independent arithmetic replay of the frozen descended Scholten control."""
import hashlib
import json
import random
from pathlib import Path

from sage.all import EllipticCurve, GF, HyperellipticCurve, PolynomialRing, matrix, vector

HERE = Path(__file__).resolve().parent
receipt = json.loads((HERE/'scholten_descended_control.json').read_text())
source = HERE/'scholten_descended_control.py'
assert receipt['source_sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
assert receipt['status'] == 'verified_one_target_descended_scholten_control'
p,r = receipt['field_p'],receipt['elliptic_order_r']
F = GF(p)
R = PolynomialRing(F,'x')
x = R.gen()
a,b = receipt['elliptic_coefficients']
E = EllipticCurve(F,[a,b])
assert E.cardinality() == r
G = E(*receipt['elliptic_generator'])
assert r*G == E(0) and G != E(0)
K3 = F.extension(x**3+F(a)*x+F(b),'zeta')
zeta = K3.gen()
P = PolynomialRing(K3,'T')
T = P.gen()
def decode(rows):
    return P([sum((K3(F(int(v)))*zeta**i for i,v in enumerate(row)),K3(0))
              for row in rows])
Hpoly = decode(receipt['cover_polynomial_coefficients_Fp3'])
assert Hpoly.is_squarefree()
J = HyperellipticCurve(Hpoly).jacobian()
N = J([decode(receipt['norm_u_coefficients_Fp3']),
       decode(receipt['norm_v_coefficients_Fp3'])])
assert N != J(0) and r*N == J(0)

x0 = K3(receipt['anchor_x'])
w0 = sum((K3(F(int(v)))*zeta**i
          for i,v in enumerate(receipt['anchor_w_coefficients'])),K3(0))
assert w0**2 == Hpoly(x0)
h = receipt['projection_cofactor']
trace = p+1-r
s0,s1 = 2,trace
for _ in range(2,7): s0,s1 = s1,trace*s1-p*s0
assert receipt['extension_order'] == p**6+1-s1 == h*r
base = []
for xx,coeffs in receipt['base_selection']:
    xx = K3(xx)
    ww = sum((K3(F(int(v)))*zeta**i for i,v in enumerate(coeffs)),K3(0))
    assert xx != x0 and ww**2 == Hpoly(xx)
    U = (T-x0)*(T-xx)
    V = P(w0+(ww-w0)*(T-x0)/(xx-x0))
    assert (V**2-Hpoly)%U == 0
    Q = h*J([U,V])
    assert Q != J(0) and r*Q == J(0) and Q not in base
    base.append(Q)
assert len(base) == receipt['base_classes']

rng = random.Random(receipt['ordinary_query_seed'])
assert receipt['ordinary_query_scalars'] == [
    rng.randrange(1,r) for _ in range(receipt['ordinary_queries'])]
assert rng.randrange(1,r) == receipt['target_scalar_fixture']

FR = GF(r)
rows, rhs = [], []
for witness in receipt['relation_witnesses']:
    indices = witness['indices']
    scalar = witness['scalar']
    assert sum((base[i] for i in indices),J(0)) == scalar*N
    row = [0]*len(base)
    for i in indices: row[i] += 1
    rows.append(row); rhs.append(scalar)
M = matrix(FR,rows)
assert M.rank() == len(base) == receipt['final_rank']
logs = vector(FR,receipt['base_logs'])
assert M*logs == vector(FR,rhs)
for i,Q in enumerate(base): assert int(logs[i])*N == Q

target = receipt['target_scalar_fixture']*N
assert hashlib.sha256(str(target).encode()).hexdigest() == receipt['target_point_sha256']
assert receipt['recovered_scalar']*N == target
assert receipt['rho_recovered_scalar']*N == target
assert receipt['recovered_scalar'] == receipt['rho_recovered_scalar']

out = {'status':'independently_verified_descended_control_receipt',
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'run_source_sha256':receipt['source_sha256'],
       'base_classes':len(base),'relation_witnesses':len(rows),
       'rank':int(M.rank()),'target_scalar':receipt['recovered_scalar'],
       'same_point_rho_scalar':receipt['rho_recovered_scalar']}
(HERE/'scholten_descended_control_verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
