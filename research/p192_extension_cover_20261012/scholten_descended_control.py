"""Complete one-target IC control on a genuine degree-three Scholten descent.

The F_23 instance uses the same algebraic construction as the P-192 branch.
Every base and target class is checked in the descended Jacobian over F_(23^3).
"""
import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, HyperellipticCurve, Integer, PolynomialRing, matrix, prod, vector

HERE = Path(__file__).resolve().parent
source_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
p, a, b, r = Integer(23), Integer(1), Integer(4), Integer(29)
F = GF(p)
E = EllipticCurve(F, [F(a), F(b)])
assert E.cardinality() == r and r.is_prime()
G = next(P for P in E.points() if P != E(0))
assert r*G == E(0)
R = PolynomialRing(F, 'x')
x = R.gen()
f = x**3+F(a)*x+F(b)
assert f.is_irreducible() and p % 3 == 2
K3 = F.extension(f, 'zeta')
zeta = K3.gen()
s = -(3*zeta**2+K3(a))/(3*zeta)
A, B = s**3, 3*s**3+3*zeta*s**2
L = GF(p**6, 'k')
root = f.change_ring(L).roots(multiplicities=False)[0]
embed = lambda c: L(c.polynomial()(root))
z, sl, AL, BL = [embed(c) for c in (zeta,s,A,B)]
eta = L(-3).sqrt()
omega = (-1+eta)/2
assert omega**2+omega+1 == 0 and omega**p**3 == omega**2
P6 = PolynomialRing(L, 'T')
T6 = P6.gen()
numerator = (AL*(eta+T6)**6 + BL*omega**2*(eta+T6)**4*(eta-T6)**2
             + BL*omega*(eta+T6)**2*(eta-T6)**4 + AL*(eta-T6)**6)
H6poly = numerator/(2*eta)**6
sigma = lambda c: c**(p**3)
assert all(sigma(c) == c for c in H6poly.list())

basis = [L(1),z,z**2]
def prime_coords(c):
    cs = c.polynomial().list()
    return vector(F, cs+[F(0)]*(6-len(cs)))
embedding_matrix = matrix(F,[prime_coords(c) for c in basis]).transpose()
assert embedding_matrix.rank() == 3
def descend(c):
    assert sigma(c) == c
    coords = embedding_matrix.solve_right(prime_coords(c))
    out = sum((K3(coords[i])*zeta**i for i in range(3)),K3(0))
    assert embed(out) == c
    return out

P3 = PolynomialRing(K3, 'T')
T3 = P3.gen()
H3poly = P3([descend(c) for c in H6poly.list()])
assert H3poly.is_squarefree()
J6 = HyperellipticCurve(H6poly).jacobian()
J3 = HyperellipticCurve(H3poly).jacobian()
tG = (L(G[0])-z)/sl-1
U = P6(omega*(eta+T6)**2-tG*(eta-T6)**2).monic()
V = P6(L(G[1])*(eta-T6)**3/(2*eta)**3) % U
assert (V**2-H6poly)%U == 0
D = J6([U,V])
Ds = J6([P6([sigma(c) for c in U.list()]),
         P6([sigma(c) for c in V.list()])])
N6 = D+Ds
assert N6 != J6(0) and r*N6 == J6(0)
N = J3([P3([descend(c) for c in N6[0].list()]),
        P3([descend(c) for c in N6[1].list()])])
assert N != J3(0) and r*N == J3(0)

# Exact Frobenius recurrence gives #E(F_(p^6)). Each projected Jacobian
# point is independently checked as nonzero r-torsion.
trace = p+1-r
s0,s1 = Integer(2),trace
for _ in range(2,7):
    s0,s1 = s1,trace*s1-p*s0
extension_order = p**6+1-s1
assert extension_order%r == 0
cofactor = extension_order//r
assert cofactor%r != 0

setup_done = time.monotonic_ns()
start = setup_done
eligible = []
for xx in F:
    rhs = H3poly(K3(xx))
    if rhs.is_square():
        ww = rhs.sqrt()
        for sign in (1,-1):
            eligible.append((K3(xx),sign*ww))
anchor = eligible[0]
base, selection = [], []
for xx,ww in eligible[1:]:
    x0,w0 = anchor
    if xx == x0:
        continue
    U3 = (T3-x0)*(T3-xx)
    V3 = P3(w0+(ww-w0)*(T3-x0)/(xx-x0))
    assert (V3**2-H3poly)%U3 == 0
    Q = cofactor*J3([U3,V3])
    assert r*Q == J3(0)
    if Q == J3(0) or Q in base:
        continue
    base.append(Q)
    selection.append((int(xx),[str(c) for c in ww.polynomial().list()]))
    if len(base) == 8:
        break
assert len(base) == 8
base_ns = time.monotonic_ns()-start

start = time.monotonic_ns()
pairs = {}
for i in range(len(base)):
    for j in range(i,len(base)):
        pairs.setdefault(str(base[i]+base[j]),[]).append((i,j))
pair_ns = time.monotonic_ns()-start

def find_three(T):
    for k in range(len(base)):
        for i,j in pairs.get(str(T-base[k]),[]):
            if base[i]+base[j]+base[k] == T:
                return i,j,k
    return None

FR = GF(r)
rng = random.Random(0x2320261015)
rows,rhs,witnesses = [],[],[]
attempts,verified,rank = 0,0,0
queried_scalars = []
start = time.monotonic_ns()
while attempts < 500 and rank < len(base):
    scalar = rng.randrange(1,int(r))
    queried_scalars.append(scalar)
    attempts += 1
    indices = find_three(scalar*N)
    if indices is None:
        continue
    assert sum((base[i] for i in indices),J3(0)) == scalar*N
    verified += 1
    row = [0]*len(base)
    for i in indices:
        row[i] += 1
    new_rank = matrix(FR,rows+[row]).rank()
    if new_rank > rank:
        rows.append(row)
        rhs.append(scalar)
        witnesses.append({'scalar':scalar,'indices':list(indices)})
        rank = new_rank
relation_ns = time.monotonic_ns()-start
assert rank == len(base)

start = time.monotonic_ns()
M = matrix(FR,rows)
logs = [int(v) for v in M.solve_right(vector(FR,rhs))]
assert M*vector(FR,logs) == vector(FR,rhs)
for log,Q in zip(logs,base):
    assert log*N == Q
la_and_replay_ns = time.monotonic_ns()-start

target_scalar = rng.randrange(1,int(r))
target = target_scalar*N
target_digest = hashlib.sha256(str(target).encode()).hexdigest()
start = time.monotonic_ns()
recovered = None
for offset in range(int(r)):
    indices = find_three(target+offset*N)
    if indices is not None:
        recovered = (sum(logs[i] for i in indices)-offset)%int(r)
        assert recovered*N == target
        break
target_online_ns = time.monotonic_ns()-start
assert recovered == target_scalar

start = time.monotonic_ns()
alpha,beta = 1,1
Q = N+target
seen = {}
for steps in range(1,10001):
    key = str(Q)
    if key in seen:
        old_alpha,old_beta = seen[key]
        denominator = (beta-old_beta)%int(r)
        if denominator:
            answer = ((old_alpha-alpha)*pow(denominator,-1,int(r)))%int(r)
            if answer*N == target:
                break
    else:
        seen[key] = (alpha,beta)
    branch = hashlib.sha256(key.encode()).digest()[0]%3
    if branch == 0:
        Q += N; alpha = (alpha+1)%int(r)
    elif branch == 1:
        Q += target; beta = (beta+1)%int(r)
    else:
        Q += Q; alpha = (2*alpha)%int(r); beta = (2*beta)%int(r)
else:
    raise AssertionError('rho cap exhausted')
rho_online_ns = time.monotonic_ns()-start
assert answer == target_scalar

out = {
    'status':'verified_one_target_descended_scholten_control',
    'source_sha256':source_sha,
    'field_p':int(p),'elliptic_coefficients':[int(a),int(b)],
    'elliptic_order_r':int(r),'elliptic_generator':[int(G[0]),int(G[1])],
    'cover_field_degree':3,'cover_genus':2,
    'cover_polynomial_coefficients_Fp3':[[str(v) for v in c.polynomial().list()] for c in H3poly.list()],
    'norm_u_coefficients_Fp3':[[str(v) for v in c.polynomial().list()] for c in N[0].list()],
    'norm_v_coefficients_Fp3':[[str(v) for v in c.polynomial().list()] for c in N[1].list()],
    'norm_class_nonzero':True,'r_times_norm_zero':True,
    'extension_order':int(extension_order),'projection_cofactor':int(cofactor),
    'eligible_affine_sign_choices':len(eligible),
    'anchor_x':int(anchor[0]),
    'anchor_w_coefficients':[str(c) for c in anchor[1].polynomial().list()],
    'base_selection':selection,'base_classes':len(base),
    'pair_index_entries':sum(len(v) for v in pairs.values()),
    'ordinary_query_seed':0x2320261015,
    'ordinary_query_scalars':queried_scalars,
    'ordinary_queries':attempts,'verified_relations':verified,
    'independent_rows':len(rows),'final_rank':int(rank),'relation_witnesses':witnesses,
    'base_logs':logs,'target_point_sha256':target_digest,
    'target_scalar_fixture':target_scalar,'recovered_scalar':recovered,
    'target_offset_queries':offset+1,
    'rho_recovered_scalar':answer,'rho_walk_steps':steps,
    'base_ns':base_ns,'pair_index_ns':pair_ns,
    'relation_collection_ns':relation_ns,'final_la_and_replay_ns':la_and_replay_ns,
    'target_online_ns':target_online_ns,'rho_online_ns':rho_online_ns,
    'timing_boundary':'same-point target-dependent intervals only; shared-host diagnostics, not isolated CPU speedup',
}
(HERE/'scholten_descended_control.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
