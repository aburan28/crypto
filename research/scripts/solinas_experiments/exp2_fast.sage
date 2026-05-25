"""Twist attack - fast version using manual point selection (avoids slow Et.random_point)."""
import time, sys

def attack(label, p, b):
    print(f"=== {label} ===", flush=True)
    p = ZZ(p); Fp = GF(p)
    E = EllipticCurve(Fp, [-3, b])
    nE = E.order()
    t = p + 1 - nE
    nT = p + 1 + t
    nT_fact = factor(nT)
    nT_prime = max(q for q,_ in nT_fact)
    h_T = nT // nT_prime
    print(f"  nT = {nT}", flush=True)
    print(f"  nT factor = {nT_fact}", flush=True)
    print(f"  h_T = {h_T}, large prime tail bits = {nT_prime.nbits()}", flush=True)

    # Build twist: pick non-residue d, twist coefficients (a*d^2, b*d^3)
    d = Fp(1)
    while d.is_square():
        d = Fp.random_element()
        if d == 0: d = Fp(1)
    a2 = -3 * d**2
    b2 = b * d**3
    Et = EllipticCurve(Fp, [a2, b2])
    Et.set_order(nT)
    print(f"  twist built with a={a2.lift().nbits()}-bit, b={b2.lift().nbits()}-bit coefs", flush=True)

    # Find a random point on twist by sampling x until y^2 = x^3 + a2*x + b2 is a QR.
    t0 = time.time()
    tries = 0
    R = None
    while R is None:
        tries += 1
        x = Fp.random_element()
        rhs = x**3 + a2*x + b2
        if rhs.is_square():
            y = rhs.sqrt()
            R = Et(x, y)
    print(f"  random point found in {tries} tries ({(time.time()-t0)*1000:.1f}ms)", flush=True)

    # Project to order-h_T subgroup
    t1 = time.time()
    Q = (nT // h_T) * R
    if Q == Et(0):
        print(f"  unlucky: projected to identity, retry"); return None
    actual_ord = Q.order()
    print(f"  Q order = {actual_ord}  (want h_T = {h_T})  in {(time.time()-t1)*1000:.1f}ms", flush=True)
    if actual_ord != h_T:
        # Project further if order is a proper divisor of h_T
        print(f"  projecting from order {actual_ord} → not full h_T; will solve mod {actual_ord}")
        h_T_eff = actual_ord
    else:
        h_T_eff = h_T

    # Simulate oracle: server has secret k, returns k*Q
    k = ZZ.random_element(2, nT - 1)
    response = k * Q

    # Attacker recovers k mod h_T_eff via PH
    t2 = time.time()
    k_rec = discrete_log(response, Q, ord=h_T_eff, operation='+')
    dt = time.time() - t2
    leak = ceil(log(h_T_eff, 2))
    ok = (k % h_T_eff) == (k_rec % h_T_eff)
    print(f"  PH on order-{h_T_eff} subgroup: {dt*1000:.2f} ms", flush=True)
    print(f"  leakage = ceil(log2 {h_T_eff}) = {int(leak)} bits", flush=True)
    print(f"  recovered correctly: {ok}", flush=True)
    print(flush=True)
    return {'label':label, 'h_T':int(h_T_eff), 'leak':int(leak), 't_ms':float(dt*1000), 'ok':bool(ok)}

print("=" * 70, flush=True)
print("Exp 2 (fast): Twist attack with manual point sampling", flush=True)
print("=" * 70, flush=True)
print(flush=True)

# Just P-256/96 since P-224/56 already confirmed
r = attack("P-256/96", 0xfff002000000000fffffffff, 7)
print("=" * 70, flush=True)
print("Summary (combining with P-224/56 from earlier run):", flush=True)
print(f"  P-224/56  h_T=33  PH=0.79ms  leaked=6 bits  ok=True", flush=True)
if r:
    print(f"  P-256/96  h_T={r['h_T']}  PH={r['t_ms']:.2f}ms  leaked={r['leak']} bits  ok={r['ok']}", flush=True)
