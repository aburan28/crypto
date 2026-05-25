"""
Experiment 1: Pohlig-Hellman DLP in F_p* on the P-224-family scale-down primes.
Plus an estimate of cost on real P-224.

Hypothesis: PH should be fast precisely because p-1 is super-smooth
(by the cyclotomic identity p - 1 = 2^k * (2^(n-k) - 1) for the P-224 form).
"""
import time

P224_FAMILY = [
    (56,  0xffffffff000001),
    (64,  0xffffffffff000001),
    (70,  0x3fffffff0000000001),
    (80,  0xffffffff000000000001),
    (84,  0xffff00000000000000001),
    (96,  0xffffffffffffffff00000001),
    (126, 0x3ffffffffffffffffffffffffffc0001),
]
P224_REAL = ("P-224 (real)", 224, 2^224 - 2^96 + 1)

def find_generator(p):
    """Find a generator of F_p*."""
    F = GF(p)
    phi = p - 1
    fact = factor(phi)
    for g in range(2, 200):
        gF = F(g)
        ok = True
        for q, _ in fact:
            if gF^(phi // q) == 1:
                ok = False; break
        if ok:
            return gF
    raise RuntimeError("no generator < 200")

def run_pohlig_hellman(p, verbose=True):
    F = GF(p)
    phi = p - 1
    fact = factor(phi)
    largest = max(q for q, _ in fact)
    t0 = time.time()
    g = find_generator(p)
    t_gen = time.time() - t0
    # Pick a random "secret" x and compute h = g^x
    x_secret = ZZ.random_element(2, phi - 1)
    h = g^x_secret
    # Solve DLP via Pohlig-Hellman (Sage's discrete_log uses PH+BSGS automatically)
    t1 = time.time()
    x_recovered = discrete_log(h, g, ord=phi)
    t_dlp = time.time() - t1
    ok = (x_recovered % phi) == (x_secret % phi)
    if verbose:
        print(f"  p ≈ 2^{p.nbits()-1}    largest_prime(p-1) = {largest}  ({ZZ(largest).nbits()} bits)")
        print(f"    DLP recovery time: {t_dlp*1000:.1f} ms     correct: {ok}")
        print(f"    factorization of p-1:")
        # show only first 6 factors + last to keep output tight
        flist = list(fact)
        for q, e in flist[:6]:
            print(f"       {q}^{e}" if e > 1 else f"       {q}")
        if len(flist) > 7:
            print(f"       ... ({len(flist) - 7} omitted) ...")
            q, e = flist[-1]
            print(f"       {q}^{e}" if e > 1 else f"       {q}")
        elif len(flist) > 6:
            q, e = flist[-1]
            print(f"       {q}^{e}" if e > 1 else f"       {q}")
    return {'n': p.nbits()-1, 'largest_pf_p_minus_1': int(largest),
            'largest_bits': int(ZZ(largest).nbits()),
            't_dlp_ms': t_dlp*1000, 'ok': ok}

print("=" * 70)
print("Experiment 1: Pohlig-Hellman in F_p* (P-224 family)")
print("=" * 70)
print()

results = []
for n, p in P224_FAMILY:
    print(f"--- P-224/{n} ---")
    r = run_pohlig_hellman(ZZ(p))
    results.append(r)
    print()

# Estimate for real P-224 (don't run; would take days)
print(f"--- {P224_REAL[0]} (estimate, not run) ---")
p_real = P224_REAL[2]
phi = p_real - 1
fact_real = factor(phi)
largest_real = max(q for q, _ in fact_real)
print(f"  largest_prime(p-1) = {largest_real}  ({ZZ(largest_real).nbits()} bits)")
print(f"  PH cost ≈ √{largest_real} ≈ 2^{ZZ(largest_real).nbits()/2:.1f} group ops")
print(f"  At 10^7 modexps/sec → ≈ {float(sqrt(largest_real) / 10^7 / 60):.1f} minutes for the bottleneck prime")
print(f"  (random 224-bit prime would have largest factor ~ 2^224 → BSGS infeasible)")

print()
print("=" * 70)
print("Summary (PH on F_p* for P-224 family):")
print("=" * 70)
print(f"  {'n':>4s}  {'largest pf(p-1)':>20s}  {'bits':>5s}  {'PH time':>10s}  ok")
for r in results:
    print(f"  {r['n']:4d}  {r['largest_pf_p_minus_1']:20d}  {r['largest_bits']:5d}  {r['t_dlp_ms']:8.1f}ms  {r['ok']}")
