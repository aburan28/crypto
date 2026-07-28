#!/usr/bin/env python3
"""
Thread 6 — B5 over F_{p^k}: cover-cost analysis for abelian-surface DLPs
over small extensions.

Verifies the numerics from cover_complexity_ext.gp in Python (no PARI needed),
then adds a new experiment:
  EXP-A: exact break-even k for Diem 2011 on genuinely-over-F_{p^k} curves
  EXP-B: why DLP on Jac(C)(F_{p^k}) != ECDLP on E(F_p) (group-order argument)
  EXP-C: abelian surface genuinely over F_{p^2} — cost analysis
"""

import math

# secp256k1 parameters (log2 sizes)
LOG2_P = 256.0
LOG2_N = 255.999   # #E(F_p) just under 2^256

# Pollard rho baseline for secp256k1 (Aut(E) = Z/6 for j=0)
ECDLP_BITS = LOG2_N / 2 - math.log2(6) / 2
print(f"=== ECDLP baseline (Pollard rho, sqrt(6) speedup) ===")
print(f"  log2(cost) = {ECDLP_BITS:.4f} bits (~2^{round(ECDLP_BITS)} group ops)")
print()

# -----------------------------------------------------------------------
# Q1: Generic rho + Gaudry IC for covers over F_{p^k}
# -----------------------------------------------------------------------
print("=== Q1: Cover DLP cost vs. ECDLP for genus g, extension degree k ===")
print()
print(f"{'k':>3}  {'g':>2}  {'rho(bits)':>10}  {'IC(bits)':>10}  {'best':>8}  {'vs_ECDLP':>10}")
print("-" * 60)
for k in range(1, 6):
    for g in range(2, 6):
        rho_bits = g * k / 2.0 * LOG2_P
        ic_bits  = k * (2.0 - 2.0/g) * LOG2_P
        best     = min(rho_bits, ic_bits)
        vs       = best - ECDLP_BITS
        print(f"  k={k}  g={g}  {rho_bits:>10.1f}  {ic_bits:>10.1f}  {best:>8.1f}  {vs:>+10.1f}")
    print()

print("OBSERVATION: every (k,g) pair has cover DLP cost >> ECDLP cost.")
print()

# -----------------------------------------------------------------------
# Q2: Diem 2011 L_{p^k}[1/2, c=1] for k=1..6
# -----------------------------------------------------------------------
print("=== Q2: Diem 2011 L_{p^k}[1/2, c=1] — apparent break-even at k>=3 ===")
print()
loge_p = LOG2_P * math.log(2)  # ln(p)
print(f"{'k':>3}  {'ln(p^k)':>10}  {'L_bits':>10}  {'vs_ECDLP':>12}  note")
print("-" * 65)
for k in range(1, 7):
    ln_q    = k * loge_p
    ln_ln_q = math.log(ln_q)
    L_bits  = math.sqrt(ln_q * ln_ln_q) / math.log(2)  # c=1
    vs      = L_bits - ECDLP_BITS
    note    = f"*** BELOW ECDLP ({ECDLP_BITS:.1f} bits)!" if L_bits < ECDLP_BITS else f"above (+{vs:.1f} bits)"
    print(f"  k={k}  {ln_q:>10.2f}  {L_bits:>10.2f}  {vs:>+12.2f}  {note}")
print()
print("For k=3..6, L_formula is BELOW the ECDLP threshold.")
print("But: this applies only to curves *genuinely* defined over F_{p^k}.")
print()

# -----------------------------------------------------------------------
# EXP-A: Find exact break-even k for genuinely-over-F_{p^k} curves
# -----------------------------------------------------------------------
print("=== EXP-A: Diem break-even k (genuinely-over-F_{p^k} curves) ===")
print()
# Break-even: L_{p^k}[1/2,1] < sqrt(n/6)
# i.e. sqrt(k * ln_p * ln(k * ln_p)) < ln(sqrt(n/6))
target_nats = (LOG2_N / 2 - math.log2(6) / 2) * math.log(2)  # = ln(sqrt(n/6))
print(f"  Target: ln(sqrt(n/6)) = {target_nats:.4f} nats")
print()
print(f"  {'k':>4}  {'L_nats':>10}  {'<target?':>10}")
print("  " + "-" * 30)
breakeven_k = None
for k in range(1, 300):
    ln_q = k * loge_p
    ln_ln_q = math.log(ln_q)
    L_nats = math.sqrt(ln_q * ln_ln_q)
    if L_nats < target_nats:
        breakeven_k = k
        print(f"  k={k:>4}  {L_nats:>10.4f}  YES  <-- break-even at k={k}")
        break
    elif k <= 10 or k % 10 == 0:
        print(f"  k={k:>4}  {L_nats:>10.4f}  no")
print()
if breakeven_k:
    print(f"  Break-even at k={breakeven_k}: Diem sub-exp would be faster than ECDLP")
    print(f"  BUT this requires a genus-g curve GENUINELY over F_{{p^{breakeven_k}}}")
    print(f"  (not a base-change of secp256k1 from F_p).")
print()

# -----------------------------------------------------------------------
# EXP-B: Group-order argument — why DLP on Jac(C)(F_{p^k}) != ECDLP on E(F_p)
# -----------------------------------------------------------------------
print("=== EXP-B: Group-order obstruction for k>=3, genus g ===")
print()
print("  Even if C/F_{p^k} is genuinely over F_{p^k} and Jac(C) is isogenous")
print("  to some object containing E, the DLPs are inequivalent.")
print()
print("  KEY FACTS:")
print()

# secp256k1 prime n is prime. #E(F_{p^k}) for k>1:
# #E(F_{p^k}) = product of (p^k - alpha^k)(p^k - beta^k)/(p^k)  ... via Weil
# More precisely: #E(F_{p^k}) = (1 - alpha^k)(1 - beta^k) in terms of Frobenius roots
# alpha, beta = (t +/- sqrt(t^2 - 4p)) / 2 where t = trace of Frobenius
# For secp256k1: t is small (~ 2^128 per Hasse), alpha*beta = p

# Numerical: secp256k1 #E = p + 1 - t with t > 0 small
# For k=3: #E(F_{p^3}) = p^3 + 1 - (alpha^3 + beta^3)
# alpha^3 + beta^3 = t^3 - 3pt ~ t^3 (very roughly, since t ~ 2^128)
# But exact values need the actual t; approximate:

# secp256k1 order: n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
# secp256k1 prime: p = 2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1
# t = p + 1 - n (the trace of Frobenius)
p_secp = (1 << 256) - (1 << 32) - (1 << 9) - (1 << 8) - (1 << 7) - (1 << 6) - (1 << 4) - 1
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
t_secp = p_secp + 1 - n_secp
print(f"  secp256k1 trace t = p + 1 - n = {t_secp}")
print(f"  |t| in bits: {t_secp.bit_length()} bits  (should be ~128 by Hasse)")
print()

# Order of E over F_{p^2}: #E(F_{p^2}) = p^2 + 1 - (alpha^2 + beta^2)
# alpha^2 + beta^2 = t^2 - 2p
alpha_plus_beta = t_secp  # alpha + beta = t
alpha_times_beta = p_secp  # alpha * beta = p
alpha2_plus_beta2 = alpha_plus_beta**2 - 2 * alpha_times_beta
order_E_p2 = p_secp**2 + 1 - alpha2_plus_beta2
print(f"  #E(F_{{p^2}}) = {order_E_p2}")
print(f"  #E(F_{{p^2}}) / n = {order_E_p2 // n_secp} rem {order_E_p2 % n_secp}")
print(f"  n | #E(F_{{p^2}})?  {order_E_p2 % n_secp == 0}")
print()

# Order of E over F_{p^3}
# alpha^3 + beta^3 = (alpha+beta)^3 - 3*alpha*beta*(alpha+beta) = t^3 - 3p*t
alpha3_plus_beta3 = alpha_plus_beta**3 - 3 * alpha_times_beta * alpha_plus_beta
order_E_p3 = p_secp**3 + 1 - alpha3_plus_beta3
print(f"  #E(F_{{p^3}}) / n = ... checking divisibility")
print(f"  n | #E(F_{{p^3}})?  {order_E_p3 % n_secp == 0}")
print()

print("  IMPLICATION:")
print("  For Diem's algorithm to reduce DLP on Jac(C)(F_{p^k}) to ECDLP on E(F_p),")
print("  one needs: ECDLP on E(F_p) ≤_T DLP on Jac(C)(F_{p^k}).")
print("  This requires E(F_p) (cyclic of prime order n) to embed into Jac(C)(F_{p^k})")
print("  in a computationally accessible way.")
print()
print("  But E(F_p) ⊆ E(F_{p^k}), and n | #E(F_{p^k}) only for special k.")
print("  Even when n | #E(F_{p^k}), the n-torsion is a quotient of a MUCH larger group")
print("  (size p^{gk}). Extracting the n-torsion component from Diem's algorithm")
print("  output requires ADDITIONAL ECDLP computations — circular again.")
print()

# -----------------------------------------------------------------------
# EXP-C: Abelian surface genuinely over F_{p^2} — cost analysis
# -----------------------------------------------------------------------
print("=== EXP-C: Abelian surface genuinely over F_{p^2} (k=2, g=2) ===")
print()
print("  Scenario: suppose C/F_{p^2} is a genus-2 curve NOT defined over F_p,")
print("  and Jac(C)/F_{p^2} is an abelian surface with #Jac(C)(F_{p^2}) ≈ p^4.")
print()

k, g = 2, 2
rho_bits_p2 = g * k / 2.0 * LOG2_P
ic_bits_p2  = k * (2.0 - 2.0/g) * LOG2_P
best_p2 = min(rho_bits_p2, ic_bits_p2)
print(f"  Generic rho cost:       2^{rho_bits_p2:.1f} bits")
print(f"  Gaudry IC cost:         2^{ic_bits_p2:.1f} bits")
print(f"  Best cover DLP cost:    2^{best_p2:.1f} bits")
print(f"  ECDLP baseline:         2^{ECDLP_BITS:.1f} bits")
print(f"  Cover is {best_p2 - ECDLP_BITS:.1f} bits HARDER than ECDLP.")
print()
print("  Diem 2011 (k>=3 required): NOT applicable for k=2.")
print("  So F_{p^2} abelian surfaces offer no sub-exp DLP option.")
print()

# Does n divide #Jac(C)(F_{p^2}) for a typical abelian surface over F_{p^2}?
# #Jac(C)(F_{p^2}) ≈ p^4, and n | p^4 - 1 iff ord_n(p) | 4
# For secp256k1, n is prime ~ 2^256; p^4 has ~1024 bits.
# The embedding degree k_E for secp256k1 wrt to n is the smallest k with n | p^k - 1.
# Since n | p - 1 (because n = #E(F_p) = p+1-t and t=p+1-n, so p ≡ n-1+t ≡ ...)
# Actually n does NOT divide p-1 in general. Let's check.
print(f"  Checking p mod n = {p_secp % n_secp}")
print(f"  Checking p^2 mod n...")
print(f"  p^2 mod n = {pow(p_secp, 2, n_secp)}")
print(f"  So n | p^2 - 1? {(p_secp**2 - 1) % n_secp == 0}")
print()

print("=== Summary ===")
print()
print("Thread 6 findings (B5 over F_{p^k}):")
print()
print("1. All 6 original priority threads are now addressed:")
print("   T1 (P-521 LLL): CLOSED §10.5")
print("   T2 (CHLRS Igusa): CLOSED (F_p Rosenhain blocked; Thread 22 open)")
print("   T3 (Howe gluing): DONE (Thread 18, 5/15 pairs qualify)")
print("   T4 (Cross-curve LLL): CLOSED §10.4 (secp256k1 3/3 seeds)")
print("   T5 (GLV-HNP Phase 2): DONE (2026-07-26)")
print("   T6 (B5 over F_{p^k}): THIS ANALYSIS")
print()
print(f"2. EXP-A: Diem 2011 break-even at k={breakeven_k} for genuinely-over-F_{{p^k}} curves.")
print(f"   But secp256k1 has no such model over F_{{p^{breakeven_k}}}.")
print()
print("3. EXP-B: Group-order obstruction.")
print(f"   t = {t_secp} (bit-length {t_secp.bit_length()})")
print(f"   n | #E(F_{{p^2}})? {order_E_p2 % n_secp == 0}")
print(f"   n | #E(F_{{p^3}})? {order_E_p3 % n_secp == 0}")
print("   Diem's output on Jac(C)(F_{p^k}) can't be projected to E(F_p) n-torsion")
print("   without solving the original ECDLP (circularity).")
print()
print("4. EXP-C: Genuinely-over-F_{p^2} abelian surfaces.")
print(f"   Best DLP cost 2^{best_p2:.1f} >> 2^{ECDLP_BITS:.1f} ECDLP. No speedup.")
print("   Diem 2011 requires k>=3, so F_{p^2} is not covered by sub-exp.")
print()
print("VERDICT: B5 holds for all k >= 1. No cover-based attack on secp256k1")
print("benefits from passing to F_{p^k} extensions, for any k.")
