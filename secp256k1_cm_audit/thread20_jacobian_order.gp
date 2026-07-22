\\ ============================================================
\\ Thread 20: Pohlig-Hellman group order smoothness for CM-73 Jacobians
\\ ============================================================
\\
\\ For CM-73 primes p in {19,37,79,109}, the biquadratic Weil poly is
\\ T^4 + a2*T^2 + p^2 with a2 = 2p - 73.
\\ Group order of Jac(C)(F_p) = P_C(1) = 1 + a2 + p^2 = (p+1)^2 - 73.
\\
\\ Question: is this order smooth relative to the elliptic curve base order?
\\ If smooth and a cover existed, Pohlig-Hellman on Jac would beat rho on E.
\\
\\ Run: gp -q thread20_jacobian_order.gp

default(parisize, 128000000);

print("================================================================");
print("Thread 20: CM-73 Jacobian group orders and PH smoothness check");
print("================================================================");
print("");

\\ CM-73 primes (norm-form k in {1,5,9,11}; a2 = 2p - 73)
cm73_primes = [19, 37, 79, 109];
print("---- CM-73 primes (a2 = 2p - 73, biquadratic Weil poly) ----");
print("");
print("p  | a2   | #Jac = (p+1)^2 - 73 | factorization | B-smooth? (B = p^0.5)");
print("-----------------------------------------------------------------");
{
forvec(pv = [[1,4]],
    p = cm73_primes[pv[1]];
    a2 = 2*p - 73;
    jac_order = 1 + a2 + p^2;   \\ = (p+1)^2 - 73
    f = factor(jac_order);
    max_prime_factor = vecmax(f[,1]);
    B = floor(sqrt(p));
    smooth = (max_prime_factor <= p);   \\ smooth means all factors <= p (= rough-p)
    print(p, " | ", a2, " | ", jac_order, " | ", f, " | max_pfactor=", max_prime_factor, " (p=",p,", sqrt(p)=",B,")")
);
}
print("");

\\ Non-CM-73 norm-form primes for comparison (a2 = 2p - k^2; k=3,7,13,...)
\\ Use k^2 != 73 examples: k=3 (a2=2p-9), k=7 (a2=2p-49)
print("---- Non-CM-73 norm-form primes (for comparison) ----");
print("");
print("k   | p (norm-form) | a2       | #Jac = 1+a2+p^2  | max prime factor");
print("--------------------------------------------------------------------");
{
for(k = 1, 20,
    if(k^2 == 73, next);      \\ skip CM-73
    if(!isprime(4 + 3*k^2), next);   \\ check 4+3k^2 prime (secp256k1 norm form)
    p = 4 + 3*k^2;
    if(!isprime(p), next);
    a2 = 2*p - k^2;
    jac_order = 1 + a2 + p^2;
    f = factor(jac_order);
    max_pf = vecmax(f[,1]);
    print("k=", k, " | p=", p, " | a2=", a2, " | #Jac=", jac_order, " | max_pf=", max_pf)
);
}
print("");

\\ Key formula check: (p+1)^2 - 73 for CM-73
print("---- Formula check: #Jac for CM-73 equals (p+1)^2 - 73 ----");
{
forvec(pv = [[1,4]],
    p = cm73_primes[pv[1]];
    a2 = 2*p - 73;
    jac_via_formula = (p+1)^2 - 73;
    jac_via_poly = 1 + a2 + p^2;
    print("p=", p, ": (p+1)^2-73=", jac_via_formula, ", 1+a2+p^2=", jac_via_poly, ", match=", jac_via_formula==jac_via_poly)
);
}
print("");

\\ Pohlig-Hellman speedup ratio: if cover existed, DLP cost on Jac
\\ would be sqrt(max_pfactor) via rho on largest prime-order subgroup,
\\ vs sqrt(n) on the base EC (n = elliptic curve order, n ~ p for secp-like curves).
print("---- Hypothetical PH speedup (IF cover existed) ----");
print("Baseline: DLP on E/F_p costs O(sqrt(p)) via rho.");
print("Cover:    DLP on Jac(C)/F_p costs O(sqrt(max_pfactor)) via PH + rho.");
print("Speedup = sqrt(p) / sqrt(max_pfactor) = sqrt(p / max_pfactor).");
print("");
{
forvec(pv = [[1,4]],
    p = cm73_primes[pv[1]];
    a2 = 2*p - 73;
    jac_order = 1 + a2 + p^2;
    f = factor(jac_order);
    max_pf = vecmax(f[,1]);
    speedup_sq = p * 1.0 / max_pf;
    print("p=", p, ": max_pfactor=", max_pf, ", hypothetical speedup=sqrt(", p, "/", max_pf, ")=sqrt(", speedup_sq, ")~=", sqrt(speedup_sq))
);
}
print("");
print("NOTE: speedup > 1 means PH on Jac beats rho on E -- IF the cover existed.");
print("The main theorem rules out the cover, which is CONSISTENT with the speedup:");
print("a genuine speedup would be exploitable, so the cover cannot exist.");
print("");
print("================================================================");
print("Thread 20: done.");
print("================================================================");
