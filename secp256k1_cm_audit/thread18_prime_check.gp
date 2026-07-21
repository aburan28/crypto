\\ Thread 18 companion: primality and GCD analysis
default(parisize, 1024000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_known = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_known = p + 1 - n_known;

s_sq = (4*p - t_known^2) / 3;
s = sqrtint(s_sq);
T = [t_known, (t_known - 3*s)/2, -(t_known + 3*s)/2, -t_known, (3*s - t_known)/2, (t_known + 3*s)/2];
N = vector(6, jj, p + 1 - T[jj]);

print("=== Six twist orders and primality ===");
{
local(primecount);
primecount = 0;
for (kk = 1, 6,
    local(prim);
    prim = isprime(N[kk]);
    print("n_", kk-1, " prime? ", prim, "  order = ", N[kk]);
    if (prim, primecount++)
);
print("Total prime orders: ", primecount, " / 6");
}
print("");

print("=== GCD for all 15 pairs ===");
{
for (ii = 1, 6,
    for (jj = ii+1, 6,
        local(gg);
        gg = gcd(N[ii], N[jj]);
        if (gg > 1,
            print("(", ii-1, ",", jj-1, "): gcd = ", gg, "  *** SHARED FACTOR ***"),
            print("(", ii-1, ",", jj-1, "): gcd = 1")
        )
    )
);
}
print("");

print("=== Factor non-prime orders ===");
{
for (kk = 1, 6,
    if (!isprime(N[kk]),
        print("n_", kk-1, " factors:");
        print("  ", factor(N[kk]))
    )
);
}
print("");

\\ Glueable pairs cross-check
print("=== Glueable pair orders (pairs that pass H1+H2+H3) ===");
{
\\ Pairs passing all 3 from the main script: (0,2),(0,3),(0,5),(1,4),(2,3)
\\ In 1-indexed: (1,3),(1,4),(1,6),(2,5),(3,4)
glueable1 = [1, 1, 1, 2, 3];
glueable2 = [3, 4, 6, 5, 4];
for (rr = 1, 5,
    local(ii, jj, g12);
    ii = glueable1[rr];
    jj = glueable2[rr];
    g12 = gcd(N[ii], N[jj]);
    print("Pair (", ii-1, ",", jj-1, "): n_", ii-1, " prime=", isprime(N[ii]),
          "  n_", jj-1, " prime=", isprime(N[jj]),
          "  gcd=", g12)
);
}
print("");

\\ Note: the p = n_0 case (secp256k1 is prime-order) is known
print("Note: n_0 = n_secp256k1 = ", n_known, " (prime by construction)");
print("Verify N[1] == n_known: ", N[1] == n_known);
