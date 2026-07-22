\\ ============================================================
\\ Thread 18: Sextic twist analysis — gcd structure and order factorization
\\ ============================================================
\\ Run: gp -q thread18_sextic_analysis.gp

default(parisize, 256000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
t = 432420386565659656852420866390673177327;
s_sq = (4*p - t^2)/3;
s = sqrtint(s_sq);

T = [t, (t - 3*s)/2, -(t + 3*s)/2, -t, (3*s - t)/2, (t + 3*s)/2];
N = vector(6, k, p + 1 - T[k]);

small_primes = [2, 3, 4, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
sp = #small_primes;

print("==== Thread 18: Small factor analysis of twist orders ====");
print("");
print("Small prime divisors of each twist order N_k (trial division):");
{
for (k=1,6,
    local(n, small_facs, rem, q, bits);
    n = N[k];
    small_facs = [];
    rem = n;
    for (qi=1, sp,
        q = small_primes[qi];
        if (rem % q == 0,
            small_facs = concat(small_facs, [q]);
            rem = rem \ q
        )
    );
    bits = floor(log(rem) / log(2));
    print("  k=", k-1, ": small divisors=", small_facs, "  cofactor~2^", bits)
);
}
print("");

print("==== GCD matrix for all 15 pairs ====");
print("");
print("  (i,j) | gcd | H2? | H3? | Glueable?");
print("  ------|-----|-----|-----|----------");
deg_pats = [[3],[1,1,1],[3],[3],[1,1,1],[3]]; \\ from main script output
{
for (i=1,6,
    for (j=i+1,6,
        local(g, h1, h2, h3, glueable);
        g = gcd(N[i], N[j]);
        h1 = (N[i] != N[j]);
        h2 = (deg_pats[i] == deg_pats[j]);
        h3 = (g == 1);
        glueable = h1 && h2 && h3;
        printf("  (%d,%d)  | %d | %s | %s | %s\n",
            i-1, j-1, g,
            if(h2,"YES","NO "),
            if(h3,"YES","NO "),
            if(glueable,"YES <- glueable","no"))
    )
);
}
print("");

print("==== Structural pattern: which twists share small factors ====");
print("");
print("3 | N_k for k in:");
{
for (k=1,6,
    if (N[k] % 3 == 0, print("  k=",k-1,": 3 | N_k"))
);
}
print("4 | N_k for k in:");
{
for (k=1,6,
    if (N[k] % 4 == 0, print("  k=",k-1,": 4 | N_k"))
);
}
print("2 | N_k for k in:");
{
for (k=1,6,
    if (N[k] % 2 == 0, print("  k=",k-1,": 2 | N_k"))
);
}
print("");

print("==== Glueable pairs summary ====");
print("");
count_yes = 0;
{
for (i=1,6,
    for (j=i+1,6,
        local(g, h1, h2, h3, glueable);
        g = gcd(N[i], N[j]);
        h1 = (N[i] != N[j]);
        h2 = (deg_pats[i] == deg_pats[j]);
        h3 = (g == 1);
        glueable = h1 && h2 && h3;
        if (glueable, {count_yes++; print("  (",i-1,",",j-1,") GLUEABLE")})
    )
);
}
print("Total glueable: ", count_yes, " / 15");
print("");
print("secp256k1 (k=0) is in how many glueable pairs?");
{
local(cnt0);
cnt0 = 0;
for (j=1,6,
    if (j == 1, next);
    local(g, h2, h3);
    g = gcd(N[1], N[j]);
    h2 = (deg_pats[1] == deg_pats[j]);
    h3 = (g == 1);
    if ((N[1] != N[j]) && h2 && h3, {cnt0++; print("  (0,",j-1,") GLUEABLE with secp256k1")})
);
print("  secp256k1 participates in ", cnt0, " glueable pairs");
}
print("");
print("Done.");
