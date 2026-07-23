\\ Quick supplement: compute gcd(n_i, n_j) for the 4 H3-failing pairs
\\ and also check for all 15 pairs.
\\ Run: gp -q howe_sextic_gcds.gp

default(parisize, 256000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_known = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_known = p + 1 - n_known;

s_sq = (4*p - t_known^2) / 3;
s = sqrtint(s_sq);

t_plus  = (t_known + 3*s) / 2;
t_minus = (t_known - 3*s) / 2;
T = [t_known, t_minus, -(t_known + 3*s)/2, -t_known, (3*s - t_known)/2, t_plus];
N = vector(6, jj, p + 1 - T[jj]);

print("Curve orders:");
for (k = 1, 6, print("  N[", k-1, "] = ", N[k]));
print("");

print("All 15 pairwise gcd(n_i, n_j):");
print("  (i,j) | gcd | factorisation");
print("  ------|-----|---------------");
{
for (i = 1, 6,
    for (j = i+1, 6,
        local(g, f);
        g = gcd(N[i], N[j]);
        f = if(g > 1, factor(g), "1");
        printf("  (%d,%d)  | %s\n", i-1, j-1, Str(g))
    )
);
}
print("");

print("Checking: is n_1 (k=1) prime?  ", isprime(N[2]));
print("Checking: is n_4 (k=4) prime?  ", isprime(N[5]));
print("Checking: isprime for all 6 orders:");
for (k = 1, 6, print("  N[", k-1, "] prime? ", isprime(N[k])));
print("");

\\ Factor the H3-failing orders
print("Factoring orders involved in H3-failures:");
print("N[1] = N[k=1]:"); print(factor(N[2]));
print("N[3] = N[k=3]:"); print(factor(N[4]));
print("N[5] = N[k=5]:"); print(factor(N[6]));
print("N[2] = N[k=2]:"); print(factor(N[3]));
