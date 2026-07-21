\\ Thread 18: Verify Howe-glueable pairs using ellcard for exact curve orders.
\\ Uses PARI's built-in ellcard (SEA algorithm) for each twist.
\\ Run: gp -q thread18_ellcard.gp

default(parisize, 1024000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_secp = p + 1 - n_secp;
s = sqrtint((4*p - t_secp^2)/3);

g_gen = znprimroot(p);
u = lift(g_gen ^ ((p - 1) / 6));

print("p = ", p);
print("t_secp = ", t_secp);
print("s = ", s);
print("");

\\ Six b-values (k=0..5)
B = vector(6);
B[1] = 7;
forstep(k = 1, 5, 1, B[k+1] = lift(Mod(7,p) * Mod(u,p)^k));

print("--- b-values ---");
for (k = 1, 6, print("k=", k-1, ": b = ", B[k]));
print("");

\\ Six CM candidate traces (abstract)
T_cands = [t_secp, (t_secp-3*s)/2, -(t_secp+3*s)/2, -t_secp, (3*s-t_secp)/2, (t_secp+3*s)/2];
print("--- N mod 4 for each candidate trace ---");
for (i = 1, 6, print("T_cands[", i, "]:  N = p+1-T = ... mod 4 = ", (p+1-T_cands[i])%4));
print("");

\\ Compute actual orders using PARI's ellcard (SEA)
print("--- Computing exact curve orders via ellcard ---");
print("(This may take up to 30s per curve for 256-bit p)");
N_exact = vector(6);
for (k = 1, 6,
    Ek = ellinit([0, B[k]], p);
    N_exact[k] = ellcard(Ek);
    print("k=", k-1, ": #E = ", N_exact[k], "  mod 4 = ", N_exact[k] % 4)
);
print("");

\\ Match exact orders to candidate traces
print("--- Matching exact orders to CM traces ---");
T_exact = vector(6);
for (k = 1, 6,
    Tk = p + 1 - N_exact[k];
    T_exact[k] = Tk;
    \\ Find which candidate this matches
    match = 0;
    for (ci = 1, 6,
        if (T_cands[ci] == Tk, match = ci)
    );
    print("k=", k-1, ": actual trace = ", Tk, "  matches T_cands[", match, "]")
);
print("");

\\ Recompute deg_pat (2-torsion patterns)
print("--- 2-torsion patterns ---");
deg_pat = vector(6);
for (k = 1, 6,
    f = factor(Mod(1,p) * (x^3 + B[k]));
    nr = matsize(f)[1];
    pat = [];
    for (r = 1, nr, pat = concat(pat, [poldegree(f[r,1])]));
    deg_pat[k] = pat;
    print("k=", k-1, ": deg-pattern = ", pat)
);
print("");

\\ Compute all 15 Howe conditions with CORRECT orders
print("================================================================");
print("All 15 pairwise Howe conditions — CORRECT orders from ellcard");
print("================================================================");
print("");
print("  (i,j) | H1: n_i≠n_j | H2: same 2-tor | H3: gcd=1 | Glueable?");
print("  ------|-------------|----------------|-----------|----------");
count_yes = 0;
glueable_pairs = [];
for (i = 1, 6,
    for (j = i+1, 6,
        ni = N_exact[i];
        nj = N_exact[j];
        h1 = (ni != nj);
        h2 = (deg_pat[i] == deg_pat[j]);
        h3 = (gcd(ni, nj) == 1);
        gl = h1 && h2 && h3;
        print("  (", i-1, ",", j-1, ")  | ",
              if(h1,"YES","NO "), " | ",
              if(h2,"YES","NO "), " | ",
              if(h3,"YES","NO "), " | ",
              if(gl,"YES <- Howe-glueable","no"));
        if (gl, count_yes++; glueable_pairs = concat(glueable_pairs, [[i-1,j-1]]))
    )
);
print("");
print("Howe-glueable pairs (ellcard-verified): ", count_yes, " / 15");
if (#glueable_pairs > 0,
    for (r = 1, #glueable_pairs, print("  (", glueable_pairs[r][1], ",", glueable_pairs[r][2], ")"))
);
print("");

\\ Check consistency with the howe_sextic_twists_all15.gp result
\\ That script found: (0,2), (0,3), (0,5), (1,4), (2,3)
print("Script's answer (assumed ordering): (0,2),(0,3),(0,5),(1,4),(2,3)");
print("This run's answer (ellcard-verified): ", glueable_pairs);
print("");
print("Done.");
