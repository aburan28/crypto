\\ thread22_z3z_qualifying_pairs.gp
\\ Thread 22: Z/3Z search for each qualifying Howe pair, proxy p=43

default(parisize, 256000000);
default(timer, 0);

p = 43;

print("Thread 22 -- proxy p=", p);
print("================================================================");
print("");

\\ CM: E0 y^2=x^3+7, t0=13, s=1 (precomputed: 4*43=172=169+3)
t0 = 13;
s0 = 1;
T_vec = [t0, (t0-3*s0)/2, -(t0+3*s0)/2, -t0, (3*s0-t0)/2, (t0+3*s0)/2];
N_vec = vector(6, k, p + 1 - T_vec[k]);
print("Traces: ", T_vec);
print("Orders: ", N_vec);
print("");

\\ b-values: b_k = 7 * u^(k-1) mod 43 where u = prim 6th root
u = lift(znprimroot(p) ^ ((p - 1) / 6));
print("u = ", u);
b_vec = vector(6, k, lift(Mod(7, p) * Mod(u, p)^(k-1)));
print("b: ", b_vec);
print("");

\\ 2-torsion degree patterns for x^3 + b_k over F_43
dp = vector(6);
{
for (kk = 1, 6,
    local(f1, fac1, nr1, pat1);
    f1 = Mod(1, p) * (x^3 + b_vec[kk]);
    fac1 = factor(f1);
    nr1 = matsize(fac1)[1];
    pat1 = [];
    for (i = 1, nr1, pat1 = concat(pat1, [poldegree(fac1[i,1])]));
    dp[kk] = pat1;
    print("  k=", kk-1, " b=", b_vec[kk], " pat=", pat1)
);
}
print("");

\\ Build list of qualifying pairs
q_ii = List();
q_jj = List();
print("Qualifying pairs:");
{
for (ii = 1, 6,
    for (jj = ii+1, 6,
        local(h1, h2, h3);
        h1 = (N_vec[ii] != N_vec[jj]);
        h2 = (dp[ii] == dp[jj]);
        h3 = (gcd(N_vec[ii], N_vec[jj]) == 1);
        if (h1 && h2 && h3,
            listput(q_ii, ii);
            listput(q_jj, jj);
            print("  (", ii-1, ",", jj-1, ") t=(", T_vec[ii], ",", T_vec[jj],
                  ") n=(", N_vec[ii], ",", N_vec[jj], ")")
        )
    )
);
}
print("Total qualifying: ", #q_ii);
print("");

\\ For each qualifying pair, search y^2 = x^6 + a*x^3 + b
print("================================================================");
print("Z/3Z search for each qualifying pair");
print("================================================================");
{
for (pidx = 1, #q_ii,
    local(qi, qj, ti, tj, tcp, nsm, nmatch);
    qi = q_ii[pidx];
    qj = q_jj[pidx];
    ti = T_vec[qi];
    tj = T_vec[qj];
    tcp = expand((x^2 - ti*x + p) * (x^2 - tj*x + p));
    print("Pair (", qi-1, ",", qj-1, "): t=(", ti, ",", tj, ")  target=", tcp);
    nsm = 0;
    nmatch = 0;
    for (a = 0, p-1,
        for (b = 0, p-1,
            local(f2, dm, cp2);
            f2 = x^6 + a*x^3 + b;
            dm = lift(poldisc(f2) * Mod(1, p));
            if (dm == 0, next);
            nsm++;
            cp2 = hyperellcharpoly(f2 * Mod(1, p));
            if (cp2 == tcp,
                nmatch++;
                if (nmatch <= 3, print("  match: a=", a, " b=", b))
            )
        )
    );
    print("  smooth=", nsm, "  matches=", nmatch);
    if (nmatch > 0,
        print("  -> POSITIVE"),
        print("  -> NEGATIVE")
    );
    print("")
);
}

print("================================================================");
print("Done.");
