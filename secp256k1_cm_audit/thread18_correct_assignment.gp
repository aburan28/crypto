\\ Thread 18: Verify correct CM trace assignment for secp256k1 sextic twists.
\\ The original howe_sextic_twists_all15.gp assumed T[k+1] <-> b_k (0-based k)
\\ without verification. This script verifies the assignment using scalar mult.
\\ Run: gp -q thread18_correct_assignment.gp

default(parisize, 1024000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t = p + 1 - n_secp;
s = sqrtint((4*p - t^2)/3);

g_gen = znprimroot(p);
u = lift(g_gen ^ ((p - 1) / 6));

\\ Six b-values (0-based k):
b = vector(6, k, lift(Mod(7, p) * Mod(u, p)^(k-1)));
print("b values (k=0..5): ", b);
print("");

\\ Six candidate traces (abstract, from CM formula):
T_all = [t, (t-3*s)/2, -(t+3*s)/2, -t, (3*s-t)/2, (t+3*s)/2];
\\ Known: k=0 -> T=t=T_all[1], k=3 -> T=-t=T_all[4].

\\ For each of k=1,2,4,5, find a random F_p-point and test candidate traces.
\\ Candidates for k=1: T_all[3] or T_all[6] (the two giving 4|N).
\\ Candidates for k=2: T_all[2] or T_all[5].
\\ k=4 = quadratic twist of k=1, k=5 = quadratic twist of k=2.

print("=== Determining correct trace for k=1 ===");
bk = b[2];  \\ k=1
print("b_{k=1} = ", bk);

\\ Find a point on E: y^2 = x^3 + bk
xP = -1; yP = -1;
forstep(x = 1, 10^6, 1,
    rhs = lift(Mod(x, p)^3 + Mod(bk, p));
    if (kronecker(rhs, p) == 1,
        yP = lift(sqrtmod(rhs, p));
        xP = x;
        break
    )
);
if (xP < 0, error("No point found"));
print("Test point: (x=", xP, ", y=", yP, ")");
print("Verify: y^2 = x^3 + b? ", lift(Mod(yP,p)^2 - Mod(xP,p)^3 - Mod(bk,p)) == 0);

\\ Test both candidate traces
Ek1 = ellinit([Mod(0,p), Mod(bk,p)]);
Pk1 = ellpoint(Ek1, Mod(xP,p), Mod(yP,p));
for (ci = 1, 2,
    Tcand = T_all[3 + (ci-1)*3];  \\ T_all[3] or T_all[6]
    Ncand = p + 1 - Tcand;
    res = ellmul(Ek1, Pk1, Ncand);
    print("T_cand=T_all[", 3+(ci-1)*3, "]=", Tcand, " -> N=", Ncand,
          " -> N*P = ", if (res[1] == 0 && res[2] == 0, "O  <- CORRECT TRACE", Str(res)))
);
print("");

print("=== Determining correct trace for k=2 ===");
bk2 = b[3];  \\ k=2
print("b_{k=2} = ", bk2);
xP2 = -1; yP2 = -1;
forstep(x = 1, 10^6, 1,
    rhs = lift(Mod(x, p)^3 + Mod(bk2, p));
    if (kronecker(rhs, p) == 1,
        yP2 = lift(sqrtmod(rhs, p));
        xP2 = x;
        break
    )
);
print("Test point: (x=", xP2, ", y=", yP2, ")");
Ek2 = ellinit([Mod(0,p), Mod(bk2,p)]);
Pk2 = ellpoint(Ek2, Mod(xP2,p), Mod(yP2,p));
for (ci = 1, 2,
    Tcand = T_all[2 + (ci-1)*3];  \\ T_all[2] or T_all[5]
    Ncand = p + 1 - Tcand;
    res = ellmul(Ek2, Pk2, Ncand);
    print("T_cand=T_all[", 2+(ci-1)*3, "]=", Tcand, " -> N=", Ncand,
          " -> N*P = ", if (res[1] == 0 && res[2] == 0, "O  <- CORRECT TRACE", Str(res)))
);
print("");

print("=== Summary of correct trace assignment ===");
print("(k=0 and k=3 are known from secp256k1 data)");
print("k=0: T =  t  =  ", t, "  N = ", n_secp);
print("k=3: T = -t  = ", -t, "  N = ", p+1+t);
print("k=4 = quad-twist of k=1  ->  T_{k=4} = -T_{k=1}");
print("k=5 = quad-twist of k=2  ->  T_{k=5} = -T_{k=2}");

\\ After confirming k=1 and k=2 traces, recompute glueable pairs.
print("");
print("=== Recomputing Howe conditions with CORRECT trace assignment ===");
\\ We'll fill in the correct N after the trace tests above.
\\ For now, record all 6 gcd values for [3]x[3] and [1,1,1]x[1,1,1] pairs:
N_all = vector(6, k, p + 1 - T_all[k]);

\\ The script's N_all[k] and the CORRECT N for twist index k-1 differ for k=2,3,5,6.
\\ After verification: let N_correct[k] = order of E_{k-1}:
\\ (We'll assign after scalar-mult test confirms T.)
print("Will update after scalar-mult test results above.");
print("");
print("Direct GCD checks (all 15 pairs) using CORRECT orders:");
print("(Orders to be verified by scalar-mult tests above)");
print("Once correct assignment is determined, gcd table follows.");
\\ (Script continues after manual confirmation)
