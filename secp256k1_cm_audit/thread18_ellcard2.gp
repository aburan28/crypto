\\ Thread 18: Compute exact orders for k=1 and k=4 using ellcard.
\\ Determine correct trace assignment; recompute all H3 gcds.
\\ No nested loops -- avoids PARI embedded-brace restriction.
\\ Run: gp -q thread18_ellcard2.gp

default(parisize, 1024000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_secp = p + 1 - n_secp;
s = sqrtint((4*p - t_secp^2)/3);

g_gen = znprimroot(p);
u = lift(g_gen ^ ((p - 1) / 6));

b0 = 7;
b1 = lift(Mod(7,p) * Mod(u,p)^1);
b2 = lift(Mod(7,p) * Mod(u,p)^2);
b3 = p - 7;
b4 = lift(Mod(7,p) * Mod(u,p)^4);
b5 = lift(Mod(7,p) * Mod(u,p)^5);

\\ Known orders
N0 = n_secp;
N3 = p + 1 + t_secp;

\\ Abstract CM traces (6 of them)
T1 = t_secp;
T2 = (t_secp - 3*s) / 2;
T3 = -(t_secp + 3*s) / 2;
T4 = -t_secp;
T5 = (3*s - t_secp) / 2;
T6 = (t_secp + 3*s) / 2;

print("=== CM trace mod-4 check ===");
print("T3 -> N = p+1-T3, mod 4 = ", (p+1-T3) % 4, "  (must be 0 for [1,1,1] twists)");
print("T6 -> N = p+1-T6, mod 4 = ", (p+1-T6) % 4);
print("T2 -> N = p+1-T2, mod 4 = ", (p+1-T2) % 4, "  (must be odd for [3] twists)");
print("T5 -> N = p+1-T5, mod 4 = ", (p+1-T5) % 4);
print("");
print("k=1 and k=4 have [1,1,1] 2-torsion (x^3+b splits completely).");
print("=> their traces must be T3 and T6 (in some order).");
print("k=2 and k=5 have [3] 2-torsion (x^3+b irreducible).");
print("=> their traces must be T2 and T5 (in some order).");
print("");

\\ Compute #E_{k=1} using ellcard (SEA algorithm, ~10s for 256-bit p)
print("=== Computing ellcard for k=1 (may take ~10-30s) ===");
E1 = ellinit([0, b1], p);
N1 = ellcard(E1);
print("k=1: #E = ", N1);
print("  trace = ", p + 1 - N1);
print("  mod 4 = ", N1 % 4);
print("  trace matches T3? ", (p + 1 - N1) == T3);
print("  trace matches T6? ", (p + 1 - N1) == T6);
print("");

\\ From k=1, derive k=4 (quadratic twist: trace -> -trace)
T_k1 = p + 1 - N1;
T_k4 = -T_k1;
N4 = p + 1 - T_k4;
print("k=4 (quad-twist of k=1): trace = ", T_k4, "  #E = ", N4);
print("");

\\ Similarly for k=2: get trace of k=2 via ellcard
print("=== Computing ellcard for k=2 (may take ~10-30s) ===");
E2 = ellinit([0, b2], p);
N2 = ellcard(E2);
print("k=2: #E = ", N2);
print("  trace = ", p + 1 - N2);
print("  trace matches T2? ", (p + 1 - N2) == T2);
print("  trace matches T5? ", (p + 1 - N2) == T5);
print("");

T_k2 = p + 1 - N2;
T_k5 = -T_k2;
N5 = p + 1 - T_k5;
print("k=5 (quad-twist of k=2): trace = ", T_k5, "  #E = ", N5);
print("");

\\ Now we have all 6 orders: N0, N1, N2, N3, N4, N5.
print("=== All 6 orders (ellcard-verified) ===");
print("k=0: #E = ", N0, "  mod 4 = ", N0 % 4);
print("k=1: #E = ", N1, "  mod 4 = ", N1 % 4);
print("k=2: #E = ", N2, "  mod 4 = ", N2 % 4);
print("k=3: #E = ", N3, "  mod 4 = ", N3 % 4);
print("k=4: #E = ", N4, "  mod 4 = ", N4 % 4);
print("k=5: #E = ", N5, "  mod 4 = ", N5 % 4);
print("");

\\ GCDs for all relevant pairs
\\ [3] x [3] pairs: (0,2),(0,3),(0,5),(2,3),(2,5),(3,5)
\\ [1,1,1] x [1,1,1] pair: (1,4)
print("=== H3 checks (gcd = 1?) for same-2-torsion pairs ===");
print("[3] x [3] pairs:");
print("  gcd(N0, N2) [pair (0,2)] = ", gcd(N0, N2));
print("  gcd(N0, N3) [pair (0,3)] = ", gcd(N0, N3));
print("  gcd(N0, N5) [pair (0,5)] = ", gcd(N0, N5));
print("  gcd(N2, N3) [pair (2,3)] = ", gcd(N2, N3));
print("  gcd(N2, N5) [pair (2,5)] = ", gcd(N2, N5));
print("  gcd(N3, N5) [pair (3,5)] = ", gcd(N3, N5));
print("[1,1,1] x [1,1,1] pair:");
print("  gcd(N1, N4) [pair (1,4)] = ", gcd(N1, N4));
print("");

print("=== Correct Howe-glueable pairs ===");
v02 = gcd(N0,N2); v03 = gcd(N0,N3); v05 = gcd(N0,N5);
v23 = gcd(N2,N3); v25 = gcd(N2,N5); v35 = gcd(N3,N5); v14 = gcd(N1,N4);
print("(0,2): ", if(v02==1,"GLUEABLE","not (gcd=",v02,")"));
print("(0,3): ", if(v03==1,"GLUEABLE","not (gcd=",v03,")"));
print("(0,5): ", if(v05==1,"GLUEABLE","not (gcd=",v05,")"));
print("(2,3): ", if(v23==1,"GLUEABLE","not (gcd=",v23,")"));
print("(2,5): ", if(v25==1,"GLUEABLE","not (gcd=",v25,")"));
print("(3,5): ", if(v35==1,"GLUEABLE","not (gcd=",v35,")"));
print("(1,4): ", if(v14==1,"GLUEABLE","not (gcd=",v14,")"));
print("");
print("Done.");
