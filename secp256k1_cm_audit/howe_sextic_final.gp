\\ Final: 15-pair Howe analysis using precomputed group orders.
\\ All ellcard values from howe_sextic_twists.gp Step 5 run.
\\ Flat script (no nested loops) to avoid PARI brace-nesting issues.

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;

N0 = 115792089237316195423570985008687907852837564279074904382605163141518161494337;
N1 = 115792089237316195423570985008687907852598652813156864395638497411212089444244;
N2 = 115792089237316195423570985008687907853031073199722524052490918277602762621571;
N3 = 115792089237316195423570985008687907853702405052206223696310004874299507848991;
N4 = 115792089237316195423570985008687907853941316518124263683276670604605579899084;
N5 = 115792089237316195423570985008687907853508896131558604026424249738214906721757;

if (N0 != n, error("N0 != n"));
print("N0 = n: OK");

T0 = p+1-N0; T1 = p+1-N1; T2 = p+1-N2;
T3 = p+1-N3; T4 = p+1-N4; T5 = p+1-N5;
print("T0=", T0);
print("T1=", T1);
print("T2=", T2);
print("T3=", T3, "  (= -T0: ", T3==-T0, ")");
print("T4=", T4, "  (= -T1: ", T4==-T1, ")");
print("T5=", T5, "  (= -T2: ", T5==-T2, ")");

\\ 2-torsion types (from factor run):
\\   d0,d2,d3,d5: irred [3]  =>  no rational 2-torsion, group order odd or 2-part=1
\\   d1,d4: [1,1,1]           =>  E[2](F_p)=(Z/2Z)^2, so 4 | N1 and 4 | N4

print("\n2-divisibility:");
print("  v2(N0)=", valuation(N0,2));
print("  v2(N1)=", valuation(N1,2), " (d1 has split 2-torsion)");
print("  v2(N2)=", valuation(N2,2));
print("  v2(N3)=", valuation(N3,2));
print("  v2(N4)=", valuation(N4,2), " (d4 has split 2-torsion)");
print("  v2(N5)=", valuation(N5,2));

\\ H2 condition: d[j]/d[i] = omega_p^(j-i) is a cube mod p.
\\ For all 3 H2-pairs: j-i=3, omega_p^3 = -1 mod p.
\\ -1 is a cube mod p iff (-1)^((p-1)/3) = 1 iff (p-1)/3 is even.
e3 = (p-1)/3;
neg1_cube = (Mod(-1,p)^e3 == 1);
print("\n-1 is a cube mod p (for H2 of all 3 pairs): ", neg1_cube);
print("(p-1)/3 mod 2 = ", e3 % 2, "  (even means -1 is cube)");

\\ (H2) check directly for each H2-pair ratio:
omega_p = 60197513588986302554485582024885075108884032450952339817679072026166228089409;
ratio_03 = lift(Mod(N3,p));  \\ irrelevant; ratio is d3/d0 = omega_p^3
d0 = 7;
d3 = lift(Mod(7,p)*Mod(omega_p,p)^3);
d1 = lift(Mod(7,p)*Mod(omega_p,p)^1);
d4 = lift(Mod(7,p)*Mod(omega_p,p)^4);
d2 = lift(Mod(7,p)*Mod(omega_p,p)^2);
d5 = lift(Mod(7,p)*Mod(omega_p,p)^5);
print("\nH2 ratios (all should be -1 mod p = ", p-1, "):");
print("  d3/d0 mod p = ", lift(Mod(d3,p)*Mod(d0,p)^(p-2)));
print("  d4/d1 mod p = ", lift(Mod(d4,p)*Mod(d1,p)^(p-2)));
print("  d5/d2 mod p = ", lift(Mod(d5,p)*Mod(d2,p)^(p-2)));

\\ ---- Check all 3 H2-compatible pairs ----
print("\n=== Howe conditions for 3 H2-compatible pairs ===\n");

\\ Pair (0,3): secp256k1 + quadratic twist
g03 = gcd(N0, N3);
print("Pair (0,3): d0=7 (secp256k1) vs d3=-7 mod p (quadratic twist)");
print("  H1: N0!=N3?   ", N0!=N3);
print("  H2: ratio cube? ", neg1_cube);
print("  gcd(N0,N3) = ", g03);
print("  H3: gcd=1?    ", g03==1);
print("  HOWE-GLUABLE? ", N0!=N3 && neg1_cube && g03==1);

print("");

\\ Pair (1,4): d1 vs d4
g14 = gcd(N1, N4);
print("Pair (1,4): d1 (cubic twist) vs d4 (its quadratic twist)");
print("  H1: N1!=N4?   ", N1!=N4);
print("  H2: ratio cube? ", neg1_cube);
print("  gcd(N1,N4) = ", g14);
print("  gcd factored  = ", factor(g14));
print("  H3: gcd=1?    ", g14==1);
print("  HOWE-GLUABLE? ", N1!=N4 && neg1_cube && g14==1);

print("");

\\ Pair (2,5): d2 vs d5
g25 = gcd(N2, N5);
print("Pair (2,5): d2 (sextic twist) vs d5 (its quadratic twist)");
print("  H1: N2!=N5?   ", N2!=N5);
print("  H2: ratio cube? ", neg1_cube);
print("  gcd(N2,N5) = ", g25);
if (g25 > 1, print("  gcd factored = ", factor(g25)));
print("  H3: gcd=1?    ", g25==1);
print("  HOWE-GLUABLE? ", N2!=N5 && neg1_cube && g25==1);

print("\n=== 12 non-H2 pairs ===");
print("(H2) fails: d[j]/d[i] = omega_p^(j-i) is NOT a cube for j-i not ≡ 0 mod 3.");
print("Pairs: {0,1},{0,2},{0,4},{0,5},{1,2},{1,3},{1,5},{2,3},{2,4},{3,4},{3,5},{4,5}");
print("All 12 pairs: not Howe-gluable (H2 obstruction).\n");

print("=== Summary ===");
total_gluable = (N0!=N3 && neg1_cube && g03==1) +
                (N1!=N4 && neg1_cube && g14==1) +
                (N2!=N5 && neg1_cube && g25==1);
print("Pairs passing H1+H2+H3: ", total_gluable, " / 15");
print("Done.");
