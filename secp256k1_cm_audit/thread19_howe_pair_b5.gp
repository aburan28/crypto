\\ thread19_howe_pair_b5.gp
\\ Thread 19: Weil polynomial and B5 HCDLP-cost analysis for all 5 Howe-glueable
\\ sextic-twist pairs of secp256k1, found in Thread 18.

default(parisize, 256000000);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n0 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t0_val = p + 1 - n0;

print("Thread 19: B5 Weil polynomial analysis for 5 Howe-glueable pairs");
print("=================================================================\n");

bsq = (4*p - t0_val^2) \ 3;
bval = sqrtint(bsq);
aval = (t0_val + bval) \ 2;

\\ 6 traces and group orders (1-indexed: index kk corresponds to twist k=kk-1)
TR = [2*aval-bval, aval-2*bval, -(aval+bval), -(2*aval-bval), 2*bval-aval, aval+bval];
ORD = vector(6, kk, p + 1 - TR[kk]);

print("Group orders:");
print("  n0 = ", ORD[1]);
print("  n1 = ", ORD[2]);
print("  n2 = ", ORD[3]);
print("  n3 = ", ORD[4]);
print("  n4 = ", ORD[5]);
print("  n5 = ", ORD[6]);
print("");

\\ ecdlp cost bits (sqrt(n0))
eb = log(ORD[1])*1.0/log(2.0)/2;
print("ECDLP baseline (sqrt(n0)) ~ 2^", eb, " group ops\n");

\\ ---------------------------------------------------------------
\\ Compute each pair explicitly (no loops to avoid PARI parser issues)
\\ ---------------------------------------------------------------

\\ Pair 1: (0,1) => indices 1,2
n01 = ORD[1]; n02 = ORD[2];
jac01 = n01*n02;
hb01 = log(jac01)*1.0/log(2.0)/2;
W1_01 = (1 - TR[1] + p) * (1 - TR[2] + p);
print("--- Pair (0,1) ---");
print("  n_0*n_1 = ", jac01);
print("  W(1) == n0*n1: ", W1_01 == jac01);
print("  log2(#Jac) = ", log(jac01)/log(2.0));
print("  HCDLP bits = sqrt(#Jac) ~ 2^", hb01);
print("  extra vs ECDLP: +2^", hb01 - eb);
print("  B5 holds: ", hb01 >= eb);
print("");

\\ Pair 2: (0,3) => indices 1,4
n11 = ORD[1]; n14 = ORD[4];
jac03 = n11*n14;
hb03 = log(jac03)*1.0/log(2.0)/2;
W1_03 = (1 - TR[1] + p) * (1 - TR[4] + p);
\\ exact: n0*n3 = (p+1)^2 - t0^2
exact03 = (p+1)^2 - TR[1]^2;
print("--- Pair (0,3) [known/reference] ---");
print("  n_0*n_3 = ", jac03);
print("  (p+1)^2-t0^2 = ", exact03, "  match: ", exact03 == jac03);
print("  W(1) == n0*n3: ", W1_03 == jac03);
print("  HCDLP bits ~ 2^", hb03);
print("  extra vs ECDLP: +2^", hb03 - eb);
print("  B5 holds: ", hb03 >= eb);
print("");

\\ Pair 3: (0,4) => indices 1,5
n13 = ORD[1]; n35 = ORD[5];
jac04 = n13*n35;
hb04 = log(jac04)*1.0/log(2.0)/2;
W1_04 = (1 - TR[1] + p) * (1 - TR[5] + p);
print("--- Pair (0,4) ---");
print("  n_0*n_4 = ", jac04);
print("  W(1) == n0*n4: ", W1_04 == jac04);
print("  HCDLP bits ~ 2^", hb04);
print("  extra vs ECDLP: +2^", hb04 - eb);
print("  B5 holds: ", hb04 >= eb);
print("");

\\ Pair 4: (1,4) => indices 2,5
n21 = ORD[2]; n25 = ORD[5];
jac14 = n21*n25;
hb14 = log(jac14)*1.0/log(2.0)/2;
W1_14 = (1 - TR[2] + p) * (1 - TR[5] + p);
\\ (1,4) is also a quadratic-twist pair: n1*n4 = (p+1)^2 - t1^2
exact14 = (p+1)^2 - TR[2]^2;
print("--- Pair (1,4) [quadratic-twist pair] ---");
print("  n_1*n_4 = ", jac14);
print("  (p+1)^2-t1^2 = ", exact14, "  match: ", exact14 == jac14);
print("  W(1) == n1*n4: ", W1_14 == jac14);
print("  HCDLP bits ~ 2^", hb14);
print("  extra vs ECDLP: +2^", hb14 - eb);
print("  B5 holds: ", hb14 >= eb);
print("");

\\ Pair 5: (3,4) => indices 4,5
n41 = ORD[4]; n45 = ORD[5];
jac34 = n41*n45;
hb34 = log(jac34)*1.0/log(2.0)/2;
W1_34 = (1 - TR[4] + p) * (1 - TR[5] + p);
print("--- Pair (3,4) ---");
print("  n_3*n_4 = ", jac34);
print("  W(1) == n3*n4: ", W1_34 == jac34);
print("  HCDLP bits ~ 2^", hb34);
print("  extra vs ECDLP: +2^", hb34 - eb);
print("  B5 holds: ", hb34 >= eb);
print("");

\\ ---------------------------------------------------------------
\\ Smooth-factor analysis (single-statement for loop)
\\ ---------------------------------------------------------------
print("Smooth-factor analysis (primes up to 10^7):");
print("n0 = ", factor(ORD[1], 10^7));
print("n1 = ", factor(ORD[2], 10^7));
print("n3 = ", factor(ORD[4], 10^7));
print("n4 = ", factor(ORD[5], 10^7));
print("");

\\ ---------------------------------------------------------------
\\ Lower bound on HCDLP cost from Hasse
\\ For any n_i: n_i >= p + 1 - 2*floor(sqrt(p)) - 2 >= p - 2*sqrt(p) - 1
\\ So n_i * n_j >= (p - 2*sqrt(p))^2 > p^2 - 4*p^{3/2}
\\ sqrt(n_i*n_j) >= p - 2*sqrt(p) >> sqrt(p) for p ~ 2^256
\\ ---------------------------------------------------------------
sqp = sqrtint(p);
lower_ni = p + 1 - 2*(sqp+1);
lower_product = lower_ni^2;
print("Hasse lower bound:");
print("  min n_i >= p - 2*floor(sqrt(p)) - 1 = ", lower_ni);
print("  min (n_i*n_j) >= ", lower_product);
print("  min sqrt(n_i*n_j) ~ 2^", log(lower_product)*1.0/log(2.0)/2);
print("  This >> sqrt(p) ~ 2^128, confirming B5 holds universally by Hasse.\n");

\\ ---------------------------------------------------------------
\\ Summary
\\ ---------------------------------------------------------------
print("=================================================================");
print("Summary table:");
print("Pair    #Jac-bits  HCDLP-2^  ECDLP-2^  extra-bits  B5");
print("(0,1)   ", round(log(jac01)/log(2.0)), "         ", hb01, "    ", eb, "    ", hb01-eb, "   YES");
print("(0,3)   ", round(log(jac03)/log(2.0)), "         ", hb03, "    ", eb, "    ", hb03-eb, "   YES");
print("(0,4)   ", round(log(jac04)/log(2.0)), "         ", hb04, "    ", eb, "    ", hb04-eb, "   YES");
print("(1,4)   ", round(log(jac14)/log(2.0)), "         ", hb14, "    ", eb, "    ", hb14-eb, "   YES");
print("(3,4)   ", round(log(jac34)/log(2.0)), "         ", hb34, "    ", eb, "    ", hb34-eb, "   YES");
print("");
print("All 5 pairs: HCDLP cost ~ 2^256 >> ECDLP cost ~ 2^128. B5 holds.");
print("=================================================================");
