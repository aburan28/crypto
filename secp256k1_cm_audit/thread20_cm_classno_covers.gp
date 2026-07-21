\\ Thread 20: CM class numbers for genus-2 Jacobians of secp256k1 glueable pairs
\\ Fixed for PARI 2.15.4 (no nested if/for bodies with braces).
\\
\\ Hypothesis: for ANY j=0 sextic twist pair, ALL biquadratic glueable Jacobians
\\ have CM field Q(sqrt(-3)) with h=1, making [P]^2=1 trivially satisfied.
\\
\\ Run: gp -q thread20_cm_classno_covers.gp

default(parisize, 2048 * 1024 * 1024);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t0 = p + 1 - n_secp;

print("=== secp256k1 CM trace ===");
print("t0 = ", t0);

\\ Cornacchia: 4p = t0^2 + 3*u0^2 (valid for j=0 CM)
u0_sq = (4*p - t0^2) / 3;
u0 = sqrtint(u0_sq);
print("u0 = ", u0);
print("Check 4p = t0^2 + 3*u0^2: ", if(4*p == t0^2 + 3*u0^2, "PASS", "FAIL"));
print("");

\\ The 6 sextic-twist traces for j=0 over F_p are:
\\   t0, -t0, ta, -ta, tb, -tb
\\   where ta = (t0 + 3*u0)/2, tb = (t0 - 3*u0)/2
\\ (t0 and 3*u0 are both odd since t0 odd, u0 odd => t0+3*u0 = even => /2 is integer)
ta = (t0 + 3*u0) / 2;
tb = (t0 - 3*u0) / 2;

print("=== Six sextic-twist traces ===");
print("k=0: t0 = ", t0, "  N=", p+1-t0);
print("k=1: ta = ", ta, "  N=", p+1-ta);
print("k=2: tb = ", tb, "  N=", p+1-tb);
print("k=3: -t0 = ", -t0, "  N=", p+1+t0);
print("k=4: -ta = ", -ta, "  N=", p+1+ta);
print("k=5: -tb = ", -tb, "  N=", p+1+tb);
print("");

\\ Cornacchia check for ta and tb: t_k^2 + 3*u_k^2 = 4p ?
diff_ta = 4*p - ta^2;
diff_tb = 4*p - tb^2;
ua = sqrtint(diff_ta / 3);
ub = sqrtint(diff_tb / 3);
print("Cornacchia check for ta: ", if(ta^2 + 3*ua^2 == 4*p, "PASS (ua=", "FAIL"));
print("  ua = ", ua, ", check: ", ta^2 + 3*ua^2 == 4*p);
print("Cornacchia check for tb: ua=", ub, ", check: ", tb^2 + 3*ub^2 == 4*p);
print("");

\\ ====== BIQUADRATIC GLUEABLE PAIRS: (0,3) and (1,4) ======
\\ For pair (i, i+3), t_j = -t_i, so the Weil poly is T^4 + a2*T^2 + p^2.
\\ sf(D) = sf(t_i^2 - 4p) = sf(-3*u_i^2) = -3.
\\ h(Q(sqrt(-3))) = 1 (class number of Z[zeta_3] ring of integers = 1).

print("=== Biquadratic pairs (0,3) and (1,4) ===");

\\ Pair (0,3): traces t0 and -t0
a2_03 = 2*p - t0^2;
D_03 = a2_03^2 - 4*p^2;
sf_03 = core(D_03);
print("Pair (0,3): a2 = ", a2_03);
print("  D = a2^2 - 4p^2 = ", D_03);
print("  sf(D) = ", sf_03);
h_03 = quadclassunit(sf_03).no;
print("  h(Q(sqrt(", sf_03, "))) = ", h_03);
print("  [P]^2=1 trivially? ", if(h_03 == 1, "YES (h=1, every ideal principal)", "NO (h>1)"));
print("");

\\ Pair (1,4): traces ta and -ta
a2_14 = 2*p - ta^2;
D_14 = a2_14^2 - 4*p^2;
sf_14 = core(D_14);
print("Pair (1,4): a2 = ", a2_14);
print("  D = a2^2 - 4p^2 = ", D_14);
print("  sf(D) = ", sf_14);
h_14 = quadclassunit(sf_14).no;
print("  h(Q(sqrt(", sf_14, "))) = ", h_14);
print("  [P]^2=1 trivially? ", if(h_14 == 1, "YES (h=1)", "NO (h>1)"));
print("");

\\ Also check pair (2,5) for completeness (even though it is NOT glueable per Thread 3)
a2_25 = 2*p - tb^2;
D_25 = a2_25^2 - 4*p^2;
sf_25 = core(D_25);
print("Pair (2,5) [NOT glueable]: sf(D) = ", sf_25,
      ", h = ", quadclassunit(sf_25).no, "  [for reference]");
print("");

\\ ====== NON-BIQUADRATIC GLUEABLE PAIRS: (0,2), (0,5), (2,3) ======
\\ For these, t_i + t_j != 0, so the Weil poly is a general degree-4.
\\ The CM field is a quartic field; class number requires bnfinit.
\\ We show the Weil polynomial explicitly and compute its discriminant.

print("=== Non-biquadratic pairs: (0,2), (0,5), (2,3) ===");

\\ Pair (0,2): traces t0, tb
s_02 = t0 + tb; prod_02 = t0 * tb;
print("Pair (0,2): t0+tb = ", s_02, ", t0*tb = ", prod_02);
print("  Weil poly: T^4 - (", s_02, ")*T^3 + (", prod_02+2*p, ")*T^2 - (",
      s_02*p, ")*T + ", p^2);
\\ Disc of quartic P via poldisc (using variable T=x)
P_02 = x^4 - s_02*x^3 + (prod_02+2*p)*x^2 - s_02*p*x + p^2;
disc_02 = poldisc(P_02);
print("  poldisc(P_02) = ", disc_02);
sf_disc_02 = core(disc_02);
print("  core(poldisc) = ", sf_disc_02);
print("  (bnfinit on quartic with 256-bit coefficients: SKIPPED — too slow)");
print("");

\\ Pair (0,5): traces t0, -tb (note: t[5] = -tb)
s_05 = t0 + (-tb); prod_05 = t0 * (-tb);
print("Pair (0,5): t0+(-tb) = ", s_05, ", t0*(-tb) = ", prod_05);
print("  Weil poly: T^4 - (", s_05, ")*T^3 + (", prod_05+2*p, ")*T^2 - (",
      s_05*p, ")*T + ", p^2);
P_05 = x^4 - s_05*x^3 + (prod_05+2*p)*x^2 - s_05*p*x + p^2;
disc_05 = poldisc(P_05);
print("  poldisc(P_05) = ", disc_05);
print("");

\\ Pair (2,3): traces tb, -t0
s_23 = tb + (-t0); prod_23 = tb * (-t0);
print("Pair (2,3): tb+(-t0) = ", s_23, ", tb*(-t0) = ", prod_23);
print("  Weil poly: T^4 - (", s_23, ")*T^3 + (", prod_23+2*p, ")*T^2 - (",
      s_23*p, ")*T + ", p^2);
P_23 = x^4 - s_23*x^3 + (prod_23+2*p)*x^2 - s_23*p*x + p^2;
disc_23 = poldisc(P_23);
print("  poldisc(P_23) = ", disc_23);
print("");

\\ ====== FINAL SUMMARY ======
print("=== SUMMARY: Thread 20 ===");
print("For secp256k1 (j=0, CM by Z[zeta_3]):");
print("  Trace:    t0 = ", t0);
print("  Cornacchia: 4p = t0^2 + 3*u0^2 with u0 = ", u0);
print("");
print("BIQUADRATIC glueable pairs (t_j = -t_i => biquadratic Weil poly):");
print("  (0,3): sf(D) = ", sf_03, ", h(CM field) = ", h_03,
      " => [P]^2=1 trivially (h=1)");
print("  (1,4): sf(D) = ", sf_14, ", h(CM field) = ", h_14,
      " => [P]^2=1 trivially (h=1)");
print("");
print("NON-BIQUADRATIC glueable pairs (general quartic CM field):");
print("  (0,2): Weil poly has poldisc = ", disc_02);
print("         Class number of quartic CM field: UNKNOWN (bnfinit infeasible)");
print("  (0,5): poldisc = ", disc_05, " [reference]");
print("  (2,3): poldisc = ", disc_23, " [reference]");
print("");
print("KEY FINDING: The Thread 17 Proposition ([P]^2=1 in Cl(Q(sqrt(sf(D)))))");
print("is TRIVIALLY satisfied for all biquadratic glueable pairs because");
print("h(Q(sqrt(-3))) = 1 — every ideal is principal, so [P]^2=1 follows");
print("automatically without any non-trivial algebraic constraint.");
print("=> The Proposition gives no obstruction to cover attacks via these pairs.");
print("=> The obstruction for secp256k1 must come from other structure (Gaudry B5).");
print("");
print("OPEN: class number of quartic CM field for non-biquadratic pairs (0,2),(0,5),(2,3).");
print("These would require bnfinit on a quartic with 256-bit coefficients, or");
print("a descent to a smaller model (Weil restriction + isogeny argument).");
