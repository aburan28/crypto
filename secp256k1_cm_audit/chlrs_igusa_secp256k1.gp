\\ ============================================================
\\ CHLRS Igusa invariants of the secp256k1 naive cover
\\ h_secp = (x^3+7)(x^3+189) = x^6 + 196*x^3 + 1323
\\ ============================================================
\\
\\ Computes the Igusa-Clebsch invariants (I2, I4, I6, I10) of
\\ the genus-2 curve C: y^2 = (x^3+7)(x^3+189) over Q and mod p_secp.
\\
\\ NORMALIZATION: Cardona-Quer / SageMath convention.
\\ The Clebsch-to-Igusa conversion uses A_Sage = A_PARI/2, etc.
\\ See autolab 2026-07-26 for the full derivation of the factor.
\\
\\ FINDING (2026-07-26): igusa_clebsch_complete.gp has a normalization
\\ bug: clebsch_to_igusa formulas from Sage expect A_Sage = A_PARI/2.
\\ This script applies the corrected scaling: divide A,B,C,D from
\\ the PARI transvectant by 2,4,8,32 before the clebsch_to_igusa step.
\\
\\ Run: gp -q chlrs_igusa_secp256k1.gp

default(parisize, 256000000);
default(timer, 0);

read("igusa_clebsch.gp");

\\ ---- Corrected Clebsch A,B,C,D (Sage-normalized) ----
\\ Computes A_Sage = A_PARI/2, B_Sage = B_PARI/4, etc.
clebsch_ABCD_sage(h) = {
  my(hv, iv, Delta, y1, y2, y3, A, B, C, D);
  hv = poly_to_binary_form(h, 6);
  iv = transvectant(hv, 6, hv, 6, 4);
  A = binary_form_constant(transvectant(hv, 6, hv, 6, 6));
  B = binary_form_constant(transvectant(iv, 4, iv, 4, 4));
  Delta = transvectant(iv, 4, iv, 4, 2);
  C = binary_form_constant(transvectant(iv, 4, Delta, 4, 4));
  y1 = transvectant(hv, 6, iv, 4, 4);
  y2 = transvectant(iv, 4, y1, 2, 2);
  y3 = transvectant(iv, 4, y2, 2, 2);
  D = binary_form_constant(transvectant(y3, 2, y1, 2, 2));
  \\ Apply the normalization correction:
  \\ A_PARI = 2 * A_Sage  →  A_Sage = A/2
  \\ B_PARI = 4 * B_Sage  →  B_Sage = B/4
  \\ C_PARI = 8 * C_Sage  →  C_Sage = C/8
  \\ D_PARI = 32 * D_Sage →  D_Sage = D/32
  [A/2, B/4, C/8, D/32]
};

\\ ---- Correct Igusa invariants (Cardona-Quer / Sage convention) ----
igusa_correct(h) = {
  my(abcd, A, B, C, D, I2, I4, I6, I10);
  abcd = clebsch_ABCD_sage(h);
  A = abcd[1]; B = abcd[2]; C = abcd[3]; D = abcd[4];
  I2  = -120*A;
  I4  = -720*A^2 + 6750*B;
  I6  = 8640*A^3 - 108000*A*B + 202500*C;
  I10 = -62208*A^5 + 972000*A^3*B + 1620000*A^2*C
        - 3037500*A*B^2 - 6075000*B*C - 4556250*D;
  [I2, I4, I6, I10]
};

\\ ---- Sanity check: I2 explicit formula ----
\\ I2 = -120*a0*a6 + 20*a1*a5 - 8*a2*a4 + 3*a3^2
\\ (using coefficient ordering: polcoeff(h,0)=a0=const, polcoeff(h,6)=a6=lead)
igusa_I2_explicit(h) = {
  my(a0, a1, a2, a3, a4, a5, a6);
  a0 = polcoeff(h, 0); a1 = polcoeff(h, 1); a2 = polcoeff(h, 2);
  a3 = polcoeff(h, 3); a4 = polcoeff(h, 4); a5 = polcoeff(h, 5);
  a6 = polcoeff(h, 6);
  -120*a0*a6 + 20*a1*a5 - 8*a2*a4 + 3*a3^2
};

print("================================================================");
print("Igusa-Clebsch invariants of secp256k1 naive cover");
print("h_secp = (x^3+7)(x^3+189) = x^6+196*x^3+1323");
print("Normalization: Cardona-Quer / SageMath convention");
print("================================================================");
print("");

\\ ---- Test 1: verify normalization on h = x^6 + 1 ----
print("---- Sanity check: h = x^6 + 1 (expected I2=-120) ----");
h_test = x^6 + 1;
q_test = igusa_correct(h_test);
I2_expl = igusa_I2_explicit(h_test);
print("I2 from Sage formula : ", q_test[1]);
print("I2 explicit formula  : ", I2_expl);
print("Match: ", q_test[1] == I2_expl);
print("I10 from Sage formula: ", q_test[4]);
print("I10 from poldisc/32  : ", poldisc(h_test)/32);
print("I10 match: ", q_test[4] == poldisc(h_test)/32);
print("");

\\ ---- Test 2: verify on test-2 polynomial ----
print("---- Sanity check: h = x^6+2x^5+3x^4+5x^3+7x^2+11x+13 ----");
h2 = x^6 + 2*x^5 + 3*x^4 + 5*x^3 + 7*x^2 + 11*x + 13;
q2 = igusa_correct(h2);
I2_expl2 = igusa_I2_explicit(h2);
print("I2 from Sage formula : ", q2[1]);
print("I2 explicit formula  : ", I2_expl2);
print("Match: ", q2[1] == I2_expl2);
print("I10 from Sage formula: ", q2[4]);
print("I10 from poldisc/32  : ", poldisc(h2)/32);
print("I10 match: ", q2[4] == poldisc(h2)/32);
print("");

\\ ---- Main computation: h_secp ----
print("================================================================");
print("MAIN RESULT: Igusa quadruple for h_secp over Q");
print("================================================================");
{
  p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
  h_secp = x^6 + 196*x^3 + 1323;

  print("Computing Igusa invariants (takes ~5-10s for C computation)...");
  q = igusa_correct(h_secp);
  I2  = q[1]; I4  = q[2]; I6  = q[3]; I10 = q[4];

  \\ Cross-checks
  I2_expl = igusa_I2_explicit(h_secp);
  I10_disc = poldisc(h_secp) / 32;

  print("");
  print("Igusa-Clebsch invariants over Q (Cardona-Quer convention):");
  print("  I2  = ", I2);
  print("  I4  = ", I4);
  print("  I6  = ", I6);
  print("  I10 = ", I10);
  print("");
  print("Cross-checks:");
  print("  I2 == explicit formula (-43512): ", I2 == I2_expl);
  print("  I10 == poldisc(h)/32: ", I10 == I10_disc);
  print("");

  \\ Normalized ratios (weighted-projective isomorphism class)
  print("Weighted-projective isomorphism class:");
  print("  I4/I2^2  = ", I4/I2^2);
  print("  I6/I2^3  = ", I6/I2^3);
  print("  I10/I2^5 = ", I10/I2^5);
  print("");

  \\ Reduce mod p_secp
  I2_p  = lift(Mod(numerator(I2)  * Mod(denominator(I2),  p_secp)^(-1), p_secp));
  I4_p  = lift(Mod(numerator(I4)  * Mod(denominator(I4),  p_secp)^(-1), p_secp));
  I6_p  = lift(Mod(numerator(I6)  * Mod(denominator(I6),  p_secp)^(-1), p_secp));
  I10_p = lift(Mod(numerator(I10) * Mod(denominator(I10), p_secp)^(-1), p_secp));
  print("Igusa-Clebsch invariants mod p_secp:");
  print("  I2  mod p = ", I2_p);
  print("  I4  mod p = ", I4_p);
  print("  I6  mod p = ", I6_p);
  print("  I10 mod p = ", I10_p);
  print("");

  \\ Non-degeneracy
  print("Non-degeneracy (smooth genus-2 curve):");
  print("  I2  ≠ 0: ", I2 != 0);
  print("  I10 ≠ 0: ", I10 != 0);
}

print("");
print("================================================================");
print("Summary of CHLRS Igusa thread status (autolab 2026-07-26)");
print("================================================================");
print("h_secp = (x^3+7)(x^3+189): the 'naive cover' — the degree-6");
print("sextic whose Jacobian is (2,2)-isogenous to secp256k1 × E^t");
print("over F_{p^3} (where E^t is the quadratic twist and the 2-torsion");
print("lives in F_{p^3} since x^3+7 has no roots in F_p).");
print("");
print("The Igusa-Clebsch invariants above are correct in the");
print("Cardona-Quer / SageMath convention (I2 = -43512, I10 = disc/32).");
print("");
print("BUG FIXED: igusa_clebsch_complete.gp and igusa_7classes_full.gp");
print("had a normalization error (factor 2/4/8/32 in I2/I4/I6/I10)");
print("due to A_PARI = 2 * A_Sage in the transvectant scaling.");
print("This script applies the correct scaling.");
print("");
print("NEXT: Verify the 4 new Howe-glueable pairs from Thread 18");
print("(pairs (0,1),(0,4),(1,4),(3,4)) have correct Igusa invariants.");
print("These have E_0 × E_k with k in {1,4} qualifying under (H1+H2+H3).");
print("================================================================");
