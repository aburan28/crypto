\\ Diagnostic: why does the Rosenhain formula give the wrong Jac isogeny class?
\\
\\ Hypothesis: the formula requires (-b)^{1/3} ∈ F_p (all 2-torsion F_p-rational).
\\ This holds iff (-b) is a cubic residue mod p, i.e., (-b)^{(p-1)/3} ≡ 1 (mod p).
\\
\\ If this fails, the cross-ratios λ_i are computed formally but the 2-torsion
\\ points don't actually lie in F_p, so the Rosenhain model is for a different curve.
\\
\\ Run: gp -q chlrs_rosenhain_diagnostic.gp

default(parisize, 64000000);

is_cubic_residue(a, p) = {
  Mod(a, p)^((p-1)/3) == Mod(1, p)
};

check_cubic(b, p) = {
  my(neg_b);
  neg_b = (-b) % p;
  printf("  p = %d, b = %d: (-b) = %d, cubic residue? %d\n",
    p, b, neg_b, is_cubic_residue(neg_b, p));
  if (!is_cubic_residue(neg_b, p),
    printf("    → 2-torsion of E: y^2=x^3+%d NOT F_p-rational (Rosenhain inapplicable)\n", b),
    printf("    → 2-torsion F_p-rational: Rosenhain formula SHOULD work\n")
  );
};

print("=== Cubic residue check for Rosenhain applicability ===");
print("");

\\ Toy case
p1 = 1009; b1 = 11;
printf("Toy case (p=1009, b=11):\n");
check_cubic(b1, p1);
printf("\n");

\\ For p ≡ 1 mod 3, find b where cubic residue DOES hold
printf("Finding example where cubic residue HOLDS for p=1009:\n");
for(b_try = 1, 30,
  if(is_cubic_residue((-b_try) % 1009, 1009),
    printf("  b=%d: (-b) = %d is a cubic residue mod 1009\n", b_try, (-b_try) % 1009);
    break()
  )
);
print("");

\\ secp256k1
printf("secp256k1 (p ≡ 1 mod 3, b=7):\n");
p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
printf("  p mod 3 = %d (must be 1 for cube roots of unity to exist)\n", p_secp % 3);
check_cubic(7, p_secp);
print("");

\\ If cubic residue holds for secp256k1:
\\ compute the actual 2-torsion x-coordinates
printf("=== If cubic residue: compute 2-torsion x-coords for secp256k1 ===\n");
if(is_cubic_residue(p_secp - 7, p_secp),
  my(alpha, om);
  \\ alpha^3 = -7 mod p, solve using Tonelli-Shanks for cube roots
  \\ PARI: nthroot mod p
  alpha = lift(Mod(-7, p_secp)^((p_secp + 2) / 9));  \\ only works if p ≡ 7 mod 9
  printf("  p mod 9 = %d\n", p_secp % 9);
  \\ Use general cube root: since p ≡ 1 mod 3, use discrete log approach
  \\ PARI's polrootsmod for x^3 + 7 mod p
  my(rts);
  rts = polrootsmod(x^3 + 7, p_secp);
  if(#rts > 0,
    printf("  2-torsion x-coords (roots of x^3+7 mod p):\n");
    for(i=1, #rts, printf("    α_%d = %d\n", i, lift(rts[i])));
    om = lift(Mod(-1, p_secp) * Mod(rts[2], p_secp)^(-1) * Mod(rts[3], p_secp)^0);
    \\ cube root of unity: om = alpha2 / alpha1
    printf("  ω = α₂/α₁ = %d\n",
      lift(Mod(rts[2], p_secp) * Mod(rts[1], p_secp)^(-1)));
  ,
    printf("  x^3 + 7 has no roots mod p_secp (2-torsion not F_p-rational)\n");
  )
,
  printf("  Cubic residue fails: x^3 + 7 = 0 has no F_p solutions\n");
  printf("  → secp256k1 Rosenhain formula NOT directly applicable\n");
);
print("");

\\ CONCLUSION
print("=== Conclusion ===");
print("The Rosenhain formula in chlrs_igusa_formula.gp uses the 2-torsion");
print("x-coordinates as cross-ratios. If these are not in F_p, the formula");
print("produces a curve whose Jac is NOT isogenous to E1 × E2 over F_p.");
print("");
print("Fix: verify cubic residue before applying the formula. If it fails,");
print("options are:");
print("  (a) Work over F_{p^3} (Rosenhain model defined over extension)");
print("  (b) Use Cardona-Quer quartic model (doesn't need 2-torsion in F_p)");
print("  (c) Use Mestre's Step 2 reconstruction from Igusa invariants directly");
