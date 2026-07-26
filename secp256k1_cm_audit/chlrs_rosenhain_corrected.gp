\\ chlrs_rosenhain_corrected.gp
\\
\\ Corrected Rosenhain formula for Howe-glued genus-2 curves.
\\ 
\\ Problem found in chlrs_igusa_formula.gp: the algebraic formula for
\\ the lambdas was wrong -- it doesn't produce the correct Mobius-normalised
\\ branch points.  This script uses the ACTUAL roots of x^3+b mod p and
\\ directly computes the Rosenhain lambdas via the standard Mobius
\\ normalisation (r1→0, r2→∞, s1→1).
\\
\\ Run: gp -q chlrs_rosenhain_corrected.gp

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Corrected Rosenhain formula for j=0 Howe-glued genus-2 curves");
print("================================================================");
print("");

\\
\\ mobius_normalize(pts, r1, r2, s1, p):
\\   Apply T(x) = -((x-r1)/(x-r2)) * ((s1-r2)/(s1-r1))
\\   so T(r1)=0, T(r2)=inf, T(s1)=1.
\\   Returns the images of the remaining 3 points.
\\
mobius_image(x_val, r1, r2, s1, p) = {
  my(num, den, scale);
  num = Mod(x_val - r1, p);
  den = Mod(x_val - r2, p);
  scale = Mod(s1 - r2, p) / Mod(s1 - r1, p);
  lift(num / den * scale)
};

\\
\\ rosenhain_lambdas_corrected(p, b1, b2):
\\   Compute Rosenhain lambdas for y^2=(x^3+b1)(x^3+b2) over F_p.
\\   Requires x^3+b1 and x^3+b2 to have 3 F_p-rational roots each.
\\   Returns [] if not satisfied.
\\
\\ Mobius: send r1→0, r2→∞, s1→1 where
\\   {r1,r2,r3} = roots(x^3+b1), {s1,s2,s3} = roots(x^3+b2).
\\
rosenhain_lambdas_corrected(p, b1, b2) = {
  my(rts1, rts2, r1, r2, r3, s1, s2, s3, l1, l2, l3);
  rts1 = polrootsmod(Mod(1,p)*x^3 + b1, p);
  rts2 = polrootsmod(Mod(1,p)*x^3 + b2, p);
  if (#rts1 < 3, print("  ERROR: x^3+b1 not fully split over F_p"); return([]));
  if (#rts2 < 3, print("  ERROR: x^3+b2 not fully split over F_p"); return([]));
  r1 = lift(rts1[1]); r2 = lift(rts1[2]); r3 = lift(rts1[3]);
  s1 = lift(rts2[1]); s2 = lift(rts2[2]); s3 = lift(rts2[3]);
  \\ Mobius sends r1->0, r2->inf, s1->1
  l1 = mobius_image(r3, r1, r2, s1, p);
  l2 = mobius_image(s2, r1, r2, s1, p);
  l3 = mobius_image(s3, r1, r2, s1, p);
  [r1, r2, r3, s1, s2, s3, l1, l2, l3]
};

verify_rosenhain(p, lambdas, t1) = {
  my(l1, l2, l3, h, cp, expected_cp, n_jac, target_n, match);
  l1 = lambdas[7]; l2 = lambdas[8]; l3 = lambdas[9];
  h = Mod(1,p)*x*(x-1)*(x-l1)*(x-l2)*(x-l3);
  cp = hyperellcharpoly(h);
  expected_cp = x^4 + (2*p - t1^2)*x^2 + p^2;
  n_jac = subst(cp, variable(cp), 1);
  target_n = (p+1-t1)*(p+1+t1);
  [cp, expected_cp, n_jac, target_n, cp == expected_cp]
};

\\ ================================================================
\\ Part 1: p=13, b=1 — 2-torsion IS F_p-rational
\\ ================================================================
print("---- Part 1: p=13, b=1, d=2 ----");
{
  p = 13; b1 = 1; b2 = 8;  \\ b2 = d^3 * b1 = 2^3 = 8
  E1 = ellinit([0,b1], p);
  E2 = ellinit([0,b2], p);
  t1 = p + 1 - ellcard(E1);
  t2 = p + 1 - ellcard(E2);
  printf("  E1: y^2=x^3+%d, t=%d\n", b1, t1);
  printf("  E2: y^2=x^3+%d, t=%d (expected %d)\n", b2, t2, -t1);

  data = rosenhain_lambdas_corrected(p, b1, b2);
  if (#data == 0, print("  SKIP"); ,
    printf("  E1 roots: {%d, %d, %d}\n", data[1], data[2], data[3]);
    printf("  E2 roots: {%d, %d, %d}\n", data[4], data[5], data[6]);
    printf("  Rosenhain lambdas: l1=%d, l2=%d, l3=%d\n", data[7], data[8], data[9]);
    res = verify_rosenhain(p, data, t1);
    printf("  Char poly: %Ps\n", res[1]);
    printf("  Expected : %Ps\n", res[2]);
    printf("  #Jac=%d, target=%d, match=%d\n", res[3], res[4], res[5]);
  );
}
print("");

\\ Try all 6 permutations of (r1, r2, s1) to find a matching one
print("---- Part 1b: exhaustive Mobius search (all 9 choices of anchor triple) ----");
{
  p = 13; b1 = 1; b2 = 8;
  rts1 = polrootsmod(Mod(1,p)*x^3 + b1, p);
  rts2 = polrootsmod(Mod(1,p)*x^3 + b2, p);
  R = [lift(rts1[1]), lift(rts1[2]), lift(rts1[3])];
  S = [lift(rts2[1]), lift(rts2[2]), lift(rts2[3])];
  E1 = ellinit([0,b1],p); t1 = p+1-ellcard(E1);
  target_n = (p+1-t1)*(p+1+t1);
  for(ri=1, 3,
    for(rj=1, 3, if(ri==rj, next());
      for(si=1, 3,
        r1 = R[ri]; r2 = R[rj]; s1 = S[si];
        r3 = R[6-ri-rj]; \\ the third r
        s_others = [S[j] | j <- [1..3], j != si];
        l1 = mobius_image(r3, r1, r2, s1, p);
        l2 = mobius_image(s_others[1], r1, r2, s1, p);
        l3 = mobius_image(s_others[2], r1, r2, s1, p);
        h = Mod(1,p)*x*(x-1)*(x-l1)*(x-l2)*(x-l3);
        cp = hyperellcharpoly(h);
        n_jac = subst(cp, variable(cp), 1);
        if(n_jac == target_n,
          printf("  MATCH: (r1=%d,r2=%d,s1=%d) -> l1=%d,l2=%d,l3=%d, #Jac=%d\n",
            r1, r2, s1, l1, l2, l3, n_jac)
        );
      );
    );
  );
  print("  (done exhaustive search)");
}
print("");

\\ ================================================================
\\ Part 2: p=1009, b=11 — 2-torsion is NOT F_p-rational (irreducible)
\\ ================================================================
print("---- Part 2: p=1009, b=11 ----");
{
  p = 1009; b1 = 11;
  rts = polrootsmod(Mod(1,p)*x^3 + b1, p);
  printf("  Roots of x^3+11 mod 1009: %d (need 3 for F_p Rosenhain)\n", #rts);
  if (#rts < 3,
    print("  2-torsion NOT F_p-rational. F_p Rosenhain formula NOT applicable.");
    print("  BLOCKED: need F_{p^3} extension or Mestre/CM approach.");
  );
}
print("");

\\ ================================================================
\\ Part 3: Find a prime p where BOTH x^3+7 and x^3+b_twist split over F_p
\\ (Would allow the corrected formula for a secp256k1-like curve)
\\ ================================================================
print("---- Part 3: Search for p where x^3+7 splits AND the twist splits ----");
print("  (b1=7, so need -7 a cubic residue mod p, and p ≡ 1 mod 3)");
{
  found = 0;
  forprime(p = 7, 5000,
    if (p % 3 != 1, next());
    if (lift(Mod(-7, p)^((p-1)/3)) != 1, next());
    rts1 = polrootsmod(Mod(1,p)*x^3 + 7, p);
    if (#rts1 < 3, next());
    \\ Find non-square d and check x^3 + 7*d^3 splits
    d = 2;
    while(kronecker(d,p) != -1, d++);
    b2 = lift(Mod(7*d^3, p));
    rts2 = polrootsmod(Mod(1,p)*x^3 + b2, p);
    if (#rts2 < 3, next());
    E1 = ellinit([0,7],p);
    t1 = p+1-ellcard(E1);
    data = rosenhain_lambdas_corrected(p, 7, b2);
    res = verify_rosenhain(p, data, t1);
    printf("  p=%d, d=%d, b2=%d: l1=%d,l2=%d,l3=%d, match=%d\n",
      p, d, b2, data[7], data[8], data[9], res[5]);
    if (res[5],
      print("  *** FORMULA VERIFIED FOR b1=7 AT p=", p, " ***");
      found = 1;
    );
    if (found, break());
  );
  if (!found, print("  No matching prime found in search range for b1=7."));
}
print("");

\\ ================================================================
\\ Part 4: secp256k1 summary
\\ ================================================================
print("================================================================");
print("Part 4: secp256k1 summary");
print("================================================================");
{
  p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
  is_cubic_res_secp = (lift(Mod(-7, p_secp)^((p_secp-1)/3)) == 1);
  n_roots_secp = #polrootsmod(Mod(1,p_secp)*x^3 + 7, p_secp);
  printf("  (-7)^{(p-1)/3} mod p_secp ≡ 1? %d\n", is_cubic_res_secp);
  printf("  Roots of x^3+7 over F_p_secp: %d\n", n_roots_secp);
  print("");
  print("  STATUS: F_p Rosenhain BLOCKED — 2-torsion lives in F_{p^3}.");
  print("  CORRECT APPROACH for secp256k1:");
  print("  (A) F_{p^3} Rosenhain: compute lambdas over F_{p^3} extension,");
  print("      then verify the descent to F_p gives the Igusa invariants.");
  print("  (B) Mestre reconstruction: given (J2,J4,J6,J10) mod p, invert to");
  print("      get a model of C. See RESEARCH_MESTRE_HOWE.md.");
  print("  (C) CHLRS direct formula: compute (J2,J4,J6,J10) from CM data of E");
  print("      using Theorem 2.5 of Cardona-Howe-Lercier-Ritzenthaler-Streng.");
}
