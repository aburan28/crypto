\\ CHLRS Rosenhain: cubic-residue pre-condition check
\\
\\ The Rosenhain formula for the Howe-glued genus-2 curve C with
\\ Jac(C) (2,2)-isogenous to (E × E^t) requires the 2-torsion of E
\\ to be F_p-rational.  For E: y^2 = x^3 + b (j=0, p ≡ 1 mod 3),
\\ the 2-torsion x-coords satisfy x^3 = -b, so this holds iff
\\ (-b)^{(p-1)/3} ≡ 1 (mod p).
\\
\\ Run: gp -q chlrs_cubic_residue_check.gp

default(parisize, 64000000);

is_cubic_res(a, p) = lift(Mod(a, p)^((p-1)/3)) == 1;

p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;

print("=== Cubic-residue check for Rosenhain formula applicability ===");
print("");

\\ Toy case from chlrs_igusa_formula.gp
p1 = 1009; b1 = 11;
printf("Toy (p=%d, b=%d): (-b)^{(p-1)/3} mod p = %d  →  cubic res? %d\n",
  p1, b1, lift(Mod(-b1, p1)^((p1-1)/3)), is_cubic_res(-b1, p1));
print("  → Rosenhain formula inapplicable for (1009, b=11) [confirms match=0]");
print("");

\\ secp256k1
printf("secp256k1 (p, b=7): cubic res? %d\n", is_cubic_res(-7, p_secp));
printf("  polrootsmod(x^3+7, p_secp): %d roots\n", #polrootsmod(x^3+7, p_secp));
print("  → 2-torsion NOT F_p-rational; Rosenhain over F_p inapplicable");
print("");

\\ Find a toy case where cubic residue HOLDS
print("=== Toy case where cubic residue HOLDS (Rosenhain sanity check) ===");
p2 = 1009;
found_b = 0;
for (b_try = 1, 200,
  if (is_cubic_res(-b_try, p2) && kronecker(b_try^2*4, p2) == 1,
    found_b = b_try;
    break
  )
);
if (found_b == 0,
  print("No suitable b found in [1,200] for p=1009");
,
  printf("Found b=%d: (-b)^{(p-1)/3} = %d (cubic residue ✓)\n",
    found_b, lift(Mod(-found_b, p2)^((p2-1)/3)));

  \\ Verify: roots of x^3 = -b mod p2
  rts = polrootsmod(x^3 + found_b, p2);
  printf("  Roots of x^3+%d mod %d: %d roots\n", found_b, p2, #rts);
  if (#rts == 3,
    alpha = lift(rts[1]);
    alpha2 = lift(rts[2]);
    alpha3 = lift(rts[3]);
    printf("  2-torsion x-coords: %d, %d, %d\n", alpha, alpha2, alpha3);
    om = lift(Mod(alpha2, p2) * Mod(alpha, p2)^(-1));
    printf("  omega = alpha2/alpha1 = %d (should be cube root of unity)\n", om);
    printf("  omega^3 mod p = %d (should be 1)\n", lift(Mod(om, p2)^3));
    printf("  omega^2+omega+1 mod p = %d (should be 0)\n",
      lift(Mod(om^2+om+1, p2)));
    print("");

    \\ Find a non-square d
    d_ns = 0;
    for (d_try = 2, 100,
      if (kronecker(d_try, p2) == -1, d_ns = d_try; break)
    );
    printf("  Non-square d = %d\n", d_ns);
    b2_val = lift(Mod(d_ns^3 * found_b, p2));
    printf("  E2 coefficient: b2 = d^3 * b = %d\n", b2_val);
    E1 = ellinit([0, found_b], p2);
    E2 = ellinit([0, b2_val], p2);
    t1 = p2 + 1 - ellcard(E1);
    t2 = p2 + 1 - ellcard(E2);
    printf("  Trace E1 = %d, trace E2 = %d (expected %d)\n", t1, t2, -t1);
    print("");

    \\ Apply Rosenhain formula
    d = Mod(d_ns, p2);
    w = Mod(om, p2);
    l1 = (w^2-1)*(d-w) / (w*(w-1)*(d-1));
    l2 = (d*w-1)*(d-w) / (w*(d-1)^2);
    l3 = (d*w^2-1)*(d-w) / (w*(d*w-1)*(d-1));
    printf("  Rosenhain lambdas: l1=%d, l2=%d, l3=%d\n",
      lift(l1), lift(l2), lift(l3));

    \\ Construct curve y^2 = x(x-1)(x-l1)(x-l2)(x-l3)
    h_poly = Mod(1,p2)*x*(x-1)*(x-l1)*(x-l2)*(x-l3);
    cp = hyperellcharpoly(h_poly);
    n_jac = subst(cp, variable(cp), 1);
    target_n = (p2+1-t1)*(p2+1+t1);
    printf("  Frobenius char poly: %Ps\n", cp);
    printf("  #Jac(C) = %d, target (p+1-t)(p+1+t) = %d, match = %d\n",
      n_jac, target_n, n_jac == target_n);
    if (n_jac == target_n,
      print("  ✓ Rosenhain formula CORRECT for this (p,b) pair"),
      print("  ✗ Rosenhain formula still wrong — deeper issue")
    );
  ,
    print("  Only ", #rts, " roots found (expected 3) — bug in cubic_res check")
  )
);
print("");

\\ Summary
print("=== Summary ===");
print("The Rosenhain formula requires (-b) to be a cubic residue mod p.");
print("For secp256k1 (b=7, p=p_secp): FAILS. x^3+7 has no F_p roots.");
print("For toy (b=11, p=1009): FAILS. Explains match=0 in chlrs_igusa_formula.gp.");
print("");
print("Next steps for CHLRS Igusa (Priority 2):");
print("  Option A: Work over F_{p^3} — compute Rosenhain lambdas in F_{p^3}");
print("            then check if the resulting curve descends to F_p.");
print("  Option B: Use Mestre's reconstruction from Igusa invariants of");
print("            (E × E^t) / Gamma_alpha directly (avoids 2-torsion rationality).");
print("  Option C: Search for b' s.t. j(E')=0, (-b') cubic res, E' F_p-isogenous");
print("            to secp256k1 (twist chain: may require base field extension).");
print("  STATUS: BLOCKED on 2-torsion rationality for F_p Rosenhain form.");
