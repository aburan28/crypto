\\ Find valid (p, b) toy case where Rosenhain formula CAN apply
\\ (cubic residue condition holds) and verify the formula is correct.
\\ Run: gp -q chlrs_valid_toy.gp

default(parisize, 64000000);

is_cubic_res(a, p) = lift(Mod(a, p)^((p-1)/3)) == 1;

print("=== Finding toy (p,b) where (-b) is a cubic residue mod p ===");
print("(Needed for F_p-rational 2-torsion, Rosenhain formula to apply)");
print("");

\\ Search over small primes p = 1 mod 3 and small b
found = 0;
forprime(p, 7, 200,
  if (p % 3 != 1, next());
  for (b = 1, 20,
    if (!is_cubic_res(-b, p), next());
    if (kronecker(4*b, p) == 1, next()); \\ skip j=1728 (b=0) or degenerate
    \\ also need a valid quadratic twist
    d_ns = -1;
    for (d = 2, p-1,
      if (kronecker(d, p) == -1, d_ns = d; break())
    );
    if (d_ns == -1, next());
    found = 1;
    printf("Found: p=%d, b=%d, d=%d\n", p, b, d_ns);
    printf("  (-b)^{(p-1)/3} = %d (should be 1)\n",
      lift(Mod(-b, p)^((p-1)/3)));
    rts = polrootsmod(x^3 + b, p);
    printf("  Roots of x^3+%d mod %d: [%d, %d, %d]\n",
      b, p, lift(rts[1]), lift(rts[2]), lift(rts[3]));
    alpha = lift(rts[1]);
    om = lift(Mod(rts[2], p) * Mod(rts[1], p)^(-1));
    printf("  omega = rts[2]/rts[1] = %d\n", om);
    printf("  omega^3 = %d (should be 1)\n", lift(Mod(om,p)^3));

    E1 = ellinit([0, b], p);
    b2 = lift(Mod(d_ns^3 * b, p));
    E2 = ellinit([0, b2], p);
    t1 = p + 1 - ellcard(E1);
    t2 = p + 1 - ellcard(E2);
    printf("  t(E1)=%d, t(E2)=%d (expected %d)\n", t1, t2, -t1);

    \\ Apply Rosenhain formula
    d_m = Mod(d_ns, p);
    w_m = Mod(om, p);
    l1 = (w_m^2-1)*(d_m-w_m) / (w_m*(w_m-1)*(d_m-1));
    l2 = (d_m*w_m-1)*(d_m-w_m) / (w_m*(d_m-1)^2);
    l3 = (d_m*w_m^2-1)*(d_m-w_m) / (w_m*(d_m*w_m-1)*(d_m-1));
    printf("  Rosenhain lambdas: l1=%d, l2=%d, l3=%d\n",
      lift(l1), lift(l2), lift(l3));

    h = Mod(1,p)*x*(x-l1)*(x-l2)*(x-1)*(x-l3);
    cp = hyperellcharpoly(h);
    n_jac = subst(cp, variable(cp), 1);
    target_n = (p+1-t1)*(p+1+t1);
    printf("  #Jac(C) = %d\n", n_jac);
    printf("  target (p+1-t1)(p+1+t1) = %d\n", target_n);
    printf("  match = %d\n", n_jac == target_n);
    if (n_jac == target_n,
      print("  RESULT: Rosenhain formula CORRECT for this (p,b)"),
      print("  RESULT: Formula still wrong - check Mobius ordering")
    );
    print("");
    break()
  );
  if (found, break())
);

print("=== Summary for secp256k1 ===");
p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
cr_secp = is_cubic_res(-7, p_secp);
printf("secp256k1: (-7) cubic residue mod p_secp? %d\n", cr_secp);
printf("  polrootsmod(x^3+7, p_secp): %d roots\n", #polrootsmod(x^3+7, p_secp));
if (cr_secp,
  print("  F_p Rosenhain formula applicable"),
  print("  F_p Rosenhain formula BLOCKED. Need F_{p^3} extension or Mestre approach.")
);
