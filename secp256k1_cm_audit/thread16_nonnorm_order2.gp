\\ thread16_nonnorm_order2.gp
\\
\\ Thread 16: Extend the "universal order-2 Frobenius" theorem to non-norm-form primes.
\\
\\ THEOREM (Thread 15 — fully general):
\\   Let C/F_p be a genus-2 curve with biquadratic Weil polynomial T^4 + a2*T^2 + p^2
\\   (equivalently: #C(F_p) = p+1, i.e. the T^3 and T coefficients vanish).
\\   Let D = a2^2 - 4p^2 = sf*m^2 with sf squarefree, m > 0.
\\   If p does NOT divide a2, then the prime P above p in K = Q(sqrt(sf))
\\   satisfies [P]^2 = 1 in Cl(K).
\\
\\ PROOF RECAP (A)-(E): uses ONLY biquadratic Weil poly + p not|a2.
\\   The secp256k1 norm-form condition (4p = 73 + 3k^2) is NOT used.
\\
\\ EXPERIMENT: Verify for 10 non-norm-form primes with biquadratic genus-2 curves.
\\
\\ Run: gp --stacksize 128000000 -q thread16_nonnorm_order2.gp

default(parisize, 128000000);
default(timer, 0);

\\ ---------------------------------------------------------------
\\ Utilities
\\ ---------------------------------------------------------------

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Is p a secp256k1 norm-form prime? i.e., 4p = 73 + 3k^2 for some k in Z.
is_norm_form(p) = {
  my(val = 4*p - 73, k2, k);
  if(val <= 0 || val % 3 != 0, return(0));
  k2 = val \ 3;
  k  = sqrtint(k2);
  k * k == k2
};

\\ Find genus-2 curve y^2 = x^5 + aa*x^3 + bb*x + cc over F_p
\\ with biquadratic Weil polynomial (T^3 coefficient = 0).
\\ Uses Mod(coeff, p) polynomial so hyperellcharpoly() (1-arg) works.
\\ Returns [aa, bb, cc, a2] or 0 if not found in `limit` attempts.
find_biquadratic(p, limit) = {
  my(f, g, chpoly, c1, a2);
  for(attempt = 1, limit,
    my(aa = random(p), bb = random(p), cc = random(p));
    f = Mod(1,p)*x^5 + Mod(aa,p)*x^3 + Mod(bb,p)*x + Mod(cc,p);
    \\ Smoothness check: gcd(f, f') should have degree 0
    g = gcd(f, deriv(f));
    if(poldegree(g) > 0, next);
    chpoly = hyperellcharpoly(f);  \\ returns Z-coefficient poly
    c1     = polcoeff(chpoly, 3);
    if(c1 == 0,
      a2 = polcoeff(chpoly, 2);
      return([aa, bb, cc, a2]));
  );
  0
};

\\ Verify [P]^2 = 1 for given (p, a2).  Prints result and returns verdict string.
verify_one(p, a2, aa, bb, cc) = {
  my(D, sf, m2, m, K, Pp, P, P2, P2_prin, h, nP, split_str, ord_str);

  D  = a2^2 - 4*p^2;
  sf = sf_part(D);

  if(D == 0 || sf == 0,
    printf("  p=%-5d a2=%-7d  D=0 (degenerate — skip)\n", p, a2);
    return("SKIP_D0"));

  m2 = D / sf;
  if(denominator(m2) != 1 || m2 <= 0,
    printf("  p=%-5d a2=%-7d  m2 not positive integer (%Ps) — skip\n", p, a2, m2);
    return("SKIP_m2"));

  m = sqrtint(m2);
  if(m * m != m2,
    printf("  p=%-5d a2=%-7d  D/sf=%d not a perfect square — skip\n", p, a2, m2);
    return("SKIP_sqrt"));

  \\ Condition (D): p does not divide a2
  if(a2 % p == 0,
    printf("  p=%-5d a2=%-7d  p|a2 — proof condition fails (skip)\n", p, a2);
    return("SKIP_p|a2"));

  K   = bnfinit(x^2 - sf, 1);
  h   = K.clgp.no;
  Pp  = idealprimedec(K, p);
  nP  = #Pp;

  if(nP == 0,
    printf("  p=%-5d a2=%-7d  sf=%-7d  p inert in Q(sqrt(sf)) — skip\n", p, a2, sf);
    return("SKIP_inert"));

  split_str = if(nP == 1, "ram", "spl");
  P         = Pp[1];
  P2        = idealpow(K, P, 2);
  P2_prin   = bnfisprincipal(K, P2);

  if(P2_prin[1] == 0 * P2_prin[1], ord_str = "YES", ord_str = "NO(!!)");

  printf(
    "p=%-5d  a2=%-7d  sf=%-8d  m=%-5d  h=%-3d  %s  curve x^5+%d*x^3+%d*x+%d  [P]^2=1: %s\n",
    p, a2, sf, m, h, split_str, aa, bb, cc, ord_str);

  ord_str
};

\\ ---------------------------------------------------------------
\\ Main sweep
\\ ---------------------------------------------------------------

print("Thread 16: Order-2 Theorem — Non-Norm-Form Prime Extension");
print("============================================================");
print("Proof (A)-(E): only biquadratic Weil poly + p not|a2 needed. No norm-form required.");
print();

{
  my(found = 0, target = 10, p = 50);
  while(found < target && p < 3000,
    p = nextprime(p);

    if(is_norm_form(p),
      printf("(p=%-4d norm-form — skip)\n", p);
      p++;
      next);

    my(lim = if(p < 200, 5000, if(p < 600, 3000, 1500)));
    my(res = find_biquadratic(p, lim));
    if(res == 0,
      printf("(p=%-4d: no biquadratic found in %d tries)\n", p, lim);
      p++;
      next);

    my(aa = res[1], bb = res[2], cc = res[3], a2 = res[4]);
    my(verdict = verify_one(p, a2, aa, bb, cc));

    if(verdict == "YES" || verdict == "NO(!!)",
      found++;
      if(verdict == "NO(!!)",
        print("  *** THEOREM VIOLATION — investigate! ***")));

    p++;
  );

  printf("\nVerified: %d of %d (target). All YES = theorem confirmed for non-norm-form primes.\n",
    found, target);
}

\\ ---------------------------------------------------------------
\\ Summary table of D, sf, m, h for the verified cases
\\ ---------------------------------------------------------------
\\ (already printed inline above)

\\ ---------------------------------------------------------------
\\ Cross-check: norm-form prime p=79 (k=9, known sf=-219, h=4)
\\ ---------------------------------------------------------------
print();
print("Cross-check: norm-form prime p=79 (expected sf=-219, h=4)");
{
  my(res = find_biquadratic(79, 10000));
  if(res != 0,
    my(a2 = res[4], D = a2^2 - 4*79^2, sf = sf_part(D));
    printf("  found a2=%d, sf=%d (expected sf=-219): %s\n",
      a2, sf, if(sf == -219, "MATCH", "different sf"));
    verify_one(79, a2, res[1], res[2], res[3]),
    print("  p=79: no biquadratic found in 10000 tries"));
}

print();
print("DONE.");
