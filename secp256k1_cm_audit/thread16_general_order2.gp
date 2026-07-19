\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify the universal order-2 Frobenius theorem BEYOND the norm-form family.
\\
\\ THEOREM (Thread 15, algebraic proof): For any prime p and any a2 with p ∤ a2 and
\\   D = a2^2 - 4p^2 = sf * m^2 (sf squarefree, m > 0), if p splits in K = Q(sqrt(sf)),
\\   then [P]^2 = 1 in Cl(K), where P is a prime of K above p.
\\
\\ PROOF: beta = (-a2 + m*sqrt(sf))/2 satisfies x^2 + a2*x + p^2 = 0 (alg. int., norm p^2).
\\   Since p ∤ a2 => (beta) ≠ (p), so (beta) = P^2 or Pbar^2 => [P]^2 = 1. QED.
\\
\\ KEY QUESTION (Thread 16): Does this hold for non-norm-form primes p?
\\   YES — the proof uses only the biquadratic form, NOT the norm-form 4p=73+3k^2.
\\   But it requires p SPLITS in Q(sqrt(sf)).
\\
\\ GEOMETRIC CONTENT:
\\   For ordinary E/F_p with trace t (|t| < 2*sqrt(p), p ∤ t), the product
\\   abelian surface A = E x E' (E' = quadratic twist) has biquadratic char. poly
\\     charpoly(A) = (T^2-tT+p)(T^2+tT+p) = T^4 + (2p-t^2)*T^2 + p^2.
\\   So a2 = 2p - t^2. Then:
\\     D = a2^2 - 4p^2 = (2p-t^2)^2 - 4p^2 = t^2*(t^2-4p) < 0 (since |t| < 2√p).
\\     sf = squarefree part of D.
\\     p ∤ a2 = 2p-t^2 iff p ∤ t^2 iff p ∤ t (E ordinary). ✓
\\     p splits in Q(sqrt(sf)): since D ≡ t^2*(t^2-4p) ≡ t^4 ≡ (t^2)^2 (mod p)
\\       and sf = D/m^2, if p ∤ m then sf ≡ (t^2/m)^2 (mod p) is a square => p splits. ✓
\\
\\ We test 10 non-norm-form primes, each with 3 different ordinary elliptic curves,
\\ producing 30 (p, a2) pairs. ALL should satisfy [P]^2 = 1.
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_order2.gp

default(parisize, 128000000);
default(timer, 0);

\\ ============================================================
\\ Utilities
\\ ============================================================

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_norm_form(p) = {
  my(disc, k2);
  disc = 4*p - 73;
  if(disc <= 0 || disc % 3 != 0, return(0));
  k2 = disc \ 3;
  issquare(k2)
};

\\ ============================================================
\\ Core verification: given prime p and a2, check [P]^2 = 1
\\ Returns: 1 = PASS, 0 = FAIL, -1 = SKIP (degenerate / p inert)
\\ ============================================================

verify_order2(p, a2, label) = {
  my(D, sf, m2, m, K, h, Pp, P, P2, P2_prin, beta_elt, norm_beta,
     I_beta, P2_iso_beta, result, ord_str);

  D = a2^2 - 4*p^2;

  if(D == 0,
    printf("  %-38s  a2=%-6d  SKIP: D=0 (degenerate)\n", label, a2);
    return(-1));

  sf = sf_part(D);
  if(sf == 0 || D % sf != 0,
    printf("  %-38s  a2=%-6d  SKIP: sf=0 or D not divisible by sf\n", label, a2);
    return(-1));

  m2 = D \ sf;
  if(m2 <= 0,
    printf("  %-38s  a2=%-6d  SKIP: D/sf=%d not positive\n", label, a2, m2);
    return(-1));

  m = sqrtint(m2);
  if(m*m != m2,
    printf("  %-38s  a2=%-6d  SKIP: D/sf=%d not perfect square\n", label, a2, m2);
    return(-1));

  if(a2 % p == 0,
    printf("  %-38s  a2=%-6d  SKIP: p|a2 (supersingular, theorem N/A)\n", label, a2);
    return(-1));

  \\ beta = (-a2 + m*sqrt(sf))/2, minpoly x^2+a2*x+p^2
  beta_elt = Mod((-a2 + m*x)/2, x^2 - sf);
  norm_beta = (a2^2 - m2*sf) / 4;   \\ = (a2^2 - D) / 4 = p^2

  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  Pp = idealprimedec(K, p);

  if(#Pp == 0,
    printf("  %-38s  a2=%-6d  sf=%-9d  SKIP: p inert in Q(sqrt(sf))\n",
      label, a2, sf);
    return(-1));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);

  I_beta = idealhnf(K, beta_elt);
  P2_iso_beta = (I_beta == P2 ||
    (#Pp > 1 && I_beta == idealpow(K, Pp[#Pp], 2)));

  if(P2_prin[1] == 0*P2_prin[1],
    ord_str = if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1],
      "1 (P princ)", "2"),
    ord_str = concat(">2 FAIL cyc=", Str(P2_prin[1]))
  );

  result = if(P2_prin[1] == 0*P2_prin[1], 1, 0);

  printf("  %-38s  a2=%-6d  sf=%-9d  m=%-4d  h=%-3d  N(b)=p^2:%s  (b)=P^2:%s  [P]^2=1:%s\n",
    label, a2, sf, m, h,
    if(norm_beta == p^2, "Y", "N"),
    if(P2_iso_beta, "Y", "N"),
    ord_str);

  result
};

\\ ============================================================
\\ Main
\\ ============================================================

print("Thread 16: Order-2 Frobenius theorem for non-norm-form primes");
print("==============================================================");
print("Claim: for ANY ordinary E/F_p (non-norm-form p),");
print("  A=E*E' has [P]^2=1 in Cl(Q(sqrt(sf(a2^2-4p^2)))), a2=2p-t_E^2.");
print("Expected: ALL cases PASS.");
print();

{
  my(test_primes = [23, 29, 41, 43, 53, 59, 61, 67, 71, 73]);
  my(test_curves = [[1,1], [2,3], [3,7]]);
  my(total = 0, passed = 0, skipped = 0);

  for(ii = 1, #test_primes,
    my(p = test_primes[ii]);
    printf("p = %d (norm-form: %s)\n", p, if(is_norm_form(p), "YES(!)", "no"));

    for(jj = 1, #test_curves,
      my(aa = test_curves[jj][1], bb = test_curves[jj][2]);

      \\ Skip singular curves mod p
      if((4*aa^3 + 27*bb^2) % p == 0, next());

      my(E = ellinit([aa, bb]));
      my(t = ellap(E, p));
      my(a2 = 2*p - t^2);
      my(label = Str("y^2=x^3+", aa, "x+", bb, " t_E=", t));

      my(r = verify_order2(p, a2, label));
      if(r == -1, skipped++,
        total++;
        if(r == 1, passed++));
    );
    print();
  );

  printf("Summary: %d/%d non-degenerate cases PASSED [P]^2=1 (skipped %d)\n",
    passed, total, skipped);
  if(passed == total && total > 0,
    printf("CONFIRMED: Order-2 theorem holds for ALL %d non-norm-form cases tested.\n", total),
    printf("WARNING: %d case(s) FAILED — investigate!\n", total - passed)
  );
}

\\ ============================================================
\\ Extension: Check larger non-norm-form primes to confirm scaling
\\ ============================================================

print();
print("--- Extension: 5 larger non-norm-form primes ---");

{
  my(large_primes = [1009, 1013, 1019, 1021, 1031]);
  my(total = 0, passed = 0);

  for(ii = 1, #large_primes,
    my(p = large_primes[ii]);
    if(is_norm_form(p), printf("p=%d NORM-FORM, skip\n", p); next());

    my(E = ellinit([1, 1]));
    my(t = ellap(E, p));
    my(a2 = 2*p - t^2);
    my(label = Str("p=", p, " t=", t));
    my(r = verify_order2(p, a2, label));
    if(r != -1, total++; if(r == 1, passed++));
  );

  printf("\nLarge-prime summary: %d/%d PASSED\n", passed, total);
}

print();
print("DONE.");
