\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify that Theorem (A)-(E) from Thread 15 is fully general —
\\ norm-form structure was never used; it holds for ANY (p, a2) with p prime, p∤a2.
\\
\\ THEOREM (generalized):
\\   Let p be a prime, a2 an integer with p∤a2, a2≠0.
\\   Set D = a2^2 - 4p^2, sf = squarefree part of D (D ≠ 0 since |a2|<2p or |a2|>2p).
\\   Write D = sf * m^2 with m ∈ Z_{>0}. Then:
\\     β = (-a2 + m*sqrt(sf)) / 2   lies in K = Q(sqrt(sf)),
\\     β is an algebraic integer (satisfies x^2 + a2*x + p^2 = 0),
\\     N_{K/Q}(β) = p^2,
\\     (β) ≠ (p) in O_K  [since p∤a2 forces this],
\\     Therefore (β) = P^2 or Pbar^2  =>  [P]^2 = 1 in Cl(K).
\\
\\ TESTS:
\\   Part A: 10 non-norm-form primes with arbitrary a2 (both D<0 and D>0).
\\   Part B: 5 non-norm-form primes with a2 derived from actual genus-2 curves
\\           (via hyperellcharpoly scan).
\\   Part C: Inertness impossibility — algebraic proof check.
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_order2.gp

default(parisize, 128000000);
default(timer, 0);

\\====================================================================
\\ Utilities
\\====================================================================

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_norm_form(p) = {
  my(v, k2);
  v = 4*p - 73;
  if(v <= 0 || v % 3 != 0, return(0));
  k2 = v \ 3;
  my(k = sqrtint(k2));
  (k*k == k2) && (k % 2 == 1)
};

\\====================================================================
\\ Core verifier: checks [P]^2=1 for given (p, a2)
\\====================================================================

verify_pair(p, a2, label) = {
  my(D, sf, m2, m, K, Pp, P, P2, P2_prin, norm_beta, result);

  if(a2 == 0 || a2 % p == 0,
    printf("  SKIP %-30s p=%d a2=%d: degenerate\n", label, p, a2);
    return(-1));

  D   = a2^2 - 4*p^2;
  if(D == 0,
    printf("  SKIP %-30s p=%d a2=%d: D=0 (a2=±2p)\n", label, p, a2);
    return(-1));

  sf  = sf_part(D);

  \\ Guard: sf=1 means D is a positive perfect square => K=Q (not quadratic).
  \\ bnfinit(x^2-1,1) would error since x^2-1 is reducible; skip these.
  if(sf == 1,
    printf("  SKIP %-30s p=%d a2=%d: D=%d is a perfect square (K=Q, not quadratic)\n",
           label, p, a2, D);
    return(-1));

  m2  = D / sf;
  if(denominator(m2) != 1 || m2 <= 0,
    printf("  ERROR %-29s D/sf=%Ps not pos int\n", label, m2);
    return(0));
  m = sqrtint(m2);
  if(m*m != m2,
    printf("  ERROR %-29s D/sf=%d not square\n", label, m2);
    return(0));

  \\ N(β) = (a2^2 - m^2*sf)/4 = (a2^2 - D)/4 = (4p^2)/4 = p^2
  norm_beta = (a2^2 - m2*sf) / 4;
  if(norm_beta != p^2,
    printf("  ERROR %-29s N(β)=%d ≠ p^2=%d\n", label, norm_beta, p^2);
    return(0));

  K    = bnfinit(x^2 - sf, 1);
  Pp   = idealprimedec(K, p);

  \\ For a quadratic field, idealprimedec always returns >=1 ideal.
  \\ #Pp==2 => split; #Pp==1,e=1 => inert; #Pp==1,e=2 => ramified.
  \\ Our theorem (p odd, p∤a2) guarantees split; anything else is a bug.
  if(#Pp != 2,
    my(e = Pp[1][3]);
    printf("  PARADOX %-28s p=%d in Q(sqrt(%d)): %s -- should not happen!\n",
           label, p, sf, if(e == 2, "RAMIFIED", "INERT"));
    return(0));

  P      = Pp[1];
  P2     = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);
  result = (P2_prin[1] == 0*P2_prin[1]);

  printf("  %-32s p=%-7d a2=%-8d sf=%-8d m=%-5d h=%-4d D%s  [P]²=1:%s\n",
    label, p, a2, sf, m, K.clgp.no,
    if(D < 0, "<0(imag)", ">0(real)"),
    if(result, "YES", "FAIL(!)"));
  if(result, 1, 0)
};

\\====================================================================
\\ Part A: 10 non-norm-form primes with arbitrary a2
\\====================================================================

print("Part A: Arbitrary (p,a2) pairs — 10 non-norm-form primes");
print("  Both D<0 (imaginary quadratic K) and D>0 (real quadratic K) tested.");
print();

{
  my(total=0, ok=0, r);
  my(cases = [
    \\ [p,    a2,    label]
    [23,   7,     "non-NF p=23 a2=7 (D<0)"],
    [29,  11,     "non-NF p=29 a2=11 (D<0)"],
    [41,  15,     "non-NF p=41 a2=15 (D<0)"],
    [43, -13,     "non-NF p=43 a2=-13 (D<0)"],
    [47,  19,     "non-NF p=47 a2=19 (D<0)"],
    [53,  21,     "non-NF p=53 a2=21 (D<0)"],
    [59, -23,     "non-NF p=59 a2=-23 (D<0)"],
    [61,  25,     "non-NF p=61 a2=25 (D<0)"],
    [67,  27,     "non-NF p=67 a2=27 (D<0)"],
    [71, -29,     "non-NF p=71 a2=-29 (D<0)"],
    \\ D > 0 cases: |a2| > 2p
    [23,  47,     "non-NF p=23 a2=47 (D>0)"],
    [29,  59,     "non-NF p=29 a2=59 (D>0)"],
    [41,  83,     "non-NF p=41 a2=83 (D>0)"],
    [13,  27,     "non-NF p=13 a2=27 (D>0)"],
    [17,  35,     "non-NF p=17 a2=35 (D>0)"]
  ]);

  for(i=1, #cases,
    my(p=cases[i][1], a2=cases[i][2], lbl=cases[i][3]);
    r = verify_pair(p, a2, lbl);
    if(r >= 0, total++; ok += r));
  printf("\nPart A: %d/%d passed\n\n", ok, total);
}

\\====================================================================
\\ Part B: a2 from actual genus-2 curves (hyperellcharpoly sweep)
\\         for non-norm-form primes
\\====================================================================

print("Part B: a2 from actual genus-2 curves for non-norm-form primes");
print("  Search: y^2 = x^6 + c*x^3 + d over F_p for various (c,d).");
print("  Keep only biquadratic Weil polys (T^3 coeff = 0).");
print();

\\ For each prime, scan a few curves looking for biquadratic Weil polys.
\\ hyperellcharpoly(f) returns T^4+a1*T^3+a2*T^2+a3*T+a4.
\\ Biquadratic: a1=0 (T^3 coeff vanishes).

find_biquadratic_curve(p) = {
  my(fld, f, cp, c1, c2, found, a2_val);
  fld = ffgen(p, 't);
  found = 0;
  \\ Try sextic y^2 = x^6 + b*x^3 + c (palindromic family)
  for(b = 1, p-1,
    for(c = 1, p-1,
      f = fld^0*x^6 + Mod(b,p)*x^3 + Mod(c,p)*fld^0;
      cp = hyperellcharpoly(f);
      c1 = polcoeff(cp, 3);   \\ T^3 coefficient
      if(c1 == 0,
        a2_val = polcoeff(cp, 2);
        found = 1;
        return([b, c, lift(a2_val)]))));
  if(!found, return([]));
};

{
  my(total=0, ok=0, r, res, a2_int, p);
  my(test_primes = [23, 29, 31, 41, 43]);

  for(i=1, #test_primes,
    p = test_primes[i];
    if(is_norm_form(p), next);
    res = find_biquadratic_curve(p);
    if(#res == 0,
      printf("  p=%d: no biquadratic curve found in sweep\n", p);
      next);
    my(b=res[1], c=res[2], a2=res[3]);
    \\ a2 is lifted from F_p, may be > p/2; centre it
    if(a2 > p\2, a2 = a2 - p);
    printf("  p=%d found: y^2=x^6+%d*x^3+%d gives a2=%d\n", p, b, c, a2);
    if(a2 == 0 || a2 % p == 0, printf("    (degenerate a2, skip)\n"); next);
    r = verify_pair(p, a2, Str("curve y²=x⁶+",b,"x³+",c));
    if(r >= 0, total++; ok += r));

  printf("\nPart B: %d/%d passed\n\n", ok, total);
}

\\====================================================================
\\ Part C: Mass sweep — 50 random (p, a2) pairs, p non-norm-form
\\====================================================================

print("Part C: Mass sweep — 40 (p,a2) pairs, p non-norm-form");
print();

{
  my(total=0, ok=0, r, primes, cnt, a2);
  primes = [23, 29, 31, 41, 43, 47, 53, 59, 61, 67,
            71, 73, 83, 89, 97, 101, 103, 107, 113, 127];

  \\ 40 pairs: 4 a2 values per prime (10 primes × 4; j=0 skipped)
  cnt = 0;
  for(i=1, 10,
    my(p = primes[i]);
    if(is_norm_form(p), next);
    forstep(j=-2, 2, 1,
      if(j == 0, next);
      a2 = p \ 3 + j * 7;
      if(a2 % p == 0 || a2 == 0, a2 = a2 + 1);
      r = verify_pair(p, a2, Str("sweep p=",p," a2=",a2));
      if(r >= 0, total++; ok += r);
      cnt++));
  printf("\nPart C: %d/%d passed\n\n", ok, total);
}

\\====================================================================
\\ Part D: Inertness impossibility check
\\         If p∤a2, then p CANNOT be inert in Q(sqrt(sf)).
\\====================================================================

print("Part D: Splitting corollary — p∤a2 (p odd) => p SPLITS in Q(sqrt(sf))");
print("  (not inert: #Pp=1,e=1,f=2; not ramified: #Pp=1,e=2,f=1)");
print();

\\ For quadratic K=Q(sqrt(sf)):
\\   idealprimedec(K,p) returns a length-2 vector => SPLIT
\\                      returns a length-1 vector, entry[3]=1,entry[4]=2 => INERT
\\                      returns a length-1 vector, entry[3]=2,entry[4]=1 => RAMIFIED
\\ Our corollary says: when p odd and p∤a2, we ALWAYS get SPLIT.

{
  my(sf, K, Pp, D, split_type, ok_count, total_count);
  ok_count = 0; total_count = 0;
  my(examples = [
    [3,   1],
    [5,   2],
    [7,   3],
    [11,  4],
    [13,  5],
    [17,  6],
    [23,  7],
    [29,  8]
  ]);
  for(i=1, #examples,
    my(p=examples[i][1], a2=examples[i][2]);
    D = a2^2 - 4*p^2;
    sf = sf_part(D);
    if(sf == 1, printf("  SKIP p=%d a2=%d: D perfect square\n", p, a2); next);
    K = bnfinit(x^2 - sf, 1);
    Pp = idealprimedec(K, p);
    if(#Pp == 2,
      split_type = "SPLIT";
      ok_count++;
    ,
      my(e = Pp[1][3], f = Pp[1][4]);
      if(e == 2, split_type = Str("RAMIFIED(e=",e,") PARADOX!"));
      if(f == 2, split_type = Str("INERT(f=",f,") PARADOX!"))
    );
    total_count++;
    printf("  p=%d a2=%d: sf=%d, #Pp=%d  => %s\n",
      p, a2, sf, #Pp, split_type));
  printf("\nPart D: %d/%d confirmed SPLIT (0 paradoxes expected)\n",
    ok_count, total_count);
}

print();
print("DONE. All checks above should show [P]^2=1:YES and Part D all SPLIT.");
