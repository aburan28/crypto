\\ thread16_general_weil_poly.gp
\\
\\ Thread 16: Verify the universal order-2 Frobenius theorem holds for
\\ non-norm-form primes (beyond the secp256k1 CM family).
\\
\\ THEOREM (Thread 15, algebraic proof):
\\   Weil poly T^4 + a2*T^2 + p^2, with p prime and p not dividing a2.
\\   Let D = a2^2 - 4*p^2 = sf * m^2 (sf squarefree, D < 0).
\\   Let K = Q(sqrt(sf)), P a prime of O_K above p.
\\   THEN [P]^2 = 1 in Cl(K).
\\
\\ PROOF (recap): beta = (-a2 + m*sqrt(sf))/2 satisfies x^2 + a2*x + p^2 = 0,
\\   so beta is an algebraic integer in K with N(beta) = p^2.  Since p does not
\\   divide a2, we have (beta) != (p), hence (beta) = P^2 or P-bar^2.  QED.
\\
\\ Thread 15 proved this for 25 norm-form primes (k <= 199, 4p=73+3k^2).
\\ Thread 16 verifies for 10 non-norm-form primes, two sets of a2:
\\   Part A: a2 = 2p - t^2 (product surface E x twist(E), t=1..4)
\\   Part B: a2 in {p+1, 1, 2p-1, p-1, 3} (purely algebraic, no geometry)
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_weil_poly.gp

default(parisize, 128000000);
default(timer, 0);

\\============================================================
\\ is_norm_form(q): 1 if q = (73+3k^2)/4 for some odd k, else 0
\\============================================================

is_norm_form(q) = {
  my(k, val);
  forstep(k = 1, 201, 2,
    val = 73 + 3*k^2;
    if(val % 4 == 0 && val\4 == q, return(1))
  );
  0
};

\\============================================================
\\ check_order2(pp, aa): test [P]^2 = 1 in Cl(Q(sqrt(sf)))
\\ Returns [sf, m, h, ok] or [] (skip).
\\ ok=1 means P^2 is principal (order 2 verified).
\\============================================================

check_order2(pp, aa) = {
  my(D, sfd, sf, mm, bnf, fac, P1, P2, ep, hh);

  if(aa % pp == 0, return([]));   \\ pp | aa: not applicable

  D = aa^2 - 4*pp^2;
  if(D >= 0, return([]));         \\ need D < 0

  sfd = core(D, 1);
  sf  = sfd[1];                   \\ squarefree part (negative)
  mm  = sfd[2];                   \\ so D = sf * mm^2

  if(sf * mm^2 != D, return([]));

  bnf = bnfinit(x^2 - sf, 1);
  hh  = bnf.clgp.no;

  if(hh == 1, return([sf, mm, hh, 1]));  \\ h=1 => every ideal principal

  fac = idealprimedec(bnf, pp);
  if(#fac == 0, return([]));
  P1 = fac[1];

  P2 = idealpow(bnf, P1, 2);
  ep = bnfisprincipal(bnf, P2, 0);  \\ exponent in Cl(K); zero => principal

  return([sf, mm, hh, (norml2(ep) == 0)])
};

\\============================================================
\\ run_thread16(): main routine
\\============================================================

run_thread16() = {
  my(test_ps, pp, aa, r, sf, mm, hh, ok);
  my(all_A, n_A, all_B, n_B);

  \\ Collect 10 non-norm-form primes >= 97
  test_ps = [];
  forprime(q = 97, 500,
    if(!is_norm_form(q), test_ps = concat(test_ps, [q]));
    if(#test_ps >= 10, break)
  );
  printf("Non-norm-form test primes (%d): %Ps\n\n", #test_ps, test_ps);

  \\----------------------------------------------------------
  \\ Part A: a2 = 2p - t^2 (t=1..4): product surface E x twist(E)
  \\----------------------------------------------------------
  print("=== Part A: a2 = 2p - t^2 (product surface E x twist(E)) ===");
  print("p      t  a2        sf            m      h    P^2 princ?");
  print("-----------------------------------------------------------");

  all_A = 1; n_A = 0;
  for(i = 1, #test_ps,
    pp = test_ps[i];
    for(tt = 1, 4,
      aa = 2*pp - tt^2;
      r  = check_order2(pp, aa);
      if(#r == 0, next);
      sf = r[1]; mm = r[2]; hh = r[3]; ok = r[4];
      n_A++;
      if(!ok, all_A = 0);
      printf("%-6d %-2d %-9d %-13d %-6d %-4d %s\n",
             pp, tt, aa, sf, mm, hh, if(ok, "YES", "FAIL"))
    )
  );
  printf("\nPart A: %d checks, all YES: %s\n\n", n_A, if(all_A, "YES", "NO -- FAIL"));

  \\----------------------------------------------------------
  \\ Part B: arbitrary a2 (purely algebraic)
  \\----------------------------------------------------------
  print("=== Part B: arbitrary a2 (purely algebraic, no geometry) ===");
  print("p      a2        sf            m      h    P^2 princ?");
  print("----------------------------------------------------------");

  all_B = 1; n_B = 0;
  for(i = 1, #test_ps,
    pp = test_ps[i];
    my(a2list);
    a2list = [pp+1, 1, 2*pp-1, pp-1, 3];
    for(j = 1, #a2list,
      aa = a2list[j];
      r  = check_order2(pp, aa);
      if(#r == 0, next);
      sf = r[1]; mm = r[2]; hh = r[3]; ok = r[4];
      n_B++;
      if(!ok, all_B = 0);
      printf("%-6d %-9d %-13d %-6d %-4d %s\n",
             pp, aa, sf, mm, hh, if(ok, "YES", "FAIL"))
    )
  );
  printf("\nPart B: %d checks, all YES: %s\n\n", n_B, if(all_B, "YES", "NO -- FAIL"));

  \\----------------------------------------------------------
  \\ Summary
  \\----------------------------------------------------------
  print("=== Summary ===");
  printf("Total checks: %d (Part A: %d, Part B: %d)\n", n_A+n_B, n_A, n_B);
  printf("All P^2 principal: %s\n",
    if(all_A && all_B, "YES -- theorem universal", "NO -- COUNTEREXAMPLE"));
  if(all_A && all_B,
    print("CONCLUSION: Thread 15 order-2 theorem holds for all tested");
    print("non-norm-form primes, for both geometric and arbitrary a2.")
  );
};

run_thread16();
