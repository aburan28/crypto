\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify [P]^2=1 in Cl(Q(sqrt(sf))) for non-norm-form primes.
\\
\\ THEOREM (Thread 15, generalized here to all biquadratic Weil polynomials):
\\ For any prime p, integer a2 with p not dividing a2, set D=a2^2-4p^2,
\\ sf=sf_part(D), m=sqrt(D/sf), beta=(-a2+m*sqrt(sf))/2 in K=Q(sqrt(sf)).
\\   (A) beta is algebraic integer: x^2+a2*x+p^2=0 (monic, Z coeff).
\\   (B) beta in K.   (C) N(beta)=p^2.
\\   (D) p not dividing a2 => (beta)!=( p) => (beta)=P^2 or Pbar^2 (split),
\\       or P principal (inert), or P^2=(p) (ramified).
\\   (E) In ALL cases [P]^2=1 in Cl(K).  QED.
\\
\\ Thread 15 proved this for 4p=73+3k^2 (secp256k1 norm-form family).
\\ Thread 16 verifies for 10 primes NOT in that family.
\\ We use a2=2p-t^2 (Frobenius trace relation for Jac ~ E x E^t).
\\
\\ Run: gp --stacksize 256000000 -q thread16_general_order2.gp

default(parisize, 256000000);
default(timer, 0);

\\ ---------------------------------------------------------------
\\ Utilities
\\ ---------------------------------------------------------------

sf_part(nn) = if(nn == 0, 0, sign(nn) * core(abs(nn)));

is_normform(pp) = {
  my(rr, kk2);
  rr = 4*pp - 73;
  if(rr <= 0 || rr % 3 != 0, return(0));
  kk2 = rr \ 3;
  issquare(kk2)
};

\\ ---------------------------------------------------------------
\\ verify_pair(pp, tt):
\\ Check [P]^2=1 for a2=2pp-tt^2, D=a2^2-4pp^2.
\\ Returns 1=YES, 0=skip/degenerate.
\\ ---------------------------------------------------------------

verify_pair(pp, tt) = {
  my(aa2, DD, sf, mm, KK, Pp, PP, PP2, PP2_prin, PP_prin, split_str, hh, is1, ord_str);

  aa2 = 2*pp - tt^2;
  DD  = aa2^2 - 4*pp^2;    \\ = -tt^2*(4pp-tt^2) <= 0

  if(DD == 0,
    printf("p=%-7d t=%-3d  D=0 (degenerate), skip\n", pp, tt);
    return(0));

  if(aa2 % pp == 0,
    printf("p=%-7d t=%-3d  p | a2, skip\n", pp, tt);
    return(0));

  sf = sf_part(DD);
  mm = sqrtint(abs(DD) \ core(abs(DD)));

  KK = bnfinit(x^2 - sf, 1);
  hh = KK.clgp.no;
  Pp = idealprimedec(KK, pp);

  if(#Pp >= 2, split_str = "SPLIT",
    if(Pp[1][3] == 2, split_str = "RAMIFIED", split_str = "INERT"));

  if(split_str != "SPLIT",
    printf("p=%-7d t=%-3d sf=%-9d mm=%-4d h=%-5d  %-10s  [P]^2=1: YES (trivially)\n",
           pp, tt, sf, mm, hh, split_str);
    return(1));

  PP       = Pp[1];
  PP2      = idealpow(KK, PP, 2);
  PP2_prin = bnfisprincipal(KK, PP2);
  is1      = (PP2_prin[1] == 0 * PP2_prin[1]);

  PP_prin  = bnfisprincipal(KK, PP);
  ord_str  = if(PP_prin[1] == 0*PP_prin[1], "ord([P])=1", "ord([P])=2");

  printf("p=%-7d t=%-3d a2=%-10d sf=%-9d mm=%-4d h=%-5d  SPLIT  %-12s  [P]^2=1:%s\n",
         pp, tt, aa2, sf, mm, hh, ord_str, if(is1, "YES", "NO(!)"));

  if(!is1, print("  *** THEOREM FAILURE ***"));
  is1
};

\\ ---------------------------------------------------------------
\\ Main: 10 non-norm-form primes, 3 traces each
\\ ---------------------------------------------------------------

print("Thread 16: Universal order-2 Frobenius for non-norm-form primes");
print("=================================================================");
print();

{
  my(PRIMES, TRACES, pp, tt, rr, total, passed);
  PRIMES = [1009, 1013, 1019, 1021, 1031, 9001, 9013, 9049, 99991, 100003];
  TRACES = [1, 3, 7];
  total = 0;
  passed = 0;

  for(ii = 1, #PRIMES,
    pp = PRIMES[ii];
    printf("\n--- p = %d  (norm-form: %s) ---\n",
           pp, if(is_normform(pp), "YES (unexpected!)", "NO"));
    if(!is_normform(pp),
      for(jj = 1, #TRACES,
        tt = TRACES[jj];
        rr = verify_pair(pp, tt);
        if(rr != 0, total = total + 1; if(rr == 1, passed = passed + 1))
      )
    )
  );

  printf("\n=================================================================\n");
  printf("Total cases: %d   Passed: %d   Failed: %d\n",
         total, passed, total - passed);
  printf("Universal order-2 theorem holds: %s\n",
         if(passed == total, "YES", "NO (check failures above)"));
  printf("=================================================================\n");
}

print();
print("DONE.");

\\ ---------------------------------------------------------------
\\ Phase 2: Search for SPLIT cases with ord([P])=2 (nontrivial).
\\ We want: prime p, integer trace tt, such that P above p in
\\ K=Q(sqrt(sf_part(D))) has [P] of order EXACTLY 2.
\\ This shows the theorem [P]^2=1 is tight (not always trivial).
\\ ---------------------------------------------------------------

print();
print("Phase 2: Searching for ord([P])=2 cases (nontrivial [P]^2=1)");
print("-------------------------------------------------------------");

search_ord2() = {
  my(pp, tt, aa2, DD, sf, mm, KK, Pp, PP, PP2, PP_prin, PP2_prin, ord2count);
  ord2count = 0;
  forprime(pp = 5, 3000,
    for(tt = 1, min(floor(2*sqrt(pp)) - 1, 15),
      aa2 = 2*pp - tt^2;
      DD  = aa2^2 - 4*pp^2;
      if(DD == 0 || aa2 % pp == 0, next());
      sf = sf_part(DD);
      if(sf > 0, next());    \\ only imaginary quadratic fields
      KK = bnfinit(x^2 - sf, 1);
      if(KK.clgp.no % 2 != 0, next());   \\ need even h for order-2 elements
      Pp = idealprimedec(KK, pp);
      if(#Pp < 2, next());   \\ only split
      PP = Pp[1];
      PP_prin  = bnfisprincipal(KK, PP);
      if(PP_prin[1] == 0*PP_prin[1], next());  \\ skip ord=1
      \\ P is non-principal: check P^2 is principal
      PP2 = idealpow(KK, PP, 2);
      PP2_prin = bnfisprincipal(KK, PP2);
      if(PP2_prin[1] != 0*PP2_prin[1], next());  \\ P^2 not principal (unexpected)
      mm = sqrtint(abs(DD) \ core(abs(DD)));
      printf("FOUND ord=2: p=%-5d t=%-3d a2=%-8d sf=%-8d mm=%-4d h=%-4d  [P]^2=1:YES  ord=2!\n",
             pp, tt, aa2, sf, mm, KK.clgp.no);
      ord2count++;
      if(ord2count >= 5, return(ord2count))   \\ stop after 5 examples
    )
  );
  if(ord2count == 0, print("No ord=2 cases found in search range."));
  ord2count
};

search_ord2();
print();
print("Phase 2 DONE.");
