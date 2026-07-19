\\ thread16_general_theorem.gp
\\
\\ Thread 16 — 2026-07-19
\\
\\ GOAL: Verify that the "universal order-2 Frobenius" theorem from Thread 15
\\ holds for NON-NORM-FORM primes, confirming it is a GENERAL theorem.
\\
\\ THEOREM (General):
\\   Let p be any prime, a2 any integer with gcd(a2,p)=1 and D = a2^2-4p^2 < 0.
\\   Let sf = squarefree_part(D), K = Q(sqrt(sf)).
\\   Then for every prime ideal P above p in O_K, [P]^2 = 1 in Cl(K).
\\
\\ PROOF (from Thread 15, applies verbatim):
\\   Let m^2 = D/sf and beta = (-a2 + m*sqrt(sf))/2.
\\   (A) beta satisfies x^2 + a2*x + p^2 = 0, so beta is in O_K.
\\   (B) N_{K/Q}(beta) = p^2.
\\   (C) Norm-p^2 ideals: P^2, Pbar^2, (p).  Case (beta)=(p) requires p|a2,
\\       contradicting gcd(a2,p)=1.
\\   (D) Therefore (beta) = P^2 (or Pbar^2), so [P]^2 = 1.  QED.
\\ The proof uses ONLY the polynomial x^2+a2*x+p^2 and gcd(a2,p)=1.
\\ No norm-form (4p = 73 + 3k^2) condition is needed.
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_theorem.gp

default(parisize, 128000000);
default(timer, 0);

\\ ============================================================
\\ Utilities
\\ ============================================================

sf_part(n) = if(n == 0, 0, core(n));

is_norm_form(p) = {
  my(v);
  if(p < 19, return(0));
  v = 4*p - 73;
  if(v < 0 || v % 3 != 0, return(0));
  issquare(v \ 3)
}

\\ ============================================================
\\ Main check: returns 1 if [P]^2=1, 0 if fails, -1 if skipped
\\ ============================================================

check_order2(p, a2) = {
  my(D, sf, m2, m, K, Pp, nP, P, P2, prin, h, nfac, ftype);

  D = a2^2 - 4*p^2;
  if(D >= 0, return(-1));
  if(gcd(a2, p) > 1, return(-1));

  sf = sf_part(D);
  m2 = D / sf;
  if(!issquare(m2), return(-1));
  m = sqrtint(m2);

  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);
  nP = #Pp;

  if(nP == 0,
    printf("  p=%d a2=%d: no prime above p (unexpected)\n", p, a2);
    return(-1));

  \\ Classify split behaviour
  if(nP == 1,
    ftype = if(Pp[1].f == 2, "inert", "ramif"),
    ftype = "split"
  );

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  prin = bnfisprincipal(K, P2);

  \\ prin[1] is the exponent vector; zero iff P^2 is principal
  my(ok = (prin[1] == 0*prin[1]));

  printf("  %s  p=%-5d a2=%-5d D=%-10d sf=%-8d m=%-5d h=%-4d %s\n",
    if(ok, "PASS", "FAIL"), p, a2, D, sf, m, h, ftype);
  if(ok, 1, 0)
}

\\ ============================================================
\\ Build list of 10 non-norm-form primes
\\ ============================================================

primes10 = [];
{
  my(pp = 2, cnt = 0);
  while(cnt < 10,
    pp = nextprime(pp + 1);
    if(!is_norm_form(pp),
      primes10 = concat(primes10, [pp]);
      cnt++
    )
  );
}

\\ ============================================================
\\ Run verification
\\ ============================================================

print("==========================================================");
print("Thread 16: Universal Order-2 Theorem (non-norm-form primes)");
print("Theorem: For any prime p, a2 with gcd(a2,p)=1 and D=a2^2-4p^2 < 0,");
print("         [P]^2 = 1 in Cl(Q(sqrt(squarefree(D)))) for P above p.");
print("==========================================================");
printf("Test primes (non-norm-form): %Ps\n\n", primes10);

{
  my(total = 0, passed = 0, n10 = #primes10, p, a2, r);

  for(i = 1, n10,
    p = primes10[i];
    printf("p = %d:\n", p);

    \\ Test a2 = 1, 2, 3, p-2  (all should have D < 0 since a2 << 2p)
    for(j = 1, 3,
      a2 = j;
      r = check_order2(p, a2);
      if(r >= 0, total++; if(r, passed++))
    );

    \\ Also test a2 near p (still D < 0 since p < 2p)
    if(p > 4,
      a2 = p - 2;
      r = check_order2(p, a2);
      if(r >= 0, total++; if(r, passed++))
    );

    print("")
  );

  printf("==========================================================\n");
  printf("Summary: %d PASS / %d tested\n", passed, total);
  if(passed == total && total > 0,
    print("ALL PASS — general theorem confirmed"),
    print("SOME FAIL")
  );
}

\\ ============================================================
\\ Bonus: explicit beta verification for p=5, a2=1
\\ ============================================================
print("");
print("--- Bonus: explicit (beta)=P^2 check for p=5, a2=1 ---");
{
  my(p5 = 5, a2 = 1, D, sf, m, K, Pp, P, P2, I_beta, prin, beta_elt);

  D  = a2^2 - 4*p5^2;    \\ = 1 - 100 = -99
  sf = sf_part(D);        \\ core(-99) = -11 (since -99 = -11*9)
  m  = sqrtint(D / sf);   \\ m = 3

  printf("p=%d a2=%d D=%d sf=%d m=%d\n", p5, a2, D, sf, m);
  printf("beta satisfies x^2 + a2*x + p^2 = x^2 + %d*x + %d = 0\n", a2, p5^2);

  K   = bnfinit(x^2 - sf, 1);
  printf("K = Q(sqrt(%d)), h = %d\n", sf, K.clgp.no);

  Pp = idealprimedec(K, p5);
  printf("Primes above 5 in O_K: %d (split=%s)\n", #Pp, if(#Pp==2,"yes","no"));

  P   = Pp[1];
  P2  = idealpow(K, P, 2);

  \\ beta as Polmod element of K (theta^2 = sf, so theta = sqrt(sf))
  \\ beta = (-a2 + m*theta)/2 = (-1 + 3*theta)/2
  beta_elt = Mod((-a2 + m*x) / 2, x^2 - sf);

  I_beta = idealhnf(K, beta_elt);
  prin   = bnfisprincipal(K, P2);

  printf("(beta) == P^2 (HNF match): %s\n",
    if(I_beta == P2 || I_beta == idealpow(K, Pp[#Pp], 2), "YES", "NO"));
  printf("[P]^2 principal: %s\n",
    if(prin[1] == 0*prin[1], "YES", "NO"));
}

print("");
print("DONE.");
quit
