\\ thread16_general_frobenius.gp
\\
\\ Thread 16: Test the algebraic Theorem from Thread 15 on NON-NORM-FORM primes.
\\
\\ THEOREM (generalised from Thread 15):
\\   For any ordinary E/F_p with trace a != 0, let A = E x E^t.
\\   Weil poly of A: T^4 + a2*T^2 + p^2  where  a2 = 2p - a^2.
\\   D = a2^2 - 4p^2 = a^2*(a^2-4p) < 0,  sf = squarefree(D),  m = sqrt(D/sf).
\\   Then [P]^2 = 1 in Cl(Q(sqrt(sf))), where P is the prime above p.
\\
\\   PROOF (no norm-form condition used):
\\   beta = (-a2 + m*sqrt(sf))/2 satisfies x^2+a2*x+p^2=0 (alg. int., norm p^2).
\\   Since |a|<p and a!=0 => p ∤ a => p ∤ a2 = 2p-a^2.
\\   => (beta) != (p) => (beta) = P^2 or Pbar^2 => [P]^2 = 1. QED.
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_frobenius.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_norm_form(p) = {
  my(val, k2);
  val = 4*p - 73;
  if(val <= 0 || val % 3 != 0, return(0));
  k2 = val / 3;
  return(issquare(k2));
};

get_trace(p, a4, a6) = {
  my(E, N);
  E = ellinit([a4, a6] * Mod(1,p));
  N = ellcard(E);
  p + 1 - N
};

\\ verify the theorem for one (p, a) pair; returns 1=pass, 0=fail, -1=skip
verify_general(p, a, a4, a6) = {
  my(a2, D, sf, m2, m, K, h, Pp, P, P2, P2_prin,
     beta_elt, I_beta, P2_iso_beta, p_div_a2, ord_str, norm_beta);

  if(a == 0, return(-1));
  if(a % p == 0, return(-1));

  a2 = 2*p - a^2;
  D  = a2^2 - 4*p^2;
  if(D >= 0, return(-1));

  sf = sf_part(D);
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1,
    printf("  ERROR D/sf noninteger: D=%d sf=%d\n", D, sf); return(0));
  m = sqrtint(m2);
  if(m*m != m2,
    printf("  ERROR D/sf=%d not perfect square\n", m2); return(0));

  p_div_a2 = (a2 % p == 0);
  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  Pp = idealprimedec(K, p);
  if(#Pp == 0,
    printf("  p=%d inert in Q(sqrt(%d)) -- UNEXPECTED\n", p, sf); return(-1));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);
  beta_elt = Mod((-a2 + m*x)/2, x^2 - sf);
  I_beta   = idealhnf(K, beta_elt);
  P2_iso_beta = (I_beta == P2);
  if(!P2_iso_beta && #Pp > 1,
    P2_iso_beta = (I_beta == idealpow(K, Pp[#Pp], 2)));

  if(P2_prin[1] == 0*P2_prin[1],
    ord_str = if(bnfisprincipal(K,P)[1] == 0*bnfisprincipal(K,P)[1], "1", "2"),
    ord_str = ">2");

  norm_beta = (a2^2 - m2*sf) / 4;

  printf("  p=%d a=%d (x^3+%d*x+%d): sf=%d m=%d h=%d nP=%d  N(b)=p^2:%s  p|a2:%s  (b)=P^2:%s  [P]^2=1:%s\n",
    p, a, a4, a6, sf, m, h, #Pp,
    if(norm_beta == p^2, "YES","NO"),
    if(p_div_a2, "YES(!)", "no"),
    if(P2_iso_beta, "YES","NO"),
    if(P2_prin[1] == 0*P2_prin[1], "YES","NO"));

  if(!P2_iso_beta || P2_prin[1] != 0*P2_prin[1] || norm_beta != p^2, return(0));
  1
};

\\ Test one prime with all curve params; returns [ok, tested]
run_prime(p, curve_params) = {
  my(cases_ok, cases_tested, ci, parms, a4, a6, a, r);
  cases_ok = 0; cases_tested = 0;
  for(ci = 1, #curve_params,
    parms = curve_params[ci];
    a4 = parms[1] % p;
    a6 = parms[2] % p;
    if((4*a4^3 + 27*a6^2) % p == 0, next);
    a = get_trace(p, a4, a6);
    r = verify_general(p, a, a4, a6);
    if(r == 1, cases_ok++; cases_tested++);
    if(r == 0, cases_tested++));
  [cases_ok, cases_tested]
};

\\ Main test driver
run_all() = {
  my(test_primes, curve_params, total_ok, total_tested, all_pass, pi, p, res);
  test_primes  = [101, 103, 107, 113, 127, 149, 163, 197, 211, 229];
  curve_params = [[0,1],[1,0],[1,1],[2,3],[3,7],[5,5],[7,11],[11,13],[0,6],[4,2]];
  total_ok = 0; total_tested = 0; all_pass = 1;

  for(pi = 1, #test_primes,
    p = test_primes[pi];
    if(!isprime(p), next);
    if(is_norm_form(p), printf("SKIP %d: IS norm-form\n", p); next);
    printf("\n--- p = %d (non-norm-form confirmed) ---\n", p);
    res = run_prime(p, curve_params);
    printf("  => %d / %d non-trivial cases passed\n", res[1], res[2]);
    if(res[2] > 0 && res[1] < res[2], all_pass = 0);
    total_ok     += res[1];
    total_tested += res[2]);

  printf("\n==========================================================================\n");
  printf("TOTAL: %d / %d cases passed\n", total_ok, total_tested);
  if(all_pass,
    print("All verified: YES -- Theorem holds for all non-norm-form primes tested"),
    print("All verified: NO -- check errors above"));
  all_pass
};

\\ Legendre check: p always splits in Q(sqrt(sf)) for ordinary E
run_legendre() = {
  my(test_primes, curve_params, fails, pi, p, ci, parms, a4, a6, a, a2, D, sf, leg);
  test_primes  = [101, 103, 107, 113, 127, 149, 163, 197, 211, 229];
  curve_params = [[0,1],[1,0],[1,1],[2,3],[3,7],[5,5],[7,11],[11,13],[0,6],[4,2]];
  fails = 0;

  for(pi = 1, #test_primes,
    p = test_primes[pi];
    if(!isprime(p), next);
    for(ci = 1, #curve_params,
      parms = curve_params[ci];
      a4 = parms[1] % p;
      a6 = parms[2] % p;
      if((4*a4^3 + 27*a6^2) % p == 0, next);
      a = get_trace(p, a4, a6);
      if(a == 0 || a % p == 0, next);
      a2 = 2*p - a^2;
      D  = a2^2 - 4*p^2;
      sf = sf_part(D);
      leg = kronecker(sf, p);
      if(leg != 1,
        printf("  FAIL p=%d a=%d sf=%d: Legendre=%d\n", p, a, sf, leg);
        fails++)));
  if(fails == 0,
    print("  All Legendre(sf,p)=+1: p ALWAYS splits in Q(sqrt(sf)) for ordinary E."),
    printf("  %d FAILURES.\n", fails));
};

\\ =====================================================================
\\ Execute
\\ =====================================================================

print("Thread 16: Generalised Theorem -- [P]^2=1 for arbitrary ordinary curves");
print("==========================================================================");
print("Testing 10 non-norm-form primes, multiple (a4,a6) per prime.");
print("Norm-form family: 4p = 73+3k^2 (secp256k1 CM relatives).");
print();

run_all();

print();
print("Bonus: verify Legendre(sf,p)=+1 for all cases (p splits in K)");
print("--------------------------------------------------------------");
run_legendre();

print();
print("DONE.");
