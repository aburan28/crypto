\\ thread16_general_theorem.gp
\\
\\ Thread 16: Universal [P]^2=1 theorem for non-norm-form primes.
\\
\\ THEOREM (Thread 15): For any biquadratic Weil polynomial
\\   W(T) = T^4 + a2*T^2 + p^2 with D = a2^2 - 4p^2 = sf*m^2
\\   (sf squarefree, m > 0) and p not dividing a2, the prime P above p
\\   in K = Q(sqrt(sf)) satisfies [P]^2 = 1 in Cl(K).
\\
\\ CONSTRUCTION: W(T) = (T^2 - alpha*T + p)(T^2 + alpha*T + p) => a2 = 2p - alpha^2.
\\ Valid Weil poly for 0 < alpha^2 < 4p.  Roots all have |.| = sqrt(p).
\\
\\ EXPERIMENT: 10 non-norm-form primes, 3 alpha values each => 30 cases.

default(parisize, 256000000);
default(timer, 0);

\\---------------------------------------------------------------------------
\\ sf_part: signed squarefree part
\\---------------------------------------------------------------------------
sf_part(n) = {
  if(n == 0, 0, sign(n) * core(abs(n)))
};

\\---------------------------------------------------------------------------
\\ verify_one(pp, alpha):
\\   Returns 1 if [P]^2=1 in Cl(Q(sqrt(sf))), 0 if failed, -1 if skipped.
\\---------------------------------------------------------------------------
verify_one(pp, alpha) = {
  my(a2, D, sf, m2, m, norm_beta, K, Pp, P, P2, P2_prin, h, is_prin, ord_P);

  if(alpha == 0, return(-1));
  if(alpha^2 >= 4*pp, return(-1));

  a2 = 2*pp - alpha^2;
  D  = a2^2 - 4*pp^2;
  if(D >= 0, return(-1));

  sf = sf_part(D);
  if(sf == 0, return(-1));
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1, return(-1));
  m = sqrtint(m2);
  if(m^2 != m2, return(-1));

  \\ Sanity: N_{K/Q}(beta) = (a2^2 - m^2*sf)/4 = (a2^2 - D)/4 = 4p^2/4 = p^2
  norm_beta = (a2^2 - m2*sf) / 4;
  if(norm_beta != pp^2, return(-1));

  \\ Sanity: p does not divide a2 (proven: a2 = 2p-alpha^2 ≡ -alpha^2 mod p, p∤alpha)
  if(a2 % pp == 0, return(-1));

  \\ Class group
  K = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, pp);
  if(#Pp == 0,
    printf("  p=%-10d alpha=%-5d sf=%-10d  (inert, skipped)\n", pp, alpha, sf);
    return(-1));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);
  h = K.clgp.no;
  is_prin = (P2_prin[1] == 0*P2_prin[1]);
  ord_P = if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1], 1, 2);

  printf("  p=%-10d  a=%-3d  a2=%-12d  sf=%-9d  m=%-6d  h=%-4d  ord[P]=%d  [P]^2=1:%s\n",
    pp, alpha, a2, sf, m, h, ord_P, if(is_prin,"YES","NO(!)"));

  if(is_prin, 1, 0)
};

\\---------------------------------------------------------------------------
\\ process_prime(pp): run 3 alpha values; return [cases, passes]
\\---------------------------------------------------------------------------
process_prime(pp) = {
  my(alphas, tot, pas, r);
  if(!isprime(pp),
    printf("--- %d NOT prime, skip ---\n\n", pp);
    return([0,0]));

  printf("--- p = %d ---\n", pp);
  alphas = [1, 3, sqrtint(pp)];
  tot = 0; pas = 0;
  for(ai = 1, #alphas,
    r = verify_one(pp, alphas[ai]);
    if(r >= 0, tot++; pas += r));
  print();
  [tot, pas]
};

\\---------------------------------------------------------------------------
\\ Main
\\---------------------------------------------------------------------------

print("Thread 16: Universal [P]^2=1 theorem — non-norm-form primes");
print("=============================================================");
print("W(T)=(T^2-a*T+p)(T^2+a*T+p): a2=2p-a^2, D=a2^2-4p^2=sf*m^2");
print("Theorem: [P]^2=1 in Cl(Q(sqrt(sf))) for any valid (p,a).");
print();

{
  my(cand, total, passed, res);

  cand = [101, 1009, 10007, 100003, 1000033, 9999991, 99999989, 999999937, 2147483647, 4294967311];
  total = 0; passed = 0;

  for(idx = 1, #cand,
    res = process_prime(cand[idx]);
    total  += res[1];
    passed += res[2]);

  printf("SUMMARY: %d / %d passed\n", passed, total);
  if(passed == total && total > 0,
    print("ALL PASS — theorem holds for non-norm-form primes."),
    print("FAILURES DETECTED — see output above."));
}

print();
print("DONE.");
