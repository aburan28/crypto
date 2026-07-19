\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify that [P]^2 = 1 in Cl(Q(sqrt(sf))) holds for primes p
\\ NOT in the norm-form family 4p = 73 + 3k^2.
\\
\\ THEOREM (Thread 15, algebraic proof A-E):
\\   For any prime p and integer a2 with D = a2^2 - 4p^2 < 0,
\\   squarefree part sf = core(D) and m = sqrt(D/sf),
\\   K = Q(sqrt(sf)), P a prime ideal above p in O_K.
\\   If p does not divide a2, then [P]^2 = 1 in Cl(K).
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_order2.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(nn) = {
  if(nn == 0, return(0));
  sign(nn) * core(abs(nn));
};

is_norm_form(p) = {
  my(val = 4*p - 73);
  if(val < 0 || val % 3 != 0, return(0));
  issquare(val \ 3);
};

\\ Returns [p, a2, sf, m, nprimes_above_p, class_number, order2_confirmed]
verify_general(p, a2) = {
  my(D, sf, m2, m, K, Pp, nP, P, P2, cl, is_zero, h);
  D = a2^2 - 4*p^2;
  if(D >= 0, error("Need D<0"));
  if(a2 % p == 0, error("p|a2"));
  sf = sf_part(D);
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1, error("D/sf not pos int"));
  m = sqrtint(m2);
  if(m*m != m2, error("D/sf not perfect square"));
  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  Pp = idealprimedec(K, p);
  nP = #Pp;
  if(nP == 0, return([p, a2, sf, m, 0, h, 1]));
  P = Pp[1];
  P2 = idealpow(K, P, 2);
  cl = bnfisprincipal(K, P2)[1];
  is_zero = (cl == 0*cl);
  [p, a2, sf, m, nP, h, is_zero];
};

run_fixed_cases() = {
  my(cases, res, all_pass, nf, p, a2);
  \\ 20 diverse (p, a2) pairs, p NOT norm-form
  cases = [[7,3],[7,5],[11,7],[13,1],[13,-9],[17,11],[23,13],[31,1],[41,17],[47,23],[53,29],[61,1],[71,43],[97,61],[101,1],[101,99],[127,83],[151,71],[199,113],[251,199]];
  print("=== Thread 16: [P]^2=1 theorem for non-norm-form primes ===");
  print("");
  print("p      a2     sf           m     #P_above_p  h    [P]^2=1?  norm_form?");
  print("----------------------------------------------------------------------");
  all_pass = 1;
  for(i = 1, #cases,
    p = cases[i][1];
    a2 = cases[i][2];
    nf = is_norm_form(p);
    res = verify_general(p, a2);
    printf("%-6d %-6d %-12d %-5d %-11d %-4d %-9s %s\n",
      res[1], res[2], res[3], res[4], res[5], res[6],
      if(res[7], "YES", "NO!"),
      if(nf, "YES(skip)", "no"));
    if(!res[7], all_pass = 0);
  );
  print("");
  if(all_pass,
    print("ALL 20 PASSED"),
    print("FAILURES DETECTED"));
  all_pass;
};

run_stress_test() = {
  my(p, a2, res, cnt, ok, fail);
  print("");
  print("=== Stress test: 30 non-norm-form primes p in [1000,5000], a2=p\\3 ===");
  print("p        a2       sf              m      #P   h    [P]^2=1?");
  print("-------------------------------------------------------------");
  cnt = 0; ok = 0; fail = 0;
  p = 999;
  while(cnt < 30,
    p = nextprime(p + 1);
    if(p > 5000, break);
    if(is_norm_form(p), next);
    a2 = p \ 3;
    if(a2 % p == 0, a2 = a2 + 1);
    if(a2^2 - 4*p^2 >= 0, next);
    res = verify_general(p, a2);
    printf("%-8d %-8d %-15d %-6d %-4d %-4d %s\n",
      res[1], res[2], res[3], res[4], res[5], res[6],
      if(res[7], "YES", "NO!"));
    if(res[7], ok++, fail++);
    cnt++;
  );
  print("");
  printf("Stress: %d passed, %d failed (of %d)\n", ok, fail, cnt);
  if(fail == 0,
    print("THEOREM CONFIRMED GENERAL for all tested non-norm-form primes."),
    print("FAILURES: review cases above."));
  [ok, fail, cnt];
};

\\ Main
{
  my(r1, r2);
  r1 = run_fixed_cases();
  r2 = run_stress_test();
  print("");
  if(r1 && r2[2] == 0,
    print("OVERALL: ALL CASES PASS. Theorem (A-E) is general."),
    print("OVERALL: SOME FAILURES."));
}
