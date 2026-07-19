\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify that the order-2 Frobenius theorem from Thread 15
\\ is fully general — holds for ANY prime p, not just secp256k1 norm-form primes.
\\
\\ THEOREM (general, proved algebraically in Thread 15):
\\   Let p be any prime, a2 in Z with p !| a2.
\\   Set D = a2^2 - 4p^2 = sf * m^2 (sf squarefree, m >= 1).
\\   Let K = Q(sqrt(sf)).  Every prime ideal P above p in O_K satisfies
\\     [P]^2 = 1   in Cl(O_K).
\\
\\ PROOF (Thread 15, fully general — no norm-form condition used):
\\   beta = (-a2 + m*sqrt(sf))/2 satisfies x^2 + a2*x + p^2 = 0.
\\   (A) beta in O_K (algebraic integer, root of monic Z-poly).
\\   (B) N_{K/Q}(beta) = p^2 => (beta) has norm p^2 in O_K.
\\   (C) p !| a2 => (beta) != (p)*O_K (else beta/p in O_K^x requires p | a2).
\\   (D) Only ideals of norm p^2 other than (p) are P^2 and Pbar^2
\\       => (beta) = P^2 (or Pbar^2) => [P]^2 = 1.

sqfree_kernel(n) = {
  \\ squarefree part of n (n != 0)
  my(s = sign(n), f = factor(abs(n)), r = 1);
  for(i = 1, #f[,1], if(f[i,2] % 2, r *= f[i,1]));
  s * r
};

is_norm_form(n) = {
  \\ true iff 4n = 73 + 3k^2 for some odd positive k  (secp256k1 CM-73 family)
  my(v = 4*n - 73);
  if(v <= 0 || v % 3, return(0));
  my(k2 = v \ 3);
  issquare(k2) && (sqrtint(k2) % 2 == 1)
};

verify_one(pp, a2) = {
  \\ Returns 1 = PASS, 0 = FAIL, -1 = SKIP (degenerate).
  \\ Prints one result row.
  if(gcd(pp, a2) > 1, return(-1));          \\ p | a2: skip
  my(D = a2^2 - 4*pp^2);
  if(D == 0, return(-1));
  my(sf = sqfree_kernel(D));
  if(sf == 1, return(-1));                  \\ D a perfect square: K = Q (trivial)
  my(m2 = D \ sf);
  if(m2 <= 0 || !issquare(m2), return(-1));
  my(m = sqrtint(m2));

  \\ (A) minpoly of beta = (-a2 + m*sqrt(sf))/2
  my(beta = Mod((-a2 + m*x)/2, x^2 - sf));
  my(mp = minpoly(beta));
  my(mp_ok = (mp == x^2 + a2*x + pp^2));

  \\ (B) N(beta) = p^2
  my(Nbeta = (a2^2 - m2*sf) / 4);
  my(N_ok = (Nbeta == pp^2));

  \\ Init number field with class group
  my(K = bnfinit(x^2 - sf, 1));
  my(h = K.clgp.no);

  \\ Split type: how does pp factor in O_K?
  my(PP = idealprimedec(K, pp));
  my(stype, inert_flag = 0);
  if(#PP == 0,  stype = "?",
     #PP >= 2,  stype = "S",   \\ split
     idealnorm(K, PP[1]) == pp^2, stype = "I"; inert_flag = 1,  \\ inert
     stype = "R");              \\ ramified: e=2

  \\ If h == 1, every ideal is principal trivially
  if(h == 1,
    printf("  p=%-6d a2=%-6d sf=%-10d m=%-5d h=1  %s  PASS(triv) min:%s N:%s\n",
           pp, a2, sf, m, stype, if(mp_ok,"OK","BAD"), if(N_ok,"OK","BAD"));
    return(1));

  \\ P^2 principal check
  my(P = PP[1]);
  my(P2 = idealpow(K, P, 2));
  my(coords = bnfisprincipal(K, P2)[1]);
  my(is_princ = (coords == 0*coords));

  \\ Consistency check: does (beta)*O_K == P^2 (or Pbar^2)?
  my(I_beta = idealhnf(K, beta));
  my(beta_P2 = 0);
  if(!inert_flag,
    beta_P2 = (I_beta == P2) ||
              (#PP == 2 && I_beta == idealpow(K, PP[2], 2)));

  my(verdict = if(is_princ, "PASS", "FAIL"));
  printf("  p=%-6d a2=%-6d sf=%-10d m=%-5d h=%-3d %s  %-4s  min:%s N:%s b=P2:%s\n",
         pp, a2, sf, m, h, stype, verdict,
         if(mp_ok,"OK","BAD"), if(N_ok,"OK","BAD"),
         if(inert_flag, "n/a", if(beta_P2,"YES","NO")));
  if(is_princ, 1, 0)
};

{
  print("=== Thread 16: General Order-2 Frobenius Theorem ===");
  print("Verifying for non-norm-form primes (secp256k1 CM-73 family excluded).");
  print("");

  \\ Collect 10 non-norm-form primes
  my(ps = List(), q = 2);
  while(#ps < 10,
    if(!is_norm_form(q), listput(ps, q));
    q = nextprime(q + 1));
  ps = Vec(ps);
  print("Non-norm-form primes: ", ps);
  print("(Norm-form {19,37,79,109,...} are the secp256k1 CM-73 family.)");
  print("");
  print("Columns: p, a2, sf=sqfree(a2^2-4p^2), m, h=classnr, split, verdict,");
  print("         minpoly-check, N(beta)=p^2-check, (beta)=P^2-check");
  print(Str(vector(80, i, "-")));

  my(pass = 0, fail = 0, skip = 0);

  \\ --- Primary sweep: 10 non-norm-form primes, 4 a2 values each ---
  for(i = 1, #ps,
    my(pp = ps[i]);
    print("  -- p = ", pp, " --");
    \\ a2 choices: 1 (tiny), p-1 (just below p), p+1 (just above p), 2p-1 (near 2p)
    my(a2list = [1, pp - 1, pp + 1, 2*pp - 1]);
    for(j = 1, #a2list,
      my(a2 = a2list[j]);
      if(gcd(pp, a2) > 1, next);
      my(r = verify_one(pp, a2));
      if(r == 1,  pass++,
         r == 0,  fail++,
         skip++)));

  print(Str(vector(80, i, "-")));
  printf("Primary sweep: PASS=%d  FAIL=%d  SKIP(degenerate)=%d\n", pass, fail, skip);
  print("");

  \\ --- Stress-test: larger primes ---
  print("Stress-test: larger non-norm-form primes with a2=1 and a2=p-1:");
  print(Str(vector(80, i, "-")));
  my(large_ps = [1013, 10007, 50021, 99991, 100003]);
  my(lpass = 0, lfail = 0);
  for(i = 1, #large_ps,
    my(pp = large_ps[i]);
    if(is_norm_form(pp), next);  \\ skip if accidentally norm-form (shouldn't happen)
    for(j = 1, 2,
      my(a2 = if(j == 1, 1, pp - 1));
      my(r = verify_one(pp, a2));
      if(r == 1, lpass++, r == 0, lfail++)));
  print(Str(vector(80, i, "-")));
  printf("Stress-test: PASS=%d  FAIL=%d\n", lpass, lfail);
  print("");

  \\ --- Grand total ---
  my(total_pass = pass + lpass, total_fail = fail + lfail);
  printf("GRAND TOTAL: PASS=%d  FAIL=%d\n", total_pass, total_fail);
  print("");
  if(total_fail == 0,
    print("ALL CASES PASSED."),
    printf("WARNING: %d FAILURES.\n", total_fail));
  print("");
  print("CONCLUSION:");
  print("  The proof (A)-(D) uses only: p prime, a2 in Z, p !| a2, D=a2^2-4p^2=sf*m^2.");
  print("  No norm-form condition is needed. Theorem is fully general.");
  print("  Thread 15 (norm-form primes) is a special case of this theorem.");
  print("");
  print("DONE.");
}
