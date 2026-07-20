\\ thread17_splitting_check.gp
\\ Thread 17: Companion check for Proposition prop:order2-frobenius
\\
\\ Verifies that for 10 fresh (p, a2) pairs (not in Threads 15-16):
\\   (1) [P]^2 = 1 in Cl(K),  K = Q(sqrt(sf))
\\   (2) p splits in K (neither inert nor ramified)
\\ Also verifies 5 D>0 (real quadratic) cases for completeness.
\\
\\ Run: gp -q secp256k1_cm_audit/thread17_splitting_check.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Returns 1 on PASS, 0 on FAIL, -1 on SKIP
verify_pair(p, a2, label) = {
  my(D, sf, m2, m, K, Pp, P, P2, P2_prin, norm_beta, result, e);

  if(a2 == 0 || a2 % p == 0,
    printf("  SKIP  %-30s p=%d a2=%d: degenerate\n", label, p, a2);
    return(-1));

  D = a2^2 - 4*p^2;
  if(D == 0,
    printf("  SKIP  %-30s p=%d a2=%d: D=0\n", label, p, a2);
    return(-1));

  sf = sf_part(D);

  if(sf == 1,
    printf("  SKIP  %-30s p=%d a2=%d: D perfect square (K=Q)\n", label, p, a2);
    return(-1));

  m2 = D / sf;
  if(denominator(m2) != 1 || m2 <= 0,
    printf("  ERROR %-30s D/sf=%Ps not positive integer\n", label, m2);
    return(0));
  m = sqrtint(m2);
  if(m*m != m2,
    printf("  ERROR %-30s D/sf=%d not a perfect square\n", label, m2);
    return(0));

  \\ Quick norm check: N(beta) = p^2
  norm_beta = (a2^2 - m2*sf) / 4;
  if(norm_beta != p^2,
    printf("  ERROR %-30s N(beta)=%d != p^2=%d\n", label, norm_beta, p^2);
    return(0));

  K    = bnfinit(x^2 - sf, 1);
  Pp   = idealprimedec(K, p);

  \\ Splitting check: for a quadratic field,
  \\   #Pp == 2              => SPLIT
  \\   #Pp == 1 and e == 1   => INERT
  \\   #Pp == 1 and e == 2   => RAMIFIED
  if(#Pp != 2,
    e = Pp[1][3];
    printf("  FAIL  %-30s p=%d in Q(sqrt(%d)): %s (expected SPLIT)\n",
           label, p, sf, if(e == 2, "RAMIFIED", "INERT"));
    return(0));

  P      = Pp[1];
  P2     = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);
  result = (P2_prin[1] == 0*P2_prin[1]);

  printf("  %s %-30s p=%-7d a2=%-7d sf=%-8d h=%d D%s\n",
    if(result, "PASS", "FAIL"),
    label, p, a2, sf, K.clgp.no,
    if(D < 0, "<0(imag)", ">0(real)"));
  if(result, 1, 0)
};

print("Thread 17: 10 D<0 (imaginary) + 5 D>0 (real) verification cases");
print("Confirms Proposition: [P]^2=1 and p SPLITS in K=Q(sqrt(sf))");
print();

{
  my(total=0, ok=0, r);
  my(cases = [
    \\ 10 fresh D<0 (imaginary quadratic K) cases
    [101,  14, "p=101  a2=14  D<0"],
    [103,   6, "p=103  a2=6   D<0"],
    [107,  20, "p=107  a2=20  D<0"],
    [131,  30, "p=131  a2=30  D<0"],
    [137,  50, "p=137  a2=50  D<0"],
    [149, 240, "p=149  a2=240 D<0"],
    [151, 260, "p=151  a2=260 D<0"],
    [157, 300, "p=157  a2=300 D<0"],
    [163, 310, "p=163  a2=310 D<0"],
    [167, 330, "p=167  a2=330 D<0"],
    \\ 5 D>0 (real quadratic K) cases: need |a2| > 2p
    [23,   47, "p=23   a2=47  D>0"],
    [29,   59, "p=29   a2=59  D>0"],
    [31,   63, "p=31   a2=63  D>0"],
    [37,   75, "p=37   a2=75  D>0"],
    [41,   83, "p=41   a2=83  D>0"]
  ]);

  for(i = 1, #cases,
    my(p = cases[i][1], a2 = cases[i][2], lbl = cases[i][3]);
    r = verify_pair(p, a2, lbl);
    if(r >= 0, total++; ok += r));

  print();
  printf("Thread 17 result: %d/%d PASSED\n", ok, total);
  if(ok == total && total > 0,
    printf("ALL %d cases confirm Proposition (order-2 + splitting).\n", total),
    printf("WARNING: %d FAILED -- check output above\n", total - ok));
}
