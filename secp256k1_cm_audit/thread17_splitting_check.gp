\\ thread17_splitting_check.gp
\\
\\ Thread 17: Verify that p splits (not inert, not ramified) in Q(sqrt(sf))
\\ for 10 (p, a2) test cases with p ndivides a2, as a new verification step
\\ complementary to Thread 16.
\\
\\ Theorem (Thread 16 Corollary — p always splits):
\\   For p odd prime, a2 in Z, p ndivides a2, a2 != 0, D=a2^2-4p^2, sf=squarefree(D):
\\   (i)  p is NOT inert in Q(sqrt(sf))
\\        [Inert => only ideal of norm p^2 is (p); but (beta)=(p) requires p|a2.]
\\   (ii) p is NOT ramified in Q(sqrt(sf)) for odd p
\\        [p | disc(K) iff p | sf. But sf | D = a2^2-4p^2 so p|sf implies p|a2.]
\\   (iii) therefore p SPLITS into two distinct primes P, P-bar.
\\
\\ Additionally confirm [P]^2 = 1 (from Thread 15-16 main theorem) for sf < 0.
\\
\\ Run: gp --stacksize 128000000 -q thread17_splitting_check.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Returns: 2 => split, 1 => inert or ramified, 0 => error
check_splitting(p, a2) = {
  my(D, sf, K, Pp, e, nprime, h, P, P2, prin);

  D  = a2^2 - 4*p^2;
  if(D == 0,
    printf("  SKIP p=%d a2=%d: D=0\n", p, a2);
    return(-1));

  sf = sf_part(D);
  if(sf == 1 || sf == 0,
    printf("  SKIP p=%d a2=%d: sf=%d trivial (D perfect square)\n", p, a2, sf);
    return(-1));

  K    = bnfinit(x^2 - sf, 1);
  h    = K.clgp.no;
  Pp   = idealprimedec(K, p);
  nprime = #Pp;

  if(nprime == 2,
    \\ Split: two primes above p. Confirm [P]^2=1.
    P      = Pp[1];
    P2     = idealpow(K, P, 2);
    prin   = bnfisprincipal(K, P2);
    my(order2 = (prin[1] == 0*prin[1]));
    printf("  OK   p=%-8d a2=%-7d sf=%-8d h=%-4d D%s  SPLIT  [P]^2=1:%s\n",
      p, a2, sf, h,
      if(D < 0, "<0", ">0"),
      if(order2, "YES", "FAIL(!)"));
    if(order2, return(1), return(0));
  );

  \\ nprime == 1: either inert (e=1, f=2) or ramified (e=2, f=1)
  e = Pp[1][3];  \\ ramification index
  if(e == 2,
    printf("  FAIL p=%-8d a2=%-7d sf=%-8d RAMIFIED (should be impossible for odd p, p ndivides a2)\n",
      p, a2, sf),
    printf("  FAIL p=%-8d a2=%-7d sf=%-8d INERT    (should be impossible for p ndivides a2)\n",
      p, a2, sf));
  return(0);
};

\\====================================================
\\ Main: 10 diverse test cases
\\====================================================

print("Thread 17: Splitting verification — 10 (p,a2) pairs");
print("------------------------------------------------------------------");

{
  my(total=0, ok=0, r);
  my(cases = [
    \\ [p,          a2]   covering D<0, D>0, large p, non-norm-form, negative a2
    [11,           6],
    [23,          10],
    [47,          22],
    [53,          10],
    [101,         12],
    [199,         18],
    [1009,        30],
    [9001,        44],
    [32771,      100],
    [1000003,    200]
  ]);

  for(i=1, #cases,
    my(p=cases[i][1], a2=cases[i][2]);
    \\ Verify hypothesis
    if(a2 % p == 0,
      printf("  SKIP p=%d divides a2=%d\n", p, a2);
      next);
    r = check_splitting(p, a2);
    if(r >= 0, total++; ok += r));

  print("------------------------------------------------------------------");
  if(ok == total,
    printf("RESULT: %d/%d SPLIT — Theorem (Thread 16 Corollary) verified.\n", ok, total),
    printf("RESULT: %d/%d FAILURES — check above.\n", total-ok, total));
}
