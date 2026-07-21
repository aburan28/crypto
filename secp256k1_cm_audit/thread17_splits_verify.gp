\\ thread17_splits_verify.gp
\\
\\ Thread 17: Verify the corollary that p ALWAYS SPLITS (not ramifies, not inert)
\\ in K = Q(sqrt(sf)) when p is an odd prime with p ∤ a₂.
\\
\\ The proof is algebraic:
\\   - Inert excluded by Theorem (A)-(E): the only ideal of norm p² over an inert
\\     prime is (p), but (β)≠(p) rules that out.
\\   - Ramified excluded for odd p: p ramifies in Q(sqrt(sf)) iff p | disc(K),
\\     and for odd p, p | disc(K) iff p | sf. Since sf = sqfree(D) and
\\     D ≡ a₂² (mod p), we have p | D iff p | a₂ — excluded by hypothesis.
\\
\\ This script tests 10 NEW (p, a₂) pairs not in Thread 16's test set,
\\ explicitly checking the splitting type (SPLIT / INERT / RAMIFIED) and
\\ confirming that [P]² = 1 in each case.
\\
\\ Run: gp -q thread17_splits_verify.gp

default(parisize, 64000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Returns: 1=SPLIT, 0=INERT, -1=RAMIFIED
split_type(K, p) = {
  my(Pp = idealprimedec(K, p));
  if(#Pp == 2, 1,
     Pp[1][3] == 2, -1,
     0)
};

verify_splits(p, a2, label) = {
  my(D, sf, m2, m, K, stype, Pp, P, P2, P2_prin, order2);

  if(a2 % p == 0 || a2 == 0,
    printf("  SKIP %s: degenerate\n", label); return(-1));

  D  = a2^2 - 4*p^2;
  sf = sf_part(D);
  if(sf == 1,
    printf("  SKIP %s: D perfect square\n", label); return(-1));

  m2 = D / sf;
  m  = sqrtint(m2);
  if(m*m != m2, printf("  ERROR %s: D/sf not square\n", label); return(0));

  K  = bnfinit(x^2 - sf, 1);
  stype = split_type(K, p);

  Pp    = idealprimedec(K, p);
  P     = Pp[1];
  P2    = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);
  order2  = (P2_prin[1] == 0*P2_prin[1]);

  printf("  %-32s p=%-6d a2=%-8d sf=%-8d  split=%-8s  [P]²=1:%s%s\n",
    label, p, a2, sf,
    if(stype==1,"SPLIT",stype==-1,"RAMIFIED","INERT"),
    if(order2,"YES","FAIL(!)"),
    if(stype!=1," <-- PARADOX!",""));

  if(stype==1 && order2, 1, 0)
};

print("Thread 17: 10 new (p,a2) pairs — splitting type verification");
print("Expected: all SPLIT, all [P]²=1");
print();

{
  my(total=0, ok=0, r);
  my(cases = [
    \\ Chosen to be far from the Thread-16 test set; mix of D<0 and D>0
    [131,  43,    "new p=131 a2=43 (D<0)"],
    [137,  55,    "new p=137 a2=55 (D<0)"],
    [149,  61,    "new p=149 a2=61 (D<0)"],
    [151, -63,    "new p=151 a2=-63 (D<0)"],
    [157,  65,    "new p=157 a2=65 (D<0)"],
    [163,  67,    "new p=163 a2=67 (D<0)"],
    \\ D>0: |a2| > 2p
    [53,  107,    "new p=53 a2=107 (D>0)"],
    [61,  123,    "new p=61 a2=123 (D>0)"],
    [71,  143,    "new p=71 a2=143 (D>0)"],
    [79,  159,    "new p=79 a2=159 (D>0)"]
  ]);

  for(i=1, #cases,
    my(p=cases[i][1], a2=cases[i][2], lbl=cases[i][3]);
    r = verify_splits(p, a2, lbl);
    if(r >= 0, total++; ok += r));

  printf("\nResult: %d/%d passed (all SPLIT + [P]²=1)\n", ok, total);
  if(ok == total,
    print("CONFIRMED: p always splits in Q(sqrt(sf)) when p odd and p∤a2."),
    print("UNEXPECTED FAILURE — check output above."));
}
