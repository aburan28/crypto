\\ thread17_splitting_verification.gp
\\
\\ Thread 17: Verify the Splitting Corollary from Threads 15-16.
\\
\\ THEOREM: Let p odd prime, a2 in Z, p∤a2, a2≠0.
\\   D=a2^2-4p^2, sf=squarefree(D). Then:
\\   (i)  p SPLITS in K=Q(sqrt(sf))  [neither inert nor ramified]
\\   (ii) [P]^2 = 0 in Cl(K)        [P any prime above p]
\\
\\ Algebraic proof:
\\   (i)  Inert impossible: N(beta)=p^2 and only ideal of norm p^2 when
\\        inert is (p), but (beta)=(p) => p|a2. Contradiction.
\\   (ii) Ramified impossible: p|disc(K)|4*sf, and sf|(a2^2-4p^2),
\\        so p|sf => p|a2^2 => p|a2. Contradiction.
\\   (iii) p splits. Ideals of norm p^2: P^2, (p)=PP̄, P̄^2.
\\         (beta)!=( p) => (beta)=P^2 or P̄^2. Both principal => [P]^2=0.
\\
\\ Run: gp --stacksize 64000000 -q thread17_splitting_verification.gp

default(parisize, 64000000);
default(timer, 0);

sf_part(n) = {
  if(n == 0, return(0));
  sign(n) * core(abs(n))
};

\\ Returns 1=PASS (splits and [P]^2=0), 0=FAIL, -1=SKIP
check_pair(pp, aa, lab) = {
  my(D, sf, bnf, dec, P2, chi, ok);
  D = aa^2 - 4*pp^2;
  sf = sf_part(D);
  if(sf == 0 || sf == 1, return(-1));
  bnf = bnfinit(x^2 - sf, 1);
  dec = idealprimedec(bnf, pp);
  if(#dec != 2,
    printf("  NOT-SPLIT %-32s sf=%d #dec=%d FAIL\n", lab, sf, #dec);
    return(0));
  P2  = idealmul(bnf, dec[1], dec[1]);
  chi = bnfisprincipal(bnf, P2);
  ok  = (norml2(chi[1]) == 0);
  printf("  SPLIT %-32s sf=%6d  h=%d  [P]^2=0: %s\n",
         lab, sf, bnf.clgp.no, if(ok, "YES PASS", "NO FAIL"));
  if(ok, 1, 0)
};

run_part_a() = {
  my(pa=0, fa=0, r);
  print("=== Part A: 10 targeted cases ==="); print("");
  r = check_pair(11,   4,   "p=11 a2=4");    if(r==1,pa++, if(r==0,fa++));
  r = check_pair(23,   8,   "p=23 a2=8");    if(r==1,pa++, if(r==0,fa++));
  r = check_pair(47,   20,  "p=47 a2=20");   if(r==1,pa++, if(r==0,fa++));
  r = check_pair(53,   10,  "p=53 a2=10 h=12"); if(r==1,pa++, if(r==0,fa++));
  r = check_pair(101,  40,  "p=101 a2=40");  if(r==1,pa++, if(r==0,fa++));
  r = check_pair(151,  60,  "p=151 a2=60");  if(r==1,pa++, if(r==0,fa++));
  r = check_pair(199,  78,  "p=199 a2=78");  if(r==1,pa++, if(r==0,fa++));
  r = check_pair(251,  100, "p=251 a2=100"); if(r==1,pa++, if(r==0,fa++));
  r = check_pair(503,  200, "p=503 a2=200"); if(r==1,pa++, if(r==0,fa++));
  r = check_pair(1009, 400, "p=1009 a2=400"); if(r==1,pa++, if(r==0,fa++));
  printf("\nPart A: %d/10 PASS  %d FAIL\n\n", pa, fa);
  [pa, fa]
};

chk_exc(pp, aa, lab) = {
  my(D, sf, bnf, dec);
  D  = aa^2 - 4*pp^2; sf = sf_part(D);
  if(sf == 0, printf("  DEGEN  %-30s\n", lab); return);
  bnf = bnfinit(x^2 - sf, 1);
  dec = idealprimedec(bnf, pp);
  if(#dec == 2,
    printf("  SPLIT    %-30s sf=%d h=%d (splits despite p|a2)\n", lab, sf, bnf.clgp.no),
  if(#dec == 1 && dec[1][3] == 2,
    printf("  RAMIFIED %-30s sf=%d\n", lab, sf),
    printf("  INERT    %-30s sf=%d\n", lab, sf)));
};

run_part_b() = {
  print("=== Part B: Violation cases where p|a2 ==="); print("");
  chk_exc(11, 11, "p=11 a2=11 p|a2");
  chk_exc(13, 13, "p=13 a2=13 p|a2");
  chk_exc(17, 34, "p=17 a2=34 p|a2");
  chk_exc(19, 19, "p=19 a2=19 p|a2");
  print("");
};

run_part_c() = {
  my(pc=0, fc=0, sc=0, D, sf, bnf, dec, P2, chi, ok);
  my(plist=[3,5,7,11,13,17,19,23,29,31,37,41,43,47,53]);
  my(a2list=[2,4,6,8,10,14,18,22]);
  print("=== Part C: Systematic sweep (15 primes x 8 a2) ==="); print("");
  for(ii = 1, #plist,
    for(jj = 1, #a2list,
      my(pp=plist[ii], aa=a2list[jj]);
      if(aa % pp == 0, sc++; next);
      D  = aa^2 - 4*pp^2; sf = sf_part(D);
      if(sf == 0 || sf == 1, sc++; next);
      bnf = bnfinit(x^2 - sf, 1);
      dec = idealprimedec(bnf, pp);
      if(#dec != 2,
        printf("  NOT-SPLIT p=%d a2=%d sf=%d\n", pp, aa, sf); fc++; next);
      P2  = idealmul(bnf, dec[1], dec[1]);
      chi = bnfisprincipal(bnf, P2);
      ok  = (norml2(chi[1]) == 0);
      if(ok, pc++,
        printf("  ORDER>2 p=%d a2=%d sf=%d h=%d\n", pp, aa, sf, bnf.clgp.no); fc++)));
  printf("Part C: %d PASS  %d FAIL  %d SKIP\n\n", pc, fc, sc);
  [pc, fc]
};

\\ === Main ===
ra = run_part_a();
run_part_b();
rc = run_part_c();

print_summary(ra, rc) = {
  print("=== SUMMARY ===");
  printf("Part A: %d/10 PASS\n", ra[1]);
  printf("Part C: %d PASS / %d FAIL\n", rc[1], rc[2]);
  if(ra[2] + rc[2] == 0,
    print("ALL PASS -- Splitting Corollary verified numerically."),
    printf("SOME FAILURES: A=%d C=%d\n", ra[2], rc[2]));
  print("Corollary: when p odd and p ndvd a2, p SPLITS in Q(sqrt(sf)) and [P]^2=0.");
};
print_summary(ra, rc);
