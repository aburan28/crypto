\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify the universal order-2 Frobenius theorem for
\\ NON-NORM-FORM (p, a2) pairs.
\\
\\ THEOREM (proved algebraically in Thread 15):
\\   Let p be prime, a2 integer with D = a2^2 - 4p^2 < 0 and p ∤ a2.
\\   Let sf = squarefree_part(D), m = sqrt(D/sf) > 0 (integer).
\\   Then beta = (-a2 + m*sqrt(sf))/2 in O_K has N(beta) = p^2.
\\   Since p ∤ a2: (beta) = P^2 for a prime P above p, so [P]^2 = 1 in Cl(K).
\\
\\ NEW (Thread 16): p ∤ a2 (with D < 0) IMPLIES p splits in K = Q(sqrt(sf)).
\\   Proof: ramified case: p|sf => p|D => p|a2^2 => p|a2 (contradicts p∤a2).
\\          inert case: (p) is the unique O_K-ideal of norm p^2; (beta)=(p) => p|a2.
\\   Therefore split hypothesis is AUTO-SATISFIED by p ∤ a2. □
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_order2.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_norm_form(p) = {
  my(r);
  r = 4*p - 73;
  if(r <= 0 || r % 3 != 0, return(0));
  r = r \ 3;
  sqrtint(r)^2 == r
};

split_type_str(K, p) = {
  my(facts, P1);
  facts = idealprimedec(K, p);
  if(#facts == 2, return("SPLIT"));
  P1 = facts[1];
  if(idealnorm(K, P1) == p, return("RAMIFIED"), return("INERT"))
};

\\ ---------------------------------------------------------------
\\ verify_general: full Thread-15-style check for arbitrary (p,a2).
\\ Returns 1 if [P]^2=1 (pass), 0 for skip/fail.
\\ ---------------------------------------------------------------
verify_general(p, a2, label) = {
  my(D, sf, m2, m, K, h, beta_elt, minpoly_beta, norm_beta,
     facts, nfacts, P, Pbar, P2, bp2, P2_ok, I_beta, P2_eq,
     stype, bP, P1_ok, ord_str);

  if(is_norm_form(p),
    printf("SKIP  %-18s p=%d is norm-form\n", label, p); return(0));
  if(!isprime(p),
    printf("SKIP  %-18s p=%d not prime\n", label, p); return(0));
  D = a2^2 - 4*p^2;
  if(D >= 0,
    printf("SKIP  %-18s D=%d>=0\n", label, D); return(0));
  if(a2 % p == 0,
    printf("SKIP  %-18s p|a2\n", label); return(0));

  sf = sf_part(D);
  m2 = D / sf;       \\ exact division; = m^2 (positive integer)
  m  = sqrtint(m2);
  if(m^2 != m2,
    printf("ERROR %-18s D/sf=%Ps not square\n", label, m2); return(0));

  beta_elt     = Mod((-a2 + m*x) / 2, x^2 - sf);
  minpoly_beta = minpoly(beta_elt);
  norm_beta    = (a2^2 - m2*sf) / 4;

  K     = bnfinit(x^2 - sf, 1);
  h     = K.clgp.no;
  stype = split_type_str(K, p);
  facts = idealprimedec(K, p);
  nfacts = #facts;
  P    = facts[1];
  Pbar = if(nfacts >= 2, facts[2], P);
  P2   = idealpow(K, P, 2);

  bp2   = bnfisprincipal(K, P2, 1);
  P2_ok = (bp2[1] == 0 * bp2[1]);

  I_beta = idealhnf(K, beta_elt);
  P2_eq  = (I_beta == P2 || I_beta == idealpow(K, Pbar, 2));

  bP    = bnfisprincipal(K, P, 1);
  P1_ok = (bP[1] == 0 * bP[1]);
  if(P2_ok,
    ord_str = if(P1_ok, "ord=1", "ord=2"),
    ord_str = "ord>2 FAIL");

  printf("%-18s p=%-5d a2=%-5d sf=%-7d m=%-5d h=%-3d %-9s Nβ=p²:%s (β)=P²:%s [P]²=1:%s [%s]\n",
    label, p, a2, sf, m, h, stype,
    if(norm_beta == p^2, "Y", "N"),
    if(P2_eq || stype != "SPLIT", "Y", "N"),
    if(P2_ok, "Y", "FAIL"),
    ord_str);

  P2_ok
};

\\ ---------------------------------------------------------------
\\ find_case: search for (p, a2) with sf_part(a2^2-4p^2) == sf_target
\\ p in [p_lo, p_max], a2 even. Returns [p, a2, m] or [0,0,0].
\\ ---------------------------------------------------------------
find_case(sf_target, p_lo, p_max) = {
  my(D, sfv, m2, mv, resa2, rp, rm);
  rp = 0; resa2 = 0; rm = 0;
  forprime(p = p_lo, p_max,
    if(rp > 0, break());
    if(is_norm_form(p), next());
    \\ quick pre-filter: p must split in Q(sqrt(sf_target))
    if(kronecker(sf_target, p) != 1, next());
    forstep(a2 = 2, 2*p - 2, 2,
      if(rp > 0, break());
      if(a2 % p == 0, next());
      D = a2^2 - 4*p^2;
      if(D >= 0, next());
      sfv = sf_part(D);
      if(sfv != sf_target, next());
      m2 = D / sf_target;    \\ = |D| / |sf_target| > 0, exact integer
      mv = sqrtint(m2);
      if(mv^2 == m2, rp = p; resa2 = a2; rm = mv)
    )
  );
  [rp, resa2, rm]
};

\\ ================================================================
print("Thread 16: Universal order-2 Frobenius — NON-NORM-FORM verification");
print("======================================================================");
print("Legend: Nβ=p²? (ii)  (β)=P²? (iv)  [P]²=1? (v) = main theorem check");
print();

\\ ----------------------------------------------------------------
\\ Part A: sf=-1, K=Q(i), h=1. Sanity check.
\\ Solve a2^2 + m^2 = 4p^2.
\\ p=5: (a2,m)=(6,8):   36+64=100=4*25  ✓
\\ p=13: (a2,m)=(10,24): 100+576=676=4*169 ✓
\\ p=29: (a2,m)=(42,40): 1764+1600=3364=4*841 ✓
\\ p=101: (a2,m)=(198,40): 39204+1600=40804=4*10201 ✓
print("--- Part A: sf=-1, Q(i), h=1 ---");
{
  my(pass=0, total=0, r);
  r=verify_general(5,   6,   "A1:p=5");    total++; if(r, pass++);
  r=verify_general(13,  10,  "A2:p=13");   total++; if(r, pass++);
  r=verify_general(29,  42,  "A3:p=29");   total++; if(r, pass++);
  r=verify_general(101, 198, "A4:p=101");  total++; if(r, pass++);
  printf("Part A: %d/%d\n\n", pass, total);
}

\\ ----------------------------------------------------------------
\\ Part B: sf=-5, K=Q(sqrt(-5)), h=2. Non-trivial [P]^2=1.
\\ Solve a2^2 + 5m^2 = 4p^2.
\\ p=41: (62,24): 3844+2880=6724=4*1681 ✓
\\ p=61: (58,48): 3364+11520=14884=4*3721 ✓
\\ p=269: (38,240): 1444+288000=289444=4*72361 ✓
print("--- Part B: sf=-5, Q(sqrt(-5)), h=2 ---");
{
  my(pass=0, total=0, r, cas);
  r=verify_general(41,  62,  "B1:p=41,sf=-5");   total++; if(r, pass++);
  r=verify_general(61,  58,  "B2:p=61,sf=-5");   total++; if(r, pass++);
  r=verify_general(269, 38,  "B3:p=269,sf=-5");  total++; if(r, pass++);
  \\ One more via search
  cas = find_case(-5, 300, 600);
  if(cas[1] > 0,
    r = verify_general(cas[1], cas[2], Str("B4:p=",cas[1],",sf=-5"));
    total++; if(r, pass++)
  );
  printf("Part B: %d/%d\n\n", pass, total);
}

\\ ----------------------------------------------------------------
\\ Part C: sf=-23, K=Q(sqrt(-23)), h=3.
\\ Solve a2^2 + 23m^2 = 4p^2.
\\ p=59: (26,24): 676+13248=13924=4*3481=4*59^2 ✓
\\ Search for more.
print("--- Part C: sf=-23, Q(sqrt(-23)), h=3 ---");
{
  my(pass=0, total=0, r, cas);
  r=verify_general(59, 26, "C1:p=59,sf=-23");  total++; if(r, pass++);
  cas = find_case(-23, 100, 800);
  if(cas[1] > 0,
    r = verify_general(cas[1], cas[2], Str("C2:p=",cas[1],",sf=-23"));
    total++; if(r, pass++)
  ,
    print("C2: no sf=-23 case found in [100,800]")
  );
  printf("Part C: %d/%d\n\n", pass, total);
}

\\ ----------------------------------------------------------------
\\ Part D: higher class numbers — search for sf=-14 (h=4) and sf=-47 (h=5).
print("--- Part D: higher h --- sf=-14 (h=4) and sf=-47 (h=5) ---");
{
  my(pass=0, total=0, r, cas);
  cas = find_case(-14, 50, 1000);
  if(cas[1] > 0,
    r = verify_general(cas[1], cas[2], Str("D1:sf=-14,p=",cas[1]));
    total++; if(r, pass++)
  ,
    print("D1: no sf=-14 case found")
  );
  cas = find_case(-47, 50, 1000);
  if(cas[1] > 0,
    r = verify_general(cas[1], cas[2], Str("D2:sf=-47,p=",cas[1]));
    total++; if(r, pass++)
  ,
    print("D2: no sf=-47 case found")
  );
  printf("Part D: %d/%d\n\n", pass, total);
}

\\ ----------------------------------------------------------------
\\ Part E: Structural equivalence (theoretical).
print("--- Part E: Structural equivalence (no computation) ---");
print("CLAIM: for D=a2^2-4p^2<0, the condition 'p∤a2' is EQUIVALENT to 'p splits in K=Q(sqrt(sf))'.");
print("  Ramified (p|sf): p|D=a2^2-4p^2 => p|a2. Contradicts p∤a2.");
print("  Inert: (p) is the unique O_K-prime of norm p^2; N(beta)=p^2 => (beta)=(p) => p|a2.");
print("  Split (Leg(sf,p)=1): No constraint forces p|a2; theorem gives [P]^2=1.");
print("COROLLARY: the split condition in the theorem is auto-satisfied by p∤a2.");
print();

\\ ----------------------------------------------------------------
\\ Part F: Systematic sweep, 20 non-norm-form primes, a2=2.
\\ For a2=2, p>2: D=-4(p^2-1)<0, p∤2; all conditions met.
\\ p^2-1 = (p-1)(p+1) always divisible by 8; sf=-squarefree_part(p^2-1), m=2*sqrt((p^2-1)/s).
print("--- Part F: 20 non-norm-form primes, a2=2 fixed ---");
{
  my(pass=0, fail=0, count=0, lbl, r);
  forprime(p=7, 10000,
    if(count >= 20, break());
    if(is_norm_form(p), next());
    lbl = Str("F", count+1, ":p=", p);
    r = verify_general(p, 2, lbl);
    if(r, pass++; count++, fail++; count++)
  );
  printf("\nPart F: %d/%d passed\n\n", pass, pass+fail);
}

print("=== DONE ===");
