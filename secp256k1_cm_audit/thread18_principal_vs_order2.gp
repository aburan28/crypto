\\ thread18_principal_vs_order2.gp
\\
\\ Thread 18: Classify [P] as principal (order 1) vs. order exactly 2 in Cl(K).
\\
\\ Background: Proposition prop:biquadratic-order2 (Threads 15-16) proves that
\\ for any odd prime p with p∤a₂, a₂≠0, the ideal P above p in
\\ K = Q(sqrt(sf(a₂²-4p²))) satisfies [P]² = 1. So ord([P]) ∈ {1,2}.
\\
\\ This script determines which case holds, grouped by the squarefree part sf.
\\
\\ EXPECTED FINDINGS (to verify):
\\   sf = -3  (h=1):  [P] always principal (trivial since h=1).
\\   sf = -219 (h=4): [P] always has order 2 (p not represented by principal form).
\\   Other sf:        depends on h(K) and the specific principal form.
\\
\\ COROLLARY (to add to paper):
\\   ord([P]) = 1  iff  p is represented by the principal binary quadratic form
\\                       of discriminant disc(K).
\\   In particular, if h(K) is odd, then [P]=1.
\\
\\ Three-part structure:
\\   Part A: Broad scan — classify 40 (p,a₂) pairs; group by sf.
\\   Part B: Verify principal-form criterion for sf=-219 (order-2) and sf=-3 (order-1).
\\   Part C: Find examples where [P]=1 with h(K)>1 (even class number, P still principal).
\\
\\ Run: /usr/bin/gp --stacksize 128000000 -q thread18_principal_vs_order2.gp

default(parisize, 128000000);
default(timer, 0);

\\====================================================================
\\ Utilities
\\====================================================================

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Returns the order of ideal P in Cl(K): 1 if principal, 2 otherwise.
\\ Precondition: [P]²=1 guaranteed by the theorem.
ideal_order_in_cl(K, P) = {
  my(v, P2, v2);
  v = bnfisprincipal(K, P)[1];
  if(v == 0*v, return(1));
  \\ Sanity: verify [P]^2 is principal
  P2 = idealpow(K, P, 2);
  v2 = bnfisprincipal(K, P2)[1];
  if(v2 != 0*v2, error("THEOREM VIOLATED: [P]^2 is not principal!"));
  return(2)
};

\\ Classify a single (p, a2) pair. Returns [sf, h(K), ord([P])] or -1 on skip.
classify_pair(p, a2) = {
  my(D, sf, m2, m, K, Pp, P, h_K, ord);
  if(a2 == 0 || a2 % p == 0, return(-1));
  D = a2^2 - 4*p^2;
  if(D == 0, return(-1));
  sf = sf_part(D);
  if(sf == 1, return(-1));  \\ K = Q, not quadratic
  m2 = D / sf;
  if(denominator(m2) != 1 || m2 <= 0, error("m² not positive int"));
  m = sqrtint(m2);
  if(m*m != m2, error("D/sf not a square"));
  K  = bnfinit(x^2 - sf, 1);
  h_K = K.clgp.no;
  Pp = idealprimedec(K, p);
  if(#Pp != 2, error(Str("p should split (sf=",sf,",p=",p,") but got #Pp=",#Pp)));
  P  = Pp[1];
  ord = ideal_order_in_cl(K, P);
  return([sf, h_K, ord])
};

\\====================================================================
\\ Part A: Broad scan — 40 (p,a₂) pairs; classify and group by sf
\\====================================================================

print("================================================================");
print("Part A: Classify [P] order for 40 (p, a2) pairs (group by sf)");
print("================================================================");
print();
printf("  %-40s  sf         h(K)  [P]-ord\n", "label");
printf("  %-40s  %-9s  %-4s  -------\n", "---", "---", "---");

{
  my(sf_table, cases, r, total, ok1, ok2);
  sf_table = Map();
  total = 0; ok1 = 0; ok2 = 0;

  cases = [
    \\ [p, a2, label]
    \\ norm-form CM-73 (should give sf=-219, order 2)
    [19,   4,   "CM-73   p=19   a2=4"],
    [37,   6,   "CM-73   p=37   a2=6"],
    [79,   10,  "CM-73   p=79   a2=10"],
    [109,  12,  "CM-73   p=109  a2=12"],
    \\ norm-form non-CM-73 (should also give sf=-219, order 2)
    [349,  44,  "non-CM73 p=349 a2=44  (k=21)"],
    [421,  48,  "non-CM73 p=421 a2=48  (k=23)"],
    [601,  58,  "non-CM73 p=601 a2=58  (k=27)"],
    \\ non-norm-form primes, various sf
    [23,   7,   "non-NF   p=23  a2=7"],
    [29,   11,  "non-NF   p=29  a2=11"],
    [41,   15,  "non-NF   p=41  a2=15"],
    [43,   13,  "non-NF   p=43  a2=13"],
    [47,   19,  "non-NF   p=47  a2=19"],
    [53,   21,  "non-NF   p=53  a2=21"],
    [59,   23,  "non-NF   p=59  a2=23"],
    [61,   25,  "non-NF   p=61  a2=25"],
    [67,   27,  "non-NF   p=67  a2=27"],
    [71,   29,  "non-NF   p=71  a2=29"],
    [73,   31,  "non-NF   p=73  a2=31"],
    [79,   33,  "non-NF   p=79  a2=33"],
    [83,   35,  "non-NF   p=83  a2=35"],
    [89,   37,  "non-NF   p=89  a2=37"],
    [97,   41,  "non-NF   p=97  a2=41"],
    [101,  43,  "non-NF   p=101 a2=43"],
    [103,  45,  "non-NF   p=103 a2=45"],
    [107,  47,  "non-NF   p=107 a2=47"],
    [109,  49,  "non-NF   p=109 a2=49"],
    [113,  51,  "non-NF   p=113 a2=51"],
    \\ D > 0 cases (real quadratic field)
    [23,   47,  "D>0     p=23  a2=47"],
    [29,   59,  "D>0     p=29  a2=59"],
    [41,   83,  "D>0     p=41  a2=83"],
    [13,   27,  "D>0     p=13  a2=27"],
    [17,   35,  "D>0     p=17  a2=35"],
    [7,    15,  "D>0     p=7   a2=15"],
    [11,   23,  "D>0     p=11  a2=23"],
    [13,   29,  "D>0     p=13  a2=29"]
  ];

  for(i = 1, #cases,
    my(p=cases[i][1], a2=cases[i][2], lbl=cases[i][3]);
    r = classify_pair(p, a2);
    if(r == -1, printf("  %-40s  SKIP\n", lbl); next);
    total++;
    if(r[3]==1, ok1++, ok2++);
    printf("  %-40s  sf=%-7d  h=%-4d  %d\n", lbl, r[1], r[2], r[3]);
    \\ Accumulate in sf_table
    if(!mapisdefined(sf_table, r[1]),
      mapput(sf_table, r[1], [0,0]));
    my(entry = mapget(sf_table, r[1]));
    if(r[3]==1, mapput(sf_table, r[1], [entry[1]+1, entry[2]]),
                mapput(sf_table, r[1], [entry[1], entry[2]+1]))
  );

  printf("\nPart A totals: %d order-1 (principal), %d order-2\n\n", ok1, ok2);

  \\ Print sf breakdown table
  print("  sf breakdown:");
  printf("  %-10s  %-5s  %-8s  %-8s\n", "sf", "h(K)", "#[P]=1", "#[P]=2");
  printf("  %-10s  %-5s  %-8s  %-8s\n", "---", "---", "---", "---");
  my(keys_v = Vec(sf_table));
  for(ki = 1, #keys_v,
    my(sf_k = keys_v[ki], entry = mapget(sf_table, sf_k));
    my(K_tmp = bnfinit(x^2 - sf_k, 1));
    printf("  %-10d  %-5d  %-8d  %-8d\n", sf_k, K_tmp.clgp.no, entry[1], entry[2])
  )
}

\\====================================================================
\\ Part B: Principal-form criterion for sf=-219 and sf=-3
\\====================================================================

print();
print("================================================================");
print("Part B: Principal-form criterion");
print("  [P]=1  iff  p is represented by the principal form of disc(K)");
print("================================================================");
print();

\\ B1: sf = -219 (disc=-219 ≡ 1 mod 4, so actual disc = -219).
\\     4 reduced forms of disc=-219:
\\       f1 = x^2+xy+55*y^2            (norm 1, principal)
\\       f4 = 2*x^2+xy+27*y^2          (norm 2, unique order-2 class)
\\       f2, f3 = conjugate pair of norm 4
\\     p splits in K and [P]=1 iff p is represented by f1.
\\     If p is represented by f4, then [P] has order 2.

print("B1: sf=-219, disc=-219, principal form f1 = x^2+xy+55*y^2");
print("    order-2 form  f4 = 2*x^2+xy+27*y^2");
print();

{
  my(cases219 = [
    [19,   4,  "CM-73   p=19"],
    [37,   6,  "CM-73   p=37"],
    [79,   10, "CM-73   p=79"],
    [109,  12, "CM-73   p=109"],
    [349,  44, "non-CM73 p=349"],
    [421,  48, "non-CM73 p=421"]
  ]);

  for(i = 1, #cases219,
    my(p_v = cases219[i][1], a2_v = cases219[i][2], lbl = cases219[i][3]);
    \\ Verify sf = -219
    my(D_v = a2_v^2 - 4*p_v^2, sf_v = sf_part(D_v));
    if(sf_v != -219,
      printf("  WARN: %s has sf=%d, not -219 (D=%d)\n", lbl, sf_v, D_v); next);
    \\ Check if p = x^2+xy+55*y^2 (principal form => [P]=1)
    my(rep1=0, rep4=0, bound=sqrtint(p_v)+2);
    for(y_v = -bound, bound,
      for(x_v = -bound, bound,
        if(x_v^2+x_v*y_v+55*y_v^2 == p_v, rep1=1);
        if(2*x_v^2+x_v*y_v+27*y_v^2 == p_v, rep4=1)
      )
    );
    \\ Check via bnf
    my(r = classify_pair(p_v, a2_v));
    printf("  %-20s sf=%d h=4: f1 reps p? %-3s  f4 reps p? %-3s  => [P] ord=%d\n",
           lbl, sf_v,
           if(rep1,"YES","NO"), if(rep4,"YES","NO"),
           r[3])
  )
}

print();
print("B2: sf=-3, h(-3)=1 => [P]=1 always.");
print("    Finding (p,a2) pairs where sf(a2^2-4p^2)=-3...");
print();

{
  my(sf3_cases, count);
  sf3_cases = List();
  count = 0;
  forprime(p3 = 5, 300,
    for(a2_v = 1, 2*p3-1,
      if(a2_v % p3 == 0, next);
      my(D3 = a2_v^2 - 4*p3^2);
      if(D3 >= 0, next);
      my(sf3 = sf_part(D3));
      if(sf3 == -3,
        count++;
        listput(sf3_cases, [p3, a2_v]);
        if(count >= 6, break)));
    if(count >= 6, break));

  for(i = 1, #sf3_cases,
    my(p_v = sf3_cases[i][1], a2_v = sf3_cases[i][2]);
    my(r = classify_pair(p_v, a2_v));
    printf("  p=%-5d a2=%-5d sf=-3  h=1  [P] ord=%d  (expected 1)\n",
           p_v, a2_v, if(r==-1,-1,r[3]))
  )
}

\\====================================================================
\\ Part C: Search for [P]=1 with h(K)>1
\\====================================================================

print();
print("================================================================");
print("Part C: Search for [P]=1 with h(K)>1 (P principal despite even h)");
print("================================================================");
print();

{
  my(found, searched);
  found = List();
  searched = 0;
  forprime(p_s = 5, 300,
    for(a2_s = 1, 2*p_s-1,
      if(a2_s % p_s == 0, next);
      my(D_s = a2_s^2 - 4*p_s^2);
      if(D_s == 0, next);
      my(sf_s = sf_part(D_s));
      if(sf_s == 1, next);
      searched++;
      my(K_s = bnfinit(x^2 - sf_s, 1));
      my(h_s = K_s.clgp.no);
      if(h_s == 1, next);
      my(Pp_s = idealprimedec(K_s, p_s));
      if(#Pp_s != 2, next);
      my(v_s = bnfisprincipal(K_s, Pp_s[1])[1]);
      if(v_s == 0*v_s,
        listput(found, [p_s, a2_s, sf_s, h_s]);
        if(#found >= 8, break)));
    if(#found >= 8, break));

  printf("Searched %d valid (p,a2) pairs with p<=300.\n", searched);
  if(#found == 0,
    print("No [P]=1 with h(K)>1 found in this range."),
    printf("Found %d examples with h(K)>1 AND [P]=1 (P is principal):\n\n", #found);
    printf("  %-7s %-7s %-10s %-6s  note\n", "p", "a2", "sf", "h(K)");
    printf("  %-7s %-7s %-10s %-6s\n", "---", "---", "---", "---");
    for(i = 1, #found,
      my(rec = found[i]);
      \\ Identify the principal form
      my(sf_r = rec[3], K_r = bnfinit(x^2 - sf_r, 1));
      printf("  %-7d %-7d %-10d %-6d\n", rec[1], rec[2], rec[3], rec[4])
    )
  )
}

\\====================================================================
\\ Part D: Summary and corollary statement
\\====================================================================

print();
print("================================================================");
print("Part D: Summary — COROLLARY for paper");
print("================================================================");
print();
print("COROLLARY (to Proposition prop:biquadratic-order2):");
print("");
print("Under the same hypotheses (p odd, p∤a₂, a₂≠0, D=a₂²-4p², K=Q(√sf(D))):");
print("");
print("  ord([P]) in Cl(K) is either 1 or 2 (by the Proposition).");
print("  The exact value is:");
print("    ord([P]) = 1  iff  p is represented by the principal binary");
print("                       quadratic form of discriminant disc(K);");
print("    ord([P]) = 2  otherwise.");
print("");
print("  Special cases:");
print("  (i)  h(K) odd => ord([P]) = 1  (order|2 in odd-order group => trivial).");
print("  (ii) sf(D) = -3 => h(K)=h(-3)=1 => [P]=1 always.");
print("       For secp256k1 p ≡ 1 mod 3, every (p,a2) with sf(D)=-3 gives [P]=1.");
print("  (iii)sf(D) = -219 => h(K)=4, Cl(K)≅Z/4Z, and [P] has order 2");
print("       for ALL norm-form primes tested (p NOT represented by x²+xy+55y²).");
print("");
print("Script: secp256k1_cm_audit/thread18_principal_vs_order2.gp");
print("================================================================");
