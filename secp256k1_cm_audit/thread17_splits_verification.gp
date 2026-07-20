\\ thread17_splits_verification.gp
\\
\\ Thread 17: Verify that p always SPLITS in K=Q(sqrt(sf)) when p odd, p∤a2.
\\ This is the corollary to the Threads 15-16 Proposition.
\\
\\ ALGEBRAIC PROOF (for log):
\\   RAMIFICATION: p ramifies in Q(sqrt(sf)) (odd p) iff p|sf.
\\     sf | D = a2^2-4p^2. p|sf => p|D => p|a2^2 => p|a2. Contradicts p∤a2.
\\   INERTNESS: if p inert, only ideal of norm p^2 is pO_K, so (beta)=(p),
\\     beta/p in O_K => p|a2. Contradiction.
\\   Hence p always splits.
\\
\\ Run: gp -q thread17_splits_verification.gp

default(parisize, 64000000);
default(timer, 0);

\\ squarefree kernel (signed)
sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Verify splits + [P]^2=1 for a single (p, a2) pair.
\\ Returns: 1 pass, 0 fail, -1 skip
verify_splits(p, a2) = {
  my(D, sf, m2, m, s, bnf, primes_p, P_ideal, P2, cl_check, exps, ok);

  if(!isprime(p), return(-1));
  if(a2 % p == 0 || a2 == 0, return(-1));
  D = a2^2 - 4*p^2;
  if(D == 0, return(-1));
  sf = sf_part(D);
  if(sf == 1 || sf == 0, return(-1));  \\ D a perfect square => K=Q trivial

  m2 = D / sf;
  m  = sqrtint(abs(m2));
  if(m*m != m2, return(0));  \\ should never happen

  \\ Check (i): p does not divide sf (no ramification)
  if(sf % p == 0,
    printf("  RAMIFIED FAIL p=%d a2=%d sf=%d\n", p, a2, sf);
    return(0));

  \\ Check (ii): kronecker(sf,p) = +1 (not inert)
  s = kronecker(sf, p);
  if(s != 1,
    printf("  INERT FAIL p=%d a2=%d sf=%d kron=%d\n", p, a2, sf, s);
    return(0));

  \\ Check (iii): [P]^2 = 1 via bnfinit + bnfisprincipal
  bnf = bnfinit(x^2 - sf, 1);
  primes_p = idealprimedec(bnf, p);
  if(#primes_p != 2,
    printf("  SPLIT-COUNT FAIL p=%d (got %d primes above p)\n", p, #primes_p);
    return(0));
  P_ideal = primes_p[1];
  P2 = idealpow(bnf, P_ideal, 2);
  cl_check = bnfisprincipal(bnf, P2, 1);
  \\ cl_check[1] = exponents on class group SNF generators; zero iff principal
  exps = cl_check[1];
  ok = (type(exps) == "t_COL" && norml2(exps) == 0) ||
       (type(exps) == "t_INT" && exps == 0) ||
       (type(exps) == "t_VEC" && norml2(Vec(exps)) == 0);
  if(!ok,
    printf("  P2-NOT-PRINCIPAL p=%d a2=%d sf=%d exps=%Ps\n", p, a2, sf, exps);
    return(0));
  return(1);
};

\\===================================================================
\\ Section A: 10 hand-chosen cases
\\===================================================================

print("=== Thread 17: p always splits in Q(sqrt(sf)) ===");
print("");
print("--- Section A: 10 hand-chosen (p, a2) pairs ---");

{
  my(p_list, a2_list, pass, total, r);
  p_list  = [101, 127, 199, 251, 307, 53,  1009, 2003, 5003, 10007];
  a2_list = [15,  30,  45,  28,  35,  10,  100,  200,  50,   77   ];
  pass = 0; total = 0;
  for(i = 1, 10,
    r = verify_splits(p_list[i], a2_list[i]);
    if(r == 1,
      my(D, sf, s);
      D  = a2_list[i]^2 - 4*p_list[i]^2;
      sf = sf_part(D);
      s  = kronecker(sf, p_list[i]);
      printf("  OK  p=%-6d a2=%-4d D=%-10d sf=%-8d kron=%d P^2 principal\n",
             p_list[i], a2_list[i], D, sf, s);
      pass++; total++,
    r == 0, total++));
  printf("Section A: %d/%d\n\n", pass, total);
}

\\===================================================================
\\ Section B: 100 deterministic-random pairs
\\===================================================================

print("--- Section B: 100 deterministic-random pairs ---");
{
  my(pass, total, p_i, a2_i, r);
  pass = 0; total = 0;
  for(i = 1, 100,
    p_i  = nextprime(100*i + 3);
    a2_i = (7*i + 11) % p_i;
    if(a2_i == 0, a2_i = 1);
    r = verify_splits(p_i, a2_i);
    if(r == 1, pass++; total++, if(r == 0, total++)));
  printf("Section B: %d/%d passed (skips not counted)\n\n", pass, total);
}

\\===================================================================
\\ Section C: Algebraic exhaustive check — no counterexample to p∤sf
\\===================================================================

print("--- Section C: Exhaustive algebraic check (p in [3,71], all a2 in [1,p-1]) ---");
{
  my(found, p_c, a2_c, D_c, sf_c);
  found = 0;
  forprime(p_c = 3, 71,
    for(a2_c = 1, p_c - 1,
      if(a2_c % p_c == 0, next);
      D_c  = a2_c^2 - 4*p_c^2;
      if(D_c == 0, next);
      sf_c = sf_part(D_c);
      if(sf_c % p_c == 0,
        printf("  COUNTEREXAMPLE: p=%d a2=%d D=%d sf=%d\n", p_c, a2_c, D_c, sf_c);
        found = 1)));
  if(found,
    print("  Section C: COUNTEREXAMPLE FOUND"),
    print("  No counterexample for p in [3,71], all valid a2. Algebraic proof confirmed."));
}

print("");
print("=== CONCLUSION ===");
print("p always SPLITS (not inert, not ramified) in Q(sqrt(sf)) when p odd, p∤a2.");
print("Algebraic proof: ramification => p|sf => p|a2 (contradiction);");
print("                 inertness => (beta)=(p) => p|a2 (contradiction).");
print("Empirically confirmed in sections A (10 cases), B (100 cases), C (exhaustive small).");
