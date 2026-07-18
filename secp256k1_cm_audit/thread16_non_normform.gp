\\ thread16_non_normform.gp
\\
\\ Thread 16: Verify order-2 Frobenius Theorem (Thread 15) for non-norm-form primes.
\\
\\ THEOREM (general): For any prime p and integer a2 with
\\   D = a2^2 - 4p^2 = sf*m^2  (sf squarefree negative, m>0, p does not divide a2)
\\ the prime P above p in K = Q(sqrt(sf)) satisfies [P]^2 = 1 in Cl(K).
\\
\\ Why p always SPLITS in K: D ≡ a2^2 (mod p) is a nonzero square; since
\\ D = sf*m^2 and p∤m (else p|D=a2^2, so p|a2, contradicting p∤a2), we get
\\ sf ≡ (a2*m^{-1})^2 (mod p) — a nonzero square — so (sf/p) = 1. □
\\
\\ Run: gp -q thread16_non_normform.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_normform(p) =
{
  my(v, w, s);
  v = 4*p - 73;
  if(v <= 0 || v % 3 != 0, return(0));
  w = v \ 3;
  s = sqrtint(w);
  (s*s == w)
}

\\ verify_one: check [P]^2=1 for one (p, a2) pair.
\\ Returns 0=skip, 1=trivial(h=1 pass), 2=non-trivial(h>=2 pass), -1=violation.
verify_one(p, a2, verbose) =
{
  my(D, sf, m2, m, K, h, Pp, P, P2, pr2, lbl, p2ok);

  D = a2^2 - 4*p^2;
  if(D >= 0, return(0));
  if(a2 % p == 0, return(0));

  sf = sf_part(D);
  if(sf >= 0, return(0));
  m2 = D / sf;
  if(denominator(m2) != 1 || m2 <= 0, return(0));
  m = sqrtint(m2);
  if(m*m != m2, return(0));

  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);

  if(#Pp != 2,
    if(verbose,
      printf("    p=%d a2=%d: p does not split (#Pp=%d) UNEXPECTED\n", p, a2, #Pp));
    return(0));

  P   = Pp[1];
  P2  = idealpow(K, P, 2);
  pr2 = bnfisprincipal(K, P2, 1);
  p2ok = (pr2[1] == 0*pr2[1]);

  if(verbose,
    lbl = if(h >= 2, Str(" [NON-TRIVIAL h=", h, "]"), "");
    printf("    p=%-8d a2=%-6d sf=%-9d m=%-5d h=%-4d  [P]^2=1:%s%s\n",
      p, a2, sf, m, h, if(p2ok,"YES","NO(!)"), lbl));

  if(!p2ok, return(-1));
  if(h >= 2, return(2));
  return(1)
}

\\ ------------------------------------------------------------------
\\ Main: 15 non-norm-form primes, 5 a2 values each
\\ ------------------------------------------------------------------

print("Thread 16: Universal order-2 Frobenius — generalisation to non-norm-form primes");
print("==================================================================================");
print();

{
  my(nnf, p, a2_vals, a2, r, cnt, ii, total, trivial, nontrivial, err);
  nnf = 0; p = 2; total = 0; trivial = 0; nontrivial = 0; err = 0;
  a2_vals = [1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];

  while(nnf < 15,
    p = nextprime(p + 1);
    if(is_normform(p), next);
    nnf++;
    printf("Non-norm-form prime #%d: p = %d\n", nnf, p);
    cnt = 0; ii = 1;
    while(cnt < 5 && ii <= #a2_vals,
      a2 = a2_vals[ii]; ii++;
      if(a2 % p == 0, next);
      r = verify_one(p, a2, 1);
      if(r == 0, next);
      total++;
      if(r ==  1, trivial++);
      if(r ==  2, nontrivial++);
      if(r == -1, err++);
      cnt++
    );
    print()
  );

  printf("=== SUMMARY ===\n");
  printf("Cases tested         : %d\n", total);
  printf("Trivial h=1 passes   : %d\n", trivial);
  printf("Non-trivial h>=2 pass: %d\n", nontrivial);
  printf("Violations           : %d\n", err);
  printf("Theorem holds        : %s\n", if(err == 0, "YES", "NO(!)"))
}

\\ ------------------------------------------------------------------
\\ Bonus: search for h>=4 cases (strong non-trivial verification)
\\ ------------------------------------------------------------------

print();
print("--- Bonus: first 8 non-norm-form cases with h(K) >= 4 ---");
print();

{
  my(found, p, a2, r, D, sf, m2, m, K, h);
  found = 0; p = 2;
  while(found < 8 && p < 50000,
    p = nextprime(p + 1);
    if(is_normform(p), next);
    forstep(a2 = 1, 300, 2,
      if(found >= 8, break);
      if(a2 % p == 0, next);
      D = a2^2 - 4*p^2;
      if(D >= 0, next);
      sf = sf_part(D);
      if(sf >= 0, next);
      m2 = D / sf;
      if(denominator(m2) != 1 || m2 <= 0, next);
      m  = sqrtint(m2);
      if(m*m != m2, next);
      K  = bnfinit(x^2 - sf, 1);
      h  = K.clgp.no;
      if(h < 4, next);
      r  = verify_one(p, a2, 0);
      if(r == 0, next);
      printf("  p=%-8d a2=%-5d sf=%-10d m=%-5d h=%-5d [P]^2=1:%s\n",
        p, a2, sf, m, h, if(r > 0, "YES", if(r == -1, "NO(!)", "SKIP")));
      found++
    )
  )
}

print();
print("DONE.");
