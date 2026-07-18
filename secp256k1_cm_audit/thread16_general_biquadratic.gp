\\ thread16_general_biquadratic.gp
\\
\\ Thread 16: Verify the general biquadratic Weil polynomial theorem
\\ for non-norm-form primes, including genuinely non-trivial h(K) > 1 cases.
\\
\\ Theorem (Thread 15, general form — no norm-form assumption):
\\   Let p prime, a2 integer with D = a2^2 - 4p^2 < 0 = sf*m^2 (sf squarefree, m>0).
\\   If p does not divide a2, then [P]^2 = 1 in Cl(K) for K = Q(sqrt(sf)).
\\
\\ Proof sketch (purely from the Weil polynomial; norm-form irrelevant):
\\   beta = (-a2 + m*sqrt(sf))/2  satisfies  x^2 + a2*x + p^2 = 0  (A,B)
\\   N_{K/Q}(beta) = p^2  =>  (beta) has ideal-norm p^2  (C)
\\   p !| a2  =>  beta/p not a unit  =>  (beta) != (p)  (D)
\\   Only options: (beta) = P^2 or Pbar^2  =>  [P]^2 = 1.  QED  (E)
\\
\\ Run: gp -q thread16_general_biquadratic.gp

default(parisize, 256000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_norm_form(p) = {
  my(tmp, kk);
  tmp = 4*p - 73;
  if(tmp < 0 || (tmp % 3) != 0, return(0));
  tmp = tmp \ 3;
  kk = sqrtint(tmp);
  return(kk*kk == tmp);
};

\\ verify_case(p, a2):
\\   Confirms [P]^2 = 1 in Cl(Q(sqrt(sf))) for given (p, a2).
\\   Returns [sfD, m, hK, is_confirmed] or a string on error/skip.
verify_case(p, a2) = {
  my(D, sfD, m2, m, K, hK, Plist, P, Psq, cl);
  D = a2^2 - 4*p^2;
  if(D >= 0, return("D>=0"));
  sfD = sf_part(D);
  if(sfD == 0, return("sfD=0"));
  m2 = D / sfD;
  if(m2 <= 0, return("m2<=0"));
  m = sqrtint(m2);
  if(m*m != m2, return("m2 not square"));
  if(a2 % p == 0, return("p|a2"));
  K = bnfinit(x^2 - sfD, 1);
  hK = K.clgp.no;
  if(hK == 1, return([sfD, m, 1, 1]));
  Plist = idealprimedec(K, p);
  if(#Plist == 0, return("no P above p"));
  P = Plist[1];
  Psq = idealpow(K, P, 2);
  cl = bnfisprincipal(K, Psq);
  \\ cl[1] = exponent vector in class group SNF; principal <=> zero vector
  return([sfD, m, hK, (cl[1] == 0*cl[1])]);
};

\\ ==============================================================
\\ Phase 1: 12 explicit non-norm-form primes (hand-chosen)
\\ ==============================================================
print("=== Thread 16: General biquadratic Weil polynomial theorem ===");
print("Checking [P]^2 = 1 in Cl(Q(sqrt(sf))) for non-norm-form primes.");
print("");
print("--- Phase 1: 12 explicit non-norm-form primes ---");
print(" p   |  t |  a2  |        D    |   sf   |  m  | h(K) | OK?");
print("-----|----|----- |-------------|--------|-----|------|----");
{
  my(cases, p, t, a2, res, all_pass);
  cases = [97,5; 101,7; 103,9; 107,3; 113,5; 127,13; 137,11; 139,7; 149,9; 151,5; 157,7; 167,13];
  all_pass = 1;
  for(ii = 1, matsize(cases)[1],
    p = cases[ii, 1];
    t = cases[ii, 2];
    a2 = 2*p - t^2;
    if(is_norm_form(p), printf(" %3d (norm-form skip)\n", p); next);
    res = verify_case(p, a2);
    if(type(res) == "t_STR",
      printf(" %3d| %2d|%5d: SKIP(%s)\n", p, t, a2, res);
      next
    );
    if(!res[4], all_pass = 0);
    printf(" %3d| %2d| %4d| %11d| %7d| %3d| %4d | %s\n",
      p, t, a2, a2^2-4*p^2, res[1], res[2], res[3],
      if(res[4], "YES", "FAIL"));
  );
  printf("Phase 1 all_pass = %d\n", all_pass);
}

\\ ==============================================================
\\ Phase 2: Systematic search for h(K) > 1 cases (p in [50,500])
\\ ==============================================================
print("");
print("--- Phase 2: Search for h(K)>1 cases (p in [50,500]) ---");
{
  my(found, p, t, a2, D, sfD, m2, m, hK, res, all_pass);
  found = List();
  forprime(p = 50, 500,
    if(is_norm_form(p), next);
    for(t = 1, 2*sqrtint(p),
      a2 = 2*p - t^2;
      D = a2^2 - 4*p^2;
      if(D >= 0, next);
      sfD = sf_part(D);
      if(sfD == 0, next);
      m2 = D / sfD;
      if(m2 <= 0, next);
      m = sqrtint(m2);
      if(m*m != m2, next);
      if(a2 % p == 0, next);
      hK = bnfinit(x^2 - sfD, 1).clgp.no;
      if(hK > 1,
        listput(found, [p, t, a2, sfD, m, hK]);
        if(#found >= 15, break(2))
      )
    )
  );
  printf("Found %d cases with h(K) > 1.\n\n", #found);
  print(" p   |  t |  a2  |    sf    |  m  | h(K) | OK?");
  print("-----|----|----- |----------|-----|------|----");
  all_pass = 1;
  for(jj = 1, #found,
    p  = found[jj][1];
    t  = found[jj][2];
    a2 = found[jj][3];
    sfD = found[jj][4];
    m  = found[jj][5];
    hK = found[jj][6];
    res = verify_case(p, a2);
    if(type(res) == "t_STR",
      printf(" %3d| %2d|%5d: SKIP(%s)\n", p, t, a2, res); next);
    if(!res[4], all_pass = 0);
    printf(" %3d| %2d| %4d| %9d| %3d| %4d | %s\n",
      p, t, a2, sfD, m, hK, if(res[4], "YES", "FAIL"));
  );
  printf("Phase 2 all_pass = %d\n", all_pass);
}

\\ ==============================================================
\\ Phase 3: Extended search targeting h(K) >= 3
\\ ==============================================================
print("");
print("--- Phase 3: Targeting h(K)>=3 cases (p in [50,2000]) ---");
{
  my(found3, p, t, a2, D, sfD, m2, m, hK, res, all_pass);
  found3 = List();
  forprime(p = 50, 2000,
    if(is_norm_form(p), next);
    for(t = 1, 2*sqrtint(p),
      a2 = 2*p - t^2;
      D = a2^2 - 4*p^2;
      if(D >= 0, next);
      sfD = sf_part(D);
      if(sfD == 0, next);
      m2 = D / sfD;
      if(m2 <= 0, next);
      m = sqrtint(m2);
      if(m*m != m2, next);
      if(a2 % p == 0, next);
      hK = bnfinit(x^2 - sfD, 1).clgp.no;
      if(hK >= 3,
        listput(found3, [p, t, a2, sfD, m, hK]);
        if(#found3 >= 10, break(2))
      )
    )
  );
  printf("Found %d cases with h(K)>=3.\n\n", #found3);
  if(#found3 > 0,
    print(" p   |  t |  a2  |    sf    |  m  | h(K) | OK?");
    print("-----|----|----- |----------|-----|------|----");
    all_pass = 1;
    for(kk = 1, #found3,
      p   = found3[kk][1];
      t   = found3[kk][2];
      a2  = found3[kk][3];
      sfD = found3[kk][4];
      m   = found3[kk][5];
      hK  = found3[kk][6];
      res = verify_case(p, a2);
      if(type(res) == "t_STR",
        printf(" %3d| %2d|%5d: SKIP(%s)\n", p, t, a2, res); next);
      if(!res[4], all_pass = 0);
      printf(" %3d| %2d| %4d| %9d| %3d| %4d | %s\n",
        p, t, a2, sfD, m, hK, if(res[4], "YES", "FAIL"));
    );
    printf("Phase 3 all_pass = %d\n", all_pass)
  ,
    print("No h(K)>=3 cases found in range.")
  )
}

print("");
print("Done.");
