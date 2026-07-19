\\ thread16_general_weil.gp
\\
\\ Thread 16: Generalize Thread 15's order-2 theorem to non-norm-form primes.
\\
\\ THEOREM (proved in Thread 15):
\\   For any prime p, integer a2 with p not| a2, and D = a2^2 - 4p^2 < 0:
\\   set sf = sqfree_part(D), m = sqrt(D/sf), K = Q(sqrt(sf)).
\\   Let P be a prime ideal above p in O_K.
\\   Then [P]^2 = 1 in Cl(K).
\\
\\ The proof uses ONLY: beta = (-a2 + m*sqrt(sf))/2 satisfies x^2+a2*x+p^2=0,
\\ N(beta)=p^2, and p not| a2 => (beta) != (p) => (beta) = P^2 or Pbar^2.
\\ No dependence on the norm-form condition 4p=73+3k^2.
\\
\\ Run: gp -q thread16_general_weil.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Returns 1 if theorem holds, 0 if it fails, -1 if test skipped
verify_case(p, a2) = {
  my(D, sf, m2, m, K, h, fK, P, P2, beta_elt, minp, nbeta, pc);

  if(a2 == 0 || a2 % p == 0, return(-1));
  D = a2^2 - 4*p^2;
  if(D >= 0, return(-1));          \\ need imaginary quadratic

  sf = sf_part(D);
  m2 = D / sf;
  if(!issquare(m2, &m), return(-1));  \\ should always hold

  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  fK = idealprimedec(K, p);
  if(#fK == 0, return(-1));

  P  = fK[1];
  P2 = idealpow(K, P, 2);

  beta_elt = Mod((-a2 + m*x)/2, x^2 - sf);
  minp  = minpoly(beta_elt);
  nbeta = (a2^2 - m2*sf) / 4;       \\ should equal p^2

  pc = bnfisprincipal(K, P2);
  my(is_prin = (pc[1] == 0*pc[1]));

  printf("p=%-6d a2=%-5d sf=%-9d m=%-5d h=%-4d nP=%d minp=%Ps N(b)=p^2:%s [P]^2=1:%s\n",
    p, a2, sf, m, h, #fK, minp,
    if(nbeta == p^2, "YES", "NO!"),
    if(is_prin, "YES", "FAIL!!!"));

  return(if(is_prin, 1, 0));
};

\\ Norm-form primes from 4p=73+3k^2 (secp256k1 family) — exclude from tests
norm_form_ps = Set([19,37,79,109,349,487,739,937,1279,2287,3187,5437,6229,7369,7669,8287,8929,10639,12889,13687,14929,15787,18979,24049,28537]);

print("Thread 16: Order-2 Frobenius theorem -- non-norm-form verification");
print("====================================================================");
print();

\\ ---------------------------------------------------------------------------
\\ Phase 1: Non-trivial class group (h >= 2)
\\ ---------------------------------------------------------------------------
print("Phase 1: h >= 2 cases (non-trivial class group)");
print("------------------------------------------------");

{
  my(total=0, passed=0, cases=[]);

  forprime(p=2, 500,
    if(setsearch(norm_form_ps, p), next);
    forstep(a2=1, min(15, 2*p-2), 1,
      if(a2 % p == 0, next);
      my(D2 = a2^2 - 4*p^2);
      if(D2 >= 0, next);
      my(sf2 = sf_part(D2));
      my(m22 = D2 / sf2);
      if(!issquare(m22), next);
      my(K2 = bnfinit(x^2 - sf2, 1));
      if(K2.clgp.no < 2, next);
      cases = concat(cases, [[p, a2]]);
      if(#cases >= 25, break(2));
    );
  );

  printf("Found %d cases with h >= 2.\n", #cases);
  for(i=1, #cases,
    my(r = verify_case(cases[i][1], cases[i][2]));
    if(r >= 0, total++; if(r, passed++));
  );
  printf("\nPhase 1 summary: %d/%d passed\n\n", passed, total);
}

\\ ---------------------------------------------------------------------------
\\ Phase 2: First 15 non-norm-form primes, a2 in {1,2,3,4,5}
\\ ---------------------------------------------------------------------------
print("Phase 2: First 15 non-norm-form primes, a2 in {1,2,3,4,5}");
print("-----------------------------------------------------------");

{
  my(total=0, passed=0, cnt=0);
  forprime(p=2, 10000,
    if(setsearch(norm_form_ps, p), next);
    cnt++;
    if(cnt > 15, break);
    for(a2=1, 5,
      my(r = verify_case(p, a2));
      if(r >= 0, total++; if(r, passed++));
    );
  );
  printf("\nPhase 2 summary: %d/%d passed\n\n", passed, total);
}

\\ ---------------------------------------------------------------------------
\\ Phase 3: Larger primes near 1000, boundary a2 = p-1
\\ ---------------------------------------------------------------------------
print("Phase 3: Large primes {997,1009,1013,1019,1021}, a2 in {1,2,3,p-1}");
print("---------------------------------------------------------------------");

{
  my(total=0, passed=0);
  my(big_ps = [997, 1009, 1013, 1019, 1021]);
  for(i=1, #big_ps,
    my(p = big_ps[i]);
    if(setsearch(norm_form_ps, p) || !isprime(p), next);
    for(j=1, 4,
      my(a2 = if(j<=3, j, p-1));
      my(r = verify_case(p, a2));
      if(r >= 0, total++; if(r, passed++));
    );
  );
  printf("\nPhase 3 summary: %d/%d passed\n\n", passed, total);
}

\\ ---------------------------------------------------------------------------
\\ Phase 4: p=1999, a2 in 1..8 — check for large class numbers
\\ ---------------------------------------------------------------------------
print("Phase 4: p=1999 (non-norm-form), a2 in 1..8");
print("---------------------------------------------");

{
  my(total=0, passed=0);
  my(p = 1999);
  if(!isprime(p) || setsearch(norm_form_ps, p),
    print("skipped"),
    for(a2=1, 8,
      my(r = verify_case(p, a2));
      if(r >= 0, total++; if(r, passed++));
    );
    printf("\nPhase 4 summary: %d/%d passed\n", passed, total));
}

print();
print("DONE.");
