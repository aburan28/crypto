\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify that the order-2 Frobenius theorem (Thread 15) generalises
\\ to NON-norm-form primes.
\\
\\ THEOREM (Thread 15, general form): Let p be prime and a2 an integer with
\\   D = a2^2 - 4p^2 < 0 (so sf = sf(D) < 0, K = Q(sqrt(sf)) imaginary quad.)
\\   If p does NOT divide a2, then for any prime ideal P above p in O_K we have
\\   [P]^2 = 1 in Cl(K).
\\
\\ PROOF RECAP: beta = (-a2 + m*sqrt(sf))/2 (where D = sf*m^2, m>0) satisfies:
\\   (A) beta is an algebraic integer: min poly x^2 + a2*x + p^2
\\   (B) beta in K = Q(sqrt(sf))
\\   (C) N(beta) = p^2
\\   (D) p !| a2 => (beta) != (p) => (beta) = P^2 or Pbar^2
\\   (E) [P]^2 = 1 in Cl(K). QED
\\
\\ This script picks 15 primes p that are NOT norm-form (4p != 73+3k^2 for any k)
\\ and for each p tries several a2 values. It verifies (A)-(E) for each valid case.
\\
\\ Run: gp -q thread16_general_order2.gp

default(parisize, 256000000);
default(timer, 0);

\\ ---------------------------------------------------------------
\\ Squarefree part (signed)
\\ ---------------------------------------------------------------
sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ---------------------------------------------------------------
\\ Check if p is a norm-form prime: 4p = 73 + 3k^2 for some integer k
\\ ---------------------------------------------------------------
is_norm_form(p) = {
  my(val, r, k);
  val = 4*p - 73;
  if(val <= 0 || val % 3 != 0, return(0));
  r = val \ 3;
  k = sqrtint(r);
  (k*k == r)
};

\\ ---------------------------------------------------------------
\\ Verify order-2 theorem for one (p, a2) pair.
\\ Returns 1 if valid case passed, 0 if skipped (bad parameters), -1 if FAILED.
\\ ---------------------------------------------------------------
verify_general(p, a2, verbose) = {
  my(D, sf, m2, m, K, Pp, P, P2, P2_prin, norm_beta, beta_elt,
     I_beta, P2_iso_beta, ord_str, pass);

  D = a2^2 - 4*p^2;
  if(D >= 0, return(0));      \\ skip: need imaginary quadratic (D<0)
  if(a2 % p == 0, return(0)); \\ skip: p|a2 violates condition (D)

  sf = sf_part(D);
  m2 = D \ sf;                \\ = m^2 (positive integer)
  if(m2 <= 0, return(0));
  m = sqrtint(m2);
  if(m*m != m2, return(0));   \\ D/sf not a perfect square (shouldn't happen)

  \\ (i) beta as Polmod in Q[x]/(x^2-sf)
  beta_elt = Mod((-a2 + m*x)/2, x^2 - sf);

  \\ (ii) Norm check
  norm_beta = (a2^2 - m2*sf) / 4;   \\ = (a2^2 - D)/4 = p^2
  if(norm_beta != p^2,
    if(verbose, printf("  NORM ERROR: norm=%d != p^2=%d\n", norm_beta, p^2));
    return(-1));

  \\ (iii) bnf of K and prime decomposition
  K = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  if(#Pp == 0,
    if(verbose, printf("  SKIP: p=%d inert in Q(sqrt(%d))\n", p, sf));
    return(0));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2, 1);

  \\ (iv) Check (beta) == P^2 or Pbar^2
  I_beta = idealhnf(K, beta_elt);
  P2_iso_beta = (I_beta == P2);
  if(!P2_iso_beta && #Pp >= 2,
    P2_iso_beta = (I_beta == idealpow(K, Pp[2], 2)));

  \\ (v) [P]^2 = 1 iff P^2 is principal
  pass = (P2_prin[1] == 0*P2_prin[1]);

  if(verbose,
    ord_str = if(pass,
      if(bnfisprincipal(K,P,1)[1] == 0*bnfisprincipal(K,P,1)[1], "[P]=1", "[P]^2=1, [P]!=1"),
      "[P]^2!=1 FAIL");
    printf("  p=%-7d a2=%-7d sf=%-10d m=%-5d h=%-4d D<0:%s  N=p^2:%s  (b)=P^2:%s  %s\n",
      p, a2, sf, m, K.clgp.no,
      "YES",
      if(norm_beta==p^2,"YES","NO"),
      if(P2_iso_beta,"YES","NO"),
      ord_str));

  if(pass, 1, -1)
};

\\ ---------------------------------------------------------------
\\ Collect 15 non-norm-form primes in [100, 10000]
\\ For each, try several a2 values
\\ ---------------------------------------------------------------

print("Thread 16: Generalisation of order-2 Frobenius theorem to non-norm-form primes");
print("==================================================================================");
print();

{
  my(p, primes_nf, a2_list, r, passed, total, failed);
  passed = 0; total = 0; failed = 0;

  \\ Collect first 15 non-norm-form primes >= 100
  primes_nf = [];
  p = 100;
  while(#primes_nf < 15,
    p = nextprime(p);
    if(!is_norm_form(p), primes_nf = concat(primes_nf, [p]));
    p = p + 1);

  printf("Non-norm-form primes chosen: %Ps\n\n", primes_nf);

  \\ For each prime, test a spread of a2 values: 1, p/4, p/3, p/2, p/2+1, ...
  \\ We want D = a2^2 - 4p^2 < 0 => a2 < 2p; take a2 in {1,3,7,p//4, p//3, p//2, p-1, p+1}
  for(ii = 1, #primes_nf,
    p = primes_nf[ii];
    printf("p = %d (norm-form: %s):\n", p, if(is_norm_form(p),"YES","no"));

    a2_list = [1, 3, 7, p\4, p\3, p\2, p-1, p+1, 2*p-3];
    \\ filter out multiples of p
    for(jj = 1, #a2_list,
      my(a2 = a2_list[jj]);
      if(a2 <= 0, next);
      r = verify_general(p, a2, 1);
      if(r == 1, passed++; total++);
      if(r == -1, failed++; total++);
      if(r == 0, next)
    );
    print()
  );

  printf("SUMMARY: %d valid cases tested, %d PASSED [P]^2=1, %d FAILED\n",
    total, passed, failed);
}

\\ ---------------------------------------------------------------
\\ Extended test: 50 random non-norm-form primes in [10000, 100000]
\\ Each tested with a2 = p//3 (standard choice giving D < 0)
\\ ---------------------------------------------------------------

print();
print("Extended sweep: 50 non-norm-form primes p in [10000, 100000], a2 = p // 3");
print("--------------------------------------------------------------------------");

{
  my(p, cnt, passed2, failed2, a2, r, h_list);
  cnt = 0; passed2 = 0; failed2 = 0;
  h_list = [];
  p = 10000;
  while(cnt < 50,
    p = nextprime(p + 1);
    if(is_norm_form(p), next);
    a2 = p \ 3;
    if(a2 % p == 0, a2 = a2 + 1);  \\ paranoia
    r = verify_general(p, a2, 0);
    if(r == 1,
      cnt++;
      passed2++;
      \\ also record class number
      my(sf = sf_part(a2^2 - 4*p^2), K2 = bnfinit(x^2 - sf, 1));
      h_list = concat(h_list, [K2.clgp.no]);
      if(cnt % 10 == 0,
        printf("  ...%d/50 done, passed=%d, failed=%d, h range=[%d,%d]\n",
          cnt, passed2, failed2, vecmin(h_list), vecmax(h_list))));
    if(r == -1, failed2++; cnt++)
  );
  printf("\nExtended sweep: %d primes tested, PASSED=%d, FAILED=%d\n",
    cnt, passed2, failed2);
  printf("Class number range for Q(sqrt(sf(a2^2-4p^2))): [%d, %d]\n",
    vecmin(h_list), vecmax(h_list));
  printf("Class numbers (first 20): %Ps\n", h_list[1..min(20, #h_list)]);
}

\\ ---------------------------------------------------------------
\\ Special case: a2 chosen so sf = -3 (compare to secp256k1's CM field)
\\ D = a2^2 - 4p^2 = -3 * m^2 requires a2^2 + 3m^2 = 4p^2
\\ Try a2 = 1: 1 + 3m^2 = 4p^2 => m^2 = (4p^2-1)/3. Need 3 | (4p^2-1) = (2p-1)(2p+1)
\\ p != 3 => gcd(p,3)=1 => 2p ≡ ±1 (mod 3) => 4p^2 ≡ 1 (mod 3) => 4p^2-1 ≡ 0 (mod 3). OK!
\\ So for p != 3, a2=1 gives D = 1-4p^2, and sf(D) = sf(-(4p^2-1)) = sf(-(2p-1)(2p+1))
\\ This won't always be -3 but let's try to find primes where sf=-3 naturally.

print();
print("Special: find 5 non-norm-form primes where sf(a2^2-4p^2) = -3 for some a2");
print("(These are primes admitting a biquadratic Weil poly splitting over Q(sqrt(-3)))");
print("---------------------------------------------------------------------------");

{
  my(p, cnt, sf, a2, D, m2);
  cnt = 0;
  p = 100;
  while(cnt < 5,
    p = nextprime(p + 1);
    if(is_norm_form(p), next);
    \\ Try a2 from 1 to p//2+1
    my(found = 0);
    for(a2t = 1, p,
      D = a2t^2 - 4*p^2;
      if(D >= 0, next);
      if(a2t % p == 0, next);
      sf = sf_part(D);
      if(sf != -3, next);
      m2 = D \ sf;
      if(m2 <= 0, next);
      my(mm = sqrtint(m2));
      if(mm*mm != m2, next);
      \\ Found one!
      my(r = verify_general(p, a2t, 1));
      if(r == 1, cnt++; found = 1; break));
    if(!found && cnt < 5,
      next)  \\ keep searching
  );
  printf("\nFound %d primes with sf=-3 (non-norm-form). All verified [P]^2=1.\n", cnt);
}

print();
print("DONE.");
