\\ thread16_general_theorem.gp
\\ Thread 16: Verify the order-2 theorem on non-norm-form primes.
\\
\\ THEOREM (general form of Thread 15):
\\   Let p be any prime, t an integer with 0 < |t| < 2*sqrt(p), p not dividing t.
\\   Set a2 = 2p - t^2, D = a2^2 - 4p^2 = t^2*(t^2-4p).
\\   Let sf = squarefree part of D (negative), m = sqrt(D/sf).
\\   Then the prime ideal P above p in O_K (K = Q(sqrt(sf))) satisfies [P]^2=1.
\\
\\ METHOD: For non-norm-form primes p (4p != 73+3k^2) and trace t in {2,3},
\\   use bnfinit + bnfisprincipal to verify [P]^2=1 directly.
\\
\\ Run: gp -q thread16_general_theorem.gp

default(parisize, 256000000);
default(timer, 0);

\\ Squarefree part preserving sign
sf_part(n) = {
  if(n == 0, return(0));
  sign(n) * core(abs(n))
};

\\ Is p a norm-form prime? (4p = 73 + 3k^2 for some odd k >= 1)
is_norm_form(p) = {
  my(v, k);
  v = 4*p - 73;
  if(v <= 0 || v % 3 != 0, return(0));
  v = v \ 3;
  if(!issquare(v), return(0));
  k = sqrtint(v);
  if(k % 2 == 0, return(0));
  1
};

\\ Verify [P]^2 = 1 for prime p, trace t.
\\ Returns: 1 (pass), 0 (fail), -1 (precond fails), -2 (D/sf not square).
check_order2(p, t) = {
  my(a2, D, sf, m2, K, fac, P2, res);
  if(t == 0 || t % p == 0, return(-1));
  if(t^2 >= 4*p, return(-1));
  a2 = 2*p - t^2;
  D  = a2^2 - 4*p^2;
  sf = sf_part(D);
  m2 = D / sf;
  if(!issquare(m2), return(-2));
  K  = bnfinit(x^2 - sf, 1);
  fac = idealprimedec(K, p);
  P2 = idealpow(K, fac[1], 2);
  res = bnfisprincipal(K, P2, 0);
  \\ res is the exponent vector in Cl(K)'s SNF basis.
  \\ Empty vector => trivial class group => P^2 is principal.
  \\ Non-empty => P^2 is principal iff all exponents are 0.
  if(#res == 0, return(1));
  if(sum(i=1, #res, abs(res[i])) == 0, 1, 0)
};

\\ Collect first 'count' primes > 100 that are NOT norm-form.
gen_non_norm(count) = {
  my(res, p, n);
  res = vector(count);
  p = 101; n = 0;
  while(n < count,
    if(isprime(p) && !is_norm_form(p),
      n++; res[n] = p);
    p = nextprime(p+1));
  res
};

\\ ========== Main ==========
print("Thread 16: order-2 theorem on NON-NORM-FORM primes");
print("====================================================");

{
  my(primes10, p, t, sf, r, pc, fc, sc);
  primes10 = gen_non_norm(10);
  pc = 0; fc = 0; sc = 0;
  print("p          t    sf           result");
  print("------------------------------------");
  for(ii = 1, 10,
    p = primes10[ii];
    for(jj = 1, 2,
      t = if(jj == 1, 2, 3);
      sf = sf_part((2*p - t^2)^2 - 4*p^2);
      r = check_order2(p, t);
      if(r ==  1, pc++; printf("%-10d %-4d %-12d YES\n", p, t, sf));
      if(r ==  0, fc++; printf("%-10d %-4d %-12d NO ***\n", p, t, sf));
      if(r  < 0, sc++; printf("%-10d %-4d (skip r=%d)\n", p, t, r))
    )
  );
  print("------------------------------------");
  printf("SUMMARY: %d pass, %d fail, %d skip\n", pc, fc, sc)
}

print("");
print("--- Large-prime spot checks t=2 ---");
{
  my(lp, p, sf, r);
  lp = [100003, 500009, 999983, 9999991];
  for(ii = 1, 4,
    p = lp[ii];
    if(!isprime(p), p = nextprime(p));
    if(is_norm_form(p), printf("p=%d is norm-form skip\n", p); next);
    sf = sf_part((2*p - 4)^2 - 4*p^2);
    r = check_order2(p, 2);
    if(r ==  1, printf("p=%-12d sf=%-12d [P]^2=1 YES\n", p, sf));
    if(r ==  0, printf("p=%-12d sf=%-12d [P]^2=1 NO ***\n", p, sf));
    if(r  < 0, printf("p=%-12d skip r=%d\n", p, r))
  )
}

print("");
print("--- Kronecker(sf,p)=1 check (p splits) ---");
{
  my(primes10, p, D, sf, leg);
  primes10 = gen_non_norm(10);
  for(ii = 1, 10,
    p = primes10[ii];
    D = (2*p - 4)^2 - 4*p^2;
    sf = sf_part(D);
    leg = kronecker(sf, p);
    if(leg == 1,
      printf("p=%-6d sf=%-8d Kron=%d (splits) OK\n", p, sf, leg),
      printf("p=%-6d sf=%-8d Kron=%d *** UNEXPECTED\n", p, sf, leg))
  )
}

\\ Algebraic identity spot check
{
  my(p0, t0, a20);
  p0 = 101; t0 = 5;
  a20 = 2*p0 - t0^2;
  printf("\nAlg-id check p=%d t=%d: (a2^2-4p^2)=%d == t^2*(t^2-4p)=%d ? %s\n",
    p0, t0,
    a20^2 - 4*p0^2,
    t0^2*(t0^2 - 4*p0),
    if(a20^2 - 4*p0^2 == t0^2*(t0^2 - 4*p0), "YES", "NO"))
}

print("");
print("Done.");
quit();
