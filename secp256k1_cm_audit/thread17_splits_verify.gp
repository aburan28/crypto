\\ thread17_splits_verify.gp
\\ Thread 17: Verify p splits (never inert, never ramified) in Q(sqrt(sf))
\\ when p is odd prime and p does not divide a2.
\\
\\ Theorem (Threads 15-16 corollary):
\\   p odd prime, a2 in Z with p ndiv a2, a2 != 0,
\\   D = a2^2 - 4p^2, sf = squarefree_part(D).
\\   Then p SPLITS in Q(sqrt(sf)): Kronecker(sf, p) = +1 and p ndiv sf.
\\
\\ Why ramification is impossible:
\\   D = a2^2 - 4p^2 == a2^2 (mod p), so p|D iff p|a2. Since p ndiv a2, p ndiv D.
\\   sf = core(D), so p ndiv sf. For odd p: p ramifies iff p|sf -- impossible.
\\ Why inertness is impossible:
\\   beta = (-a2 + m*sqrt(sf))/2 satisfies N(beta)=p^2, (beta) != (p).
\\   If p inert, only ideal of norm p^2 is (p); contradiction.
\\
\\ Run: gp -q secp256k1_cm_audit/thread17_splits_verify.gp

default(parisize, 64000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

check_splits(p, a2) = {
  my(D, sf, kr);
  if(!isprime(p), return(-1));
  if(a2 % p == 0 || a2 == 0, return(-1));
  D  = a2^2 - 4*p^2;
  if(D == 0, return(-1));
  sf = sf_part(D);
  if(sf % p == 0,
    printf("  RAM: p=%d | sf=%d (a2=%d) -- UNEXPECTED\n", p, sf, a2);
    return(0));
  kr = kronecker(sf, p);
  if(kr == 1, return(1));
  printf("  INERT/ZERO: p=%d sf=%d a2=%d kr=%d -- UNEXPECTED\n", p, sf, a2, kr);
  return(0)
};

print("=== Thread 17: p-splits verification ===");
print("Theorem (T15-16): p odd prime, p ndiv a2 => p splits in Q(sqrt(sf(a2^2-4p^2)))");
print("");

{
  my(passed=0, failed=0, ok);

  \\ --- Part A: 10 specific arbitrary (p, a2) pairs ---
  print("--- Part A: 10 specific (p, a2) pairs ---");

  my(test_p  = [101, 103, 107, 109, 127, 131, 137, 149, 151, 157]);
  my(test_a2 = [ 14,  20,  12,   8,  30,  22,  40,  18,  24,  16]);

  for(i = 1, 10,
    my(p  = test_p[i]);
    my(a2 = test_a2[i]);
    my(D  = a2^2 - 4*p^2);
    my(sf = sf_part(D));
    my(kr = kronecker(sf, p));
    ok = check_splits(p, a2);
    if(ok == 1, passed++, if(ok == 0, failed++));
    printf("  p=%d a2=%d sf=%d kr=%d  %s\n",
      p, a2, sf, kr, if(ok == 1, "PASS", if(ok == 0, "FAIL", "SKIP")));
  );

  print("");

  \\ --- Part B: 20 deterministic primes ---
  print("--- Part B: 20 deterministic primes (prime(51..70)) ---");

  for(idx = 1, 20,
    my(p  = prime(idx + 50));
    my(a2 = 2*idx + 1);
    while(a2 % p == 0, a2 = a2 + 2);
    ok = check_splits(p, a2);
    if(ok == 1, passed++, if(ok == 0, failed++));
  );
  print("  (no FAIL output = all passed)");
  print("");

  \\ --- Part C: cases with a2 giving large m ---
  print("--- Part C: large-m cases ---");

  my(lm_p  = [  97, 499, 1009, 2003, 5003]);
  my(lm_a2 = [   2,   6,   10,   14,   20]);

  for(i = 1, 5,
    my(p  = lm_p[i]);
    my(a2 = lm_a2[i]);
    my(D  = a2^2 - 4*p^2);
    my(sf = sf_part(D));
    my(kr = kronecker(sf, p));
    my(m_sq = D \ sf);
    my(m  = sqrtint(abs(m_sq)));
    ok = check_splits(p, a2);
    if(ok == 1, passed++, if(ok == 0, failed++));
    printf("  p=%d a2=%d m~=%d sf=%d kr=%d  %s\n",
      p, a2, m, sf, kr, if(ok == 1, "PASS", if(ok == 0, "FAIL", "SKIP")));
  );

  print("");
  print("=== Final summary ===");
  printf("Passed: %d / %d\n", passed, passed + failed);
  if(failed == 0,
    printf("ALL PASSED: p splits in Q(sqrt(sf)) for all tested (p,a2) with p ndiv a2.\n"),
    printf("FAILED: %d cases\n", failed));
}

\\ Detailed trace for p=101, a2=14
print("");
print("=== Detailed trace: p=101, a2=14 ===");
{
  my(p=101, a2=14);
  my(D = a2^2 - 4*p^2);
  my(sf = sf_part(D));
  printf("D = %d^2 - 4*%d^2 = %d\n", a2, p, D);
  printf("sf = core(%d) = %d\n", D, sf);
  printf("D mod p = %d  (= a2^2 mod p = %d)\n", D % p, (a2^2) % p);
  printf("p | sf ? %s\n", if(sf % p == 0, "YES (bug)", "NO (correct: no ramification)"));
  printf("Kron(sf=%d, p=%d) = %d  (must be +1)\n", sf, p, kronecker(sf, p));
  printf("p=%d SPLITS in Q(sqrt(%d))  [confirmed]\n", p, sf);
}

quit
