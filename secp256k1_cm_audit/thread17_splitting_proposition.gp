\\ thread17_splitting_proposition.gp
\\
\\ Thread 17: Verify "p always splits" corollary of Threads 15-16.
\\
\\ COROLLARY: Let p be an odd prime, a2 in Z with p!|a2, a2!=0.
\\   Set D = a2^2 - 4p^2, sf = squarefree part of D.
\\   Then kron(sf,p) = +1 (p splits in Q(sqrt(sf))).
\\   NOT inert (kron=-1) and NOT ramified (sf%p!=0).
\\
\\ Run: gp -q thread17_splitting_proposition.gp

default(parisize, 64000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Returns 1 if p splits, 0 if failure, -1 if skipped
check_splits(p, a2) = {
  my(D, sf, leg, ram, ok);
  D = a2^2 - 4*p^2;
  if(D == 0, return(-1));
  sf = sf_part(D);
  leg = kronecker(sf, p);
  ram = (sf % p == 0);
  ok = (leg == 1) && !ram;
  printf("  %s  p=%-7d a2=%-8d sf=%-12d  kron=%-2d  ram=%d\n",
    if(ok,"PASS","FAIL"), p, a2, sf, leg, ram);
  if(!ok, printf("  *** VIOLATION: p=%d a2=%d sf=%d leg=%d ram=%d\n", p, a2, sf, leg, ram));
  ok
};

\\=================================================================
\\ Part A: 15 explicit (p, a2) pairs from Thread 16 results
\\=================================================================

print("=== Part A: 15 pairs from Thread 16 ===");

{
  my(total=0, ok=0, r);
  my(cases = [
    [23,  10], [23,  -4],
    [29,  12], [29,  -6],
    [31,   8], [31, -14],
    [41,  16], [41, -20],
    [43,  18], [43, -22],
    [47,  22], [47,  -8],
    [53,  10], [53, -26],
    [59,  14]
  ]);
  for(i=1, #cases,
    my(p=cases[i][1], a2=cases[i][2]);
    r = check_splits(p, a2);
    if(r >= 0, total++; ok += r));
  printf("\nPart A: %d / %d passed\n\n", ok, total);
}

\\=================================================================
\\ Part B: norm-form primes (4p = 73 + 3k^2, k odd) with a2 = 2p-73
\\=================================================================

print("=== Part B: norm-form primes, a2 = 2p-73 ===");

{
  my(total=0, ok=0, r, p, a2);
  forstep(k = 1, 31, 2,
    my(v = 73 + 3*k^2);
    if(v % 4 != 0, next);
    p = v \ 4;
    if(!isprime(p), next);
    a2 = 2*p - 73;
    total++;
    ok += check_splits(p, a2));
  printf("\nPart B: %d / %d passed\n\n", ok, total);
}

\\=================================================================
\\ Part C: 10 larger primes p~10^6 with explicit a2
\\=================================================================

print("=== Part C: 10 larger primes ===");

{
  my(total=0, ok=0, r);
  my(cases = [
    [1000003,  55432],
    [1000033, -88760],
    [1000037, 123456],
    [1000079, -200002],
    [1000099, 444000],
    [1000117, -777777],
    [1000121, 314159],
    [1000133, -271828],
    [1000151, 161803],
    [1000159, -500000]
  ]);
  for(i=1, #cases,
    my(p=cases[i][1], a2=cases[i][2]);
    if(a2 % p == 0, next);
    total++;
    ok += check_splits(p, a2));
  printf("\nPart C: %d / %d passed\n\n", ok, total);
}

\\=================================================================
\\ Part D: Exhaustive sweep p in [3,97], a2 in [1,p-1]
\\   Confirm: kron(sf,p) != -1 (no inertness) and sf%p != 0 (no ramification)
\\=================================================================

print("=== Part D: Exhaustive sweep p in [3,97], a2 in [1,p-1] ===");

{
  my(primes_small, found_inert, found_ram, sweep_count, D, sf, leg);
  primes_small = primes(25);  \\ first 25 primes: 2,3,5,...,97
  found_inert = 0; found_ram = 0; sweep_count = 0;
  for(pi=2, #primes_small,  \\ skip p=2 (even prime, handled separately)
    my(p = primes_small[pi]);
    for(a2=1, p-1,
      D = a2^2 - 4*p^2;
      sf = sf_part(D);
      leg = kronecker(sf, p);
      sweep_count++;
      if(leg == -1,
        found_inert++;
        printf("  INERT: p=%d a2=%d sf=%d\n", p, a2, sf));
      if(sf % p == 0,
        found_ram++;
        printf("  RAMIF: p=%d a2=%d sf=%d\n", p, a2, sf))));
  printf("  Swept %d pairs (p odd prime <= 97, a2 in [1,p-1])\n", sweep_count);
  printf("  Inert  violations: %d  =>  %s\n", found_inert, if(found_inert==0,"NONE (CONFIRMED)","FAILURES"));
  printf("  Ramif  violations: %d  =>  %s\n", found_ram,   if(found_ram==0,  "NONE (CONFIRMED)","FAILURES"));
  printf("\nPart D: %s\n\n",
    if(found_inert==0 && found_ram==0,
       "PASS — p always splits (not inert, not ramified) for p<=97",
       "FAIL — check output above"));
}

\\=================================================================
\\ Summary
\\=================================================================

print("=== SUMMARY ===");
print("Corollary: for odd prime p with p!|a2, p SPLITS in Q(sqrt(sf(a2^2-4p^2))).");
print("Verified by: Parts A (15 pairs) + B (norm-form) + C (large) + D (exhaustive sweep).");

quit
