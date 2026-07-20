\\ thread17_split_verification.gp
\\ Thread 17: Verify p always SPLITS (never inert, never ramified) in Q(sqrt(sf))
\\ when p is odd prime and p does not divide a2, D = a2^2 - 4p^2 = sf*m^2.
\\
\\ Corollary to Thread 15-16 Theorem:
\\   For odd prime p with p not dividing a2, a2 != 0, D != 0:
\\   kronecker(disc(Q(sqrt(sf))), p) = +1  [p splits].
\\
\\ Proof sketch of Corollary:
\\   INERT impossible:   p inert => only ideal of norm p^2 is (p);
\\                       but (beta) has norm p^2 and (beta) != (p) (since p !| a2). Contradiction.
\\   RAMIFIED impossible: p ramifies in Q(sqrt(sf)) iff p | disc.
\\                        For odd p, disc | 4*sf, so p | disc iff p | sf.
\\                        sf | D = a2^2 - 4p^2 => p | sf iff p | a2^2 iff p | a2. Excluded.
\\   SPLIT is the only remaining case. □
\\
\\ This script tests 30 (p, a2) pairs not covered in Thread 16.

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Test one pair: return 1=SPLITS, 0=error, -1=skip
test_split(p, a2) = {
  my(D, sf, disc_K, kron_sym);
  if(a2 % p == 0,
    printf("  SKIP p=%d a2=%d (p|a2)\n", p, a2);
    return(-1));
  D = a2^2 - 4*p^2;
  if(D == 0,
    printf("  SKIP p=%d a2=%d (D=0)\n", p, a2);
    return(-1));
  sf = sf_part(D);
  if(sf == 1 || sf == -1,
    printf("  SKIP p=%d a2=%d (D perfect square, K=Q)\n", p, a2);
    return(-1));
  disc_K = if(sf % 4 == 1, sf, 4*sf);
  kron_sym = kronecker(disc_K, p);
  printf("  p=%-7d a2=%-7d D=%-10d sf=%-7d disc=%-10d kron=%2d %s\n",
    p, a2, D, sf, disc_K, kron_sym,
    if(kron_sym ==  1, "[SPLITS OK]",
    if(kron_sym == -1, "[INERT -- THEOREM VIOLATION!]",
                       "[RAMIFIED -- THEOREM VIOLATION!]")));
  if(kron_sym == 1, 1, 0)
};

print("=== Thread 17: p-splits verification (30 new pairs) ===");
print();

{
  my(total=0, ok=0, r);
  my(cases = [
    \\ --- Batch 1: D < 0, imaginary quadratic, primes not in Thread 16 ---
    [101,  10],
    [103,  14],
    [107,  20],
    [151,  30],
    [157,  24],
    [163,  22],
    [197,  18],
    [199,  28],
    [211,  40],
    [223,  16],
    [227,  36],
    [229,  42],
    [233,  12],
    [239,  50],
    [241,  44],
    \\ --- Batch 2: D > 0, real quadratic, a2 > 2p ---
    [11,   24],
    [13,   28],
    [17,   36],
    [19,   40],
    [23,   48],
    [29,   60],
    [31,   64],
    [37,   76],
    [41,   84],
    [43,   88],
    \\ --- Batch 3: large primes ---
    [1009, 100],
    [1013, 120],
    [1019,  90],
    [2003, 200],
    [3001, 150]
  ]);

  for(i = 1, #cases,
    r = test_split(cases[i][1], cases[i][2]);
    if(r >= 0, total++; ok += r));

  print();
  printf("=== SUMMARY: %d/%d SPLITS ===\n", ok, total);
  if(ok == total,
    print("ALL PASS. Corollary confirmed: p always splits when p does not divide a2."),
    print("FAILURES DETECTED -- review above."));
}
