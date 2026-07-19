\\ thread17_p_splits.gp
\\
\\ Thread 17: Verify that p always SPLITS in Q(sqrt(sf)) when p∤a2 (p odd prime).
\\
\\ Proof outline (algebraic):
\\   - p ramifies in Q(sqrt(sf)) iff p | disc(Q(sqrt(sf))) iff p | sf
\\     (for odd prime p, since disc = sf or 4sf, and p odd => p|4sf iff p|sf).
\\     But sf = squarefree-part(a2^2 - 4p^2), and p | sf => p | a2^2 => p | a2.
\\     Excluded by hypothesis. So p cannot ramify.
\\   - p is inert in Q(sqrt(sf)) iff kronecker(sf, p) = -1.
\\     Proof by contradiction: if p inert, (p) is the unique prime of norm p^2,
\\     so (beta) = (p), giving beta = p*u (unit), Trace(beta) = -a2 = p*Trace(u),
\\     hence p | a2. Excluded. So p cannot be inert.
\\   => p SPLITS: kronecker(sf, p) = +1 (and p does not divide sf).
\\
\\ This script verifies the splitting claim numerically for 10 explicit (p, a2) pairs
\\ covering diverse cases: D<0, D>0, large sf, small sf, various class numbers.
\\
\\ Run: gp -q thread17_p_splits.gp

default(parisize, 64000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

check_split(p, a2) = {
  my(D, sf, kron_val, disc_K, ramified, inert, splits, ok);

  D = a2^2 - 4*p^2;
  if(D == 0, error("D=0 degenerate"));
  if(a2 % p == 0, error("p | a2 degenerate"));

  sf = sf_part(D);

  \\ Check p does not ramify: p should not divide sf
  ramified = (sf % p == 0);

  \\ Check splitting via Legendre/Kronecker symbol: kronecker(sf, p)
  \\ +1 => splits, -1 => inert, 0 => ramified (p|sf)
  kron_val = kronecker(sf, p);

  inert   = (kron_val == -1);
  splits  = (kron_val == 1);
  ok = !ramified && splits;

  printf("  p=%-7d a2=%-10d sf=%-8d kron(sf,p)=%-2d | %s\n",
    p, a2, sf, kron_val,
    if(ok, "SPLITS OK",
      if(ramified, "RAMIFIED (bug!)", "INERT (bug!)")));
  ok
};

print("Thread 17: Verify p splits in Q(sqrt(sf)) for 10 explicit (p,a2) pairs\n");
print("All pairs satisfy: p prime, a2 != 0, p ∤ a2\n");

results = List();

\\ --- D < 0 (imaginary quadratic) cases ---
print("--- D < 0 (imaginary quadratic K) ---");
listput(results, check_split(23, 4));     \\ D = 16-2116 < 0
listput(results, check_split(53, 10));    \\ D = 100-11236 < 0, h(K)=12 (large class group)
listput(results, check_split(101, 6));    \\ D = 36-40804 < 0
listput(results, check_split(257, 14));   \\ D = 196-264196 < 0
listput(results, check_split(1009, 44));  \\ D = 1936-4072324 < 0

\\ --- D > 0 (real quadratic) cases ---
print("--- D > 0 (real quadratic K) ---");
listput(results, check_split(7, 16));     \\ D = 256-196 = 60 > 0
listput(results, check_split(11, 24));    \\ D = 576-484 = 92 > 0
listput(results, check_split(13, 28));    \\ D = 784-676 = 108 > 0
listput(results, check_split(19, 40));    \\ D = 1600-1444 = 156 > 0
listput(results, check_split(47, 98));    \\ D = 9604-8836 = 768 > 0

\\ Summary
total = #results;
passed = sum(i=1, total, results[i]);
printf("\nResult: %d / %d pairs have p SPLITS in Q(sqrt(sf))\n", passed, total);
msg = if(passed == total, "ALL PASS", "FAIL");
printf("%s: %d/%d splits verified.\n", msg, passed, total);
