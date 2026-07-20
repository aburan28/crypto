\\ thread17_p_splits.gp
\\
\\ Thread 17: Verify that p SPLITS (Kronecker (sf/p) = +1, never inert = -1)
\\ in Q(sqrt(sf)) whenever p prime, a2 integer, p∤a2, a2≠0.
\\
\\ The algebraic proof:
\\   D = a2^2 - 4p^2, sf = squarefree part of D.
\\   (1) p ramifies in Q(sqrt(sf)) iff p|sf.
\\       p|sf => p|D = a2^2 - 4p^2 => p|a2^2 => p|a2.  Contradiction.
\\       So p does NOT ramify.
\\   (2) p inert iff Kronecker (sf/p) = -1.
\\       If p inert, N((p)) = p^2 and (p) is prime in O_K.
\\       But beta satisfies N(beta) = p^2 and (beta) ≠ (p) [since p∤a2].
\\       So (beta) and (p) are two distinct ideals of norm p^2 in O_K.
\\       If p is inert, p*O_K is the unique prime above p with norm p^2
\\       => every element of norm p^2 generates p*O_K => (beta) = (p). Contradiction.
\\       So p does NOT become inert.
\\   (3) p does not ramify (1) and does not become inert (2) => p SPLITS.
\\
\\ EMPIRICAL CHECK: 10 (p, a2) pairs, both norm-form and non-norm-form.
\\ Expected: Kronecker((sf/p)) = +1 in all 10 cases.
\\
\\ Run: gp -q thread17_p_splits.gp

default(parisize, 64000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

check_splits(p, a2, label) = {
  my(D, sf, kr, outcome);
  if(a2 % p == 0 || a2 == 0,
    printf("  SKIP  %s: p|a2\n", label); return);
  D  = a2^2 - 4*p^2;
  if(D == 0, printf("  SKIP  %s: D=0\n", label); return);
  sf = sf_part(D);
  \\ p|sf => p ramifies: already excluded by hypothesis
  \\ Kronecker symbol (sf/p): +1 = splits, -1 = inert, 0 = ramifies
  kr = kronecker(sf, p);
  if(kr == 1,  outcome = "SPLIT  ✓");
  if(kr == -1, outcome = "INERT  ✗");
  if(kr == 0,  outcome = "RAMIFY ✗");
  printf("  p=%-8d a2=%-8d sf=%-8d kr=%2d  %s  [%s]\n",
    p, a2, sf, kr, outcome, label);
  kr
};

print("=== Thread 17: p always splits in Q(sqrt(sf)) when p∤a2 ===");
print();

results = List();

\\ Test group A: norm-form primes (secp256k1 family, p ≡ 1 mod 6)
print("-- Group A: secp256k1 norm-form primes --");
listput(results, check_splits(7,   3,   "norm-form k=1"));
listput(results, check_splits(19,  5,   "norm-form k=3"));
listput(results, check_splits(37,  7,   "norm-form k=5"));
listput(results, check_splits(79,  9,   "norm-form k=7 CM-73"));
listput(results, check_splits(109, 11,  "norm-form k=9 CM-73"));
print();

\\ Test group B: non-norm-form primes with various a2
print("-- Group B: non-norm-form primes, arbitrary a2 --");
listput(results, check_splits(11,  4,   "p=11 a2=4"));
listput(results, check_splits(23,  6,   "p=23 a2=6 (D>0 real quad)"));
listput(results, check_splits(41,  8,   "p=41 a2=8"));
listput(results, check_splits(53,  10,  "p=53 a2=10 h=12"));
listput(results, check_splits(47,  22,  "p=47 a2=22 large m/p"));
print();

\\ Test group C: large p (crypto-scale sanity)
print("-- Group C: larger primes --");
listput(results, check_splits(1009, 30,  "p=1009"));
listput(results, check_splits(2003, 44,  "p=2003"));
listput(results, check_splits(9973, 100, "p=9973"));
print();

\\ Tally
splits = sum(i=1, #results, results[i] == 1);
other  = #results - splits;
printf("=== Summary: %d/%d split (all expected) ===\n", splits, #results);
if(other == 0,
  print("ALL PASS: p splits in Q(sqrt(sf)) whenever p∤a2. Theorem (16) confirmed."),
  printf("FAIL: %d non-split cases — investigate!\n", other));
