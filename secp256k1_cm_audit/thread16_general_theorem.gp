\\ thread16_general_theorem.gp
\\
\\ Thread 16: Verify the universal order-2 Frobenius theorem is GENERAL —
\\ it holds for ANY biquadratic Weil polynomial T^4 + a2*T^2 + p^2, not
\\ only the secp256k1 norm-form family 4p = 73 + 3k^2.
\\
\\ THEOREM (general form, proved algebraically in Thread 15):
\\   Let p be a prime, a2 an integer with |a2| < 2p, and set
\\     D = a2^2 - 4p^2 = sf * m^2   (sf squarefree, m in Z>0).
\\   Let K = Q(sqrt(sf)) and P a prime ideal above p in O_K.
\\   If p does not divide a2, then [P]^2 = 1 in Cl(K).
\\
\\ PROOF RECAP (purely algebraic, no norm-form needed):
\\   beta = (-a2 + m*sqrt(sf)) / 2  satisfies  x^2 + a2*x + p^2 = 0.
\\   (A) Monic, integer coefficients => algebraic integer => beta in O_K.
\\   (B) N(beta) = p^2 => (beta) has norm p^2 in O_K.
\\   (C) Ideals of norm p^2 in O_K: P^2, Pbar^2, or (p).
\\   (D) p nmid a2 => beta/p not in O_K => (beta) != (p).
\\   (E) Hence (beta) = P^2 or Pbar^2, so [P]^2 = 1.   QED.
\\
\\ TEST STRATEGY:
\\   Generate (p, a2) pairs from two sources:
\\     (S1) Product curves: E x E^twist over F_p, where E has trace a1,
\\          E^twist has trace -a1, so characteristic poly = T^4+(2p-a1^2)T^2+p^2.
\\          These are explicitly realizable; a2 = 2p - a1^2.
\\     (S2) Synthetic: arbitrary a2 with |a2| < 2p, p prime, p nmid a2.
\\          No actual curve needed — theorem is purely algebraic.
\\
\\ For each case: compute D, sf, m, K = Q(sqrt(sf)), P above p.
\\   Check (i) m^2 * sf = D, (ii) p nmid a2, (iii) P^2 principal in O_K.
\\
\\ Run: gp --stacksize 256000000 -q thread16_general_theorem.gp

default(parisize, 256000000);
default(timer, 0);

\\ ---------------------------------------------------------------
\\ Core verification function
\\ ---------------------------------------------------------------

{verify_biquadratic(p, a2, label) =
  my(D, sf, m, K, P_list, P, P2, cl, is_principal, h_K, ord_P);

  \\ --- Basic sanity ---
  if (!isprime(p), error("p not prime: ", p));
  if (abs(a2) >= 2*p, printf("%s: SKIP (|a2|>=2p)\n", label); return(0));
  if (a2 == 0, printf("%s: SKIP (a2=0, degenerate)\n", label); return(0));

  \\ D = a2^2 - 4p^2 must be < 0 (imaginary quadratic)
  D = a2^2 - 4*p^2;
  if (D >= 0, printf("%s: SKIP (D>=0)\n", label); return(0));

  \\ Squarefree decomposition: D = sf * m^2
  sf = core(D);   \\ core() gives squarefree kernel, preserving sign
  if (D % sf != 0 || !issquare(D / sf),
    printf("%s: ERROR D/sf not perfect square\n", label); return(0));
  m = sqrtint(abs(D / sf));
  if (sf * m^2 != D,
    printf("%s: ERROR D != sf*m^2 (D=%d, sf=%d, m=%d)\n", label, D, sf, m);
    return(0));

  \\ Condition (D): p must not divide a2
  if (a2 % p == 0,
    printf("%s: SKIP (p | a2, exceptional case)\n", label); return(0));

  \\ Construct K = Q(sqrt(sf))
  K = bnfinit(x^2 - sf, 1);
  h_K = K.no;

  \\ Trivial case: h_K = 1 => every ideal is principal => [P]^2 = 1 automatically
  if (h_K == 1,
    printf("%s: PASS (h_K=1, trivially principal) p=%d a2=%d sf=%d m=%d\n",
           label, p, a2, sf, m);
    return(1));

  \\ Find prime ideals above p in K
  P_list = idealprimedec(K, p);
  if (#P_list == 0,
    printf("%s: ERROR no primes above p in K\n", label); return(0));
  P = P_list[1];

  \\ Compute P^2 and check principality
  P2 = idealpow(K, P, 2);
  cl = bnfisprincipal(K, P2);

  \\ cl[1] is the exponent vector in terms of class group generators.
  \\ It is zero iff P^2 is principal.
  is_principal = (norml2(cl[1]) == 0);

  \\ Compute the order of [P] in Cl(K) using ideal order (brute force, small h_K)
  ord_P = 0;
  if (h_K <= 200,
    my(Pcur = P);
    for (j = 1, h_K,
      my(clj = bnfisprincipal(K, Pcur));
      if (norml2(clj[1]) == 0, ord_P = j; break);
      Pcur = idealmul(K, Pcur, P)));

  printf("%s: %s  p=%d a2=%d sf=%d m=%d h_K=%d ord([P])=%d\n",
         label, if(is_principal, "PASS", "FAIL !!!"),
         p, a2, sf, m, h_K, ord_P);

  return(is_principal);
}

\\ ---------------------------------------------------------------
\\ Test cases (S1): product-curve pairs  (a2 = 2p - a1^2)
\\ NOT from the secp256k1 norm-form family 4p = 73 + 3k^2.
\\ Chosen: random small/medium primes with diverse a1 values.
\\ ---------------------------------------------------------------

print("=== Thread 16: General Biquadratic Theorem Verification ===");
print("=== Source S1: product-curve cases (a2 = 2p - a1^2) ===");
print("");

{
  my(cases_s1, p, a1, a2, pass_count, total);
  pass_count = 0; total = 0;

  \\ [p, a1] pairs.  a2 = 2p - a1^2.
  \\ All primes below are NOT of the form (73 + 3k^2)/4.
  cases_s1 = [
    [101,   5],   \\ a2 = 202-25   = 177
    [101,  10],   \\ a2 = 202-100  =  102  --> but 101|a2? 102=101+1, no
    [107,   7],   \\ a2 = 214-49   = 165
    [127,   9],   \\ a2 = 254-81   = 173
    [199,  13],   \\ a2 = 398-169  = 229
    [311,  17],   \\ a2 = 622-289  = 333
    [499,  21],   \\ a2 = 998-441  = 557
    [997,  31],   \\ a2 = 1994-961 = 1033
    [2003, 43],   \\ a2 = 4006-1849= 2157
    [4999, 67],   \\ a2 = 9998-4489= 5509
    [9973, 97],   \\ a2 = 19946-9409=10537
    [49999, 223], \\ a2 = 99998-49729=50269
    [99991, 311], \\ a2 = 199982-96721=103261
    [100003, 317],\\ a2 = 200006-100489=99517
    [999983, 997] \\ a2 = 1999966-994009=1005957
  ];

  for (i = 1, #cases_s1,
    p = cases_s1[i][1];
    a1 = cases_s1[i][2];
    a2 = 2*p - a1^2;
    total++;
    if (verify_biquadratic(p, a2,
          Str("S1[", i, "] p=", p, " a1=", a1)),
      pass_count++));

  printf("\nS1 summary: %d/%d passed\n\n", pass_count, total);
}

\\ ---------------------------------------------------------------
\\ Test cases (S2): synthetic a2 values (not from product curves)
\\ These test the purely algebraic nature of the theorem.
\\ ---------------------------------------------------------------

print("=== Source S2: synthetic a2 values (purely algebraic) ===");
print("");

{
  my(cases_s2, p, a2, pass_count, total);
  pass_count = 0; total = 0;

  \\ [p, a2] pairs chosen freely with |a2| < 2p, p prime, p nmid a2.
  \\ These need not correspond to any actual curve or product.
  cases_s2 = [
    [53,    47],   \\ random: D = 47^2 - 4*53^2 = 2209 - 11236 = -9027
    [71,    53],   \\ D = 2809 - 20164 = -17355
    [83,    61],   \\ D = 3721 - 27556 = -23835
    [101,   79],   \\ D = 6241 - 40804 = -34563
    [113,   89],   \\ D = 7921 - 51076 = -43155
    [131,   103],  \\ D = 10609 - 68644 = -58035
    [157,   127],  \\ D = 16129 - 98596 = -82467
    [173,   137],  \\ D = 18769 - 119716= -100947
    [179,   151],  \\ D = 22801 - 128164= -105363
    [223,   199],  \\ D = 39601 - 198916= -159315
    [1009,  887],  \\ D = 786769 - 4072324= -3285555
    [10007, 9001], \\ D = 81018001 - 400280196 = -319262195
    [100003, 98765],\\ large p, a2 close to 2p
    [999983, 999001]\\ very large
  ];

  for (i = 1, #cases_s2,
    p = cases_s2[i][1];
    a2 = cases_s2[i][2];
    total++;
    if (verify_biquadratic(p, a2,
          Str("S2[", i, "] p=", p, " a2=", a2)),
      pass_count++));

  printf("\nS2 summary: %d/%d passed\n\n", pass_count, total);
}

\\ ---------------------------------------------------------------
\\ Test cases (S3): edge cases — negative a2, large class numbers
\\ ---------------------------------------------------------------

print("=== Source S3: edge cases (negative a2, a2 near 0, large h_K) ===");
print("");

{
  my(cases_s3, p, a2, pass_count, total);
  pass_count = 0; total = 0;

  cases_s3 = [
    [101, -79],    \\ negative a2: same D as S2[4] but different sign => same sf
    [101, -5],     \\ small negative a2
    [127, -9],     \\ negative version of S1[4]
    [257, 1],      \\ a2 very small positive
    [257, -1],     \\ a2 = -1
    [251, 3],      \\ a2 = 3 (tiny)
    [1013, 7],     \\ a2 = 7, large p => D very negative
    [10007, 11],   \\ a2 = 11
    [100003, 13]   \\ a2 = 13, very large p => expect large h_K
  ];

  for (i = 1, #cases_s3,
    p = cases_s3[i][1];
    a2 = cases_s3[i][2];
    total++;
    if (verify_biquadratic(p, a2,
          Str("S3[", i, "] p=", p, " a2=", a2)),
      pass_count++));

  printf("\nS3 summary: %d/%d passed\n\n", pass_count, total);
}

\\ ---------------------------------------------------------------
\\ Comparison: re-verify 5 known norm-form cases as control
\\ ---------------------------------------------------------------

print("=== Control: 5 known norm-form cases from Thread 15 ===");
print("");

{
  my(norm_form, p, k, a2_nf, pass_count, total);
  pass_count = 0; total = 0;

  \\ (k, p) from norm-form family: 4p = 73 + 3k^2
  \\ a2 from the isogeny-graph SPLIT condition (stored in Thread15 results)
  norm_form = [
    [1,   19,   35],   \\ k=1, p=19, a2=35 (from Thread15 table)
    [5,   37,   -1],   \\ k=5, p=37, a2=-1
    [9,   79,   -85],  \\ k=9, p=79, a2=-85
    [11,  109,  -145], \\ k=11, p=109, a2=-145
    [21,  349,  -385]  \\ k=21, p=349, a2=-385
  ];

  for (i = 1, #norm_form,
    k = norm_form[i][1];
    p = norm_form[i][2];
    a2_nf = norm_form[i][3];
    total++;
    if (verify_biquadratic(p, a2_nf,
          Str("NF k=", k, " p=", p)),
      pass_count++));

  printf("\nControl summary: %d/%d passed\n\n", pass_count, total);
}

\\ ---------------------------------------------------------------
\\ Grand total
\\ ---------------------------------------------------------------

print("=== Done. All PASS = theorem confirmed general. ===");
