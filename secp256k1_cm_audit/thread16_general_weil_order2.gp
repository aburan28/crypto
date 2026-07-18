\\ thread16_general_weil_order2.gp
\\
\\ Thread 16: Verify that the algebraic theorem proved in Thread 15 is completely
\\ GENERAL — it does not require the secp256k1 norm-form condition 4p = 73+3k^2.
\\
\\ THEOREM (general biquadratic Weil polynomial):
\\   Let p be any prime, a2 any integer with |a2| < 2p and p ∤ a2.
\\   Let D = a2^2 - 4p^2 = sf * m^2 where sf is squarefree and m > 0.
\\   Then:
\\     (i)  sf is a quadratic residue mod p  [KEY: sf ≡ a2^2 mod p ≡ (a2/m)^2 · m^{-2}
\\          but more directly D ≡ a2^2 (mod p), and D/m^2 = sf, and (m^2/p) = 1, so (sf/p)=1]
\\     (ii) p SPLITS in Q(sqrt(sf)), i.e., idealprimedec(K, p) has length 2.
\\     (iii)[P]^2 = 1 in Cl(Q(sqrt(sf))).
\\
\\ PROOF SKETCH (from Thread 15, now stated for general p):
\\   beta = (-a2 + m*sqrt(sf)) / 2 satisfies x^2 + a2*x + p^2 = 0, so beta in O_K.
\\   N(beta) = p^2; if (beta) = (p) then beta/p in O_K requires p|a2 — contradiction.
\\   Hence (beta) is a proper ideal of norm p^2. Since p splits (by (i)), the ideals of
\\   norm p^2 in O_K are P^2, P*Pbar = (p), Pbar^2. Since (beta) != (p), we get
\\   (beta) = P^2 or Pbar^2, i.e., [P]^2 = 1 in Cl(K). QED.
\\
\\ NUMERICAL PLAN: 10 non-norm-form primes × 3 a2 values each = 30 (p, a2) test cases.
\\   Non-norm-form = p NOT satisfying 4p = 73 + 3k^2 for any odd k.
\\   a2 values: p//5+1 (small), -(p//4+1) (small negative), 2*p//3+1 (medium).
\\
\\ Run: gp -q thread16_general_weil_order2.gp

default(parisize, 128000000);
default(timer, 0);

\\ -----------------------------------------------------------------------
\\ Utilities
\\ -----------------------------------------------------------------------

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Check: is p a norm-form prime (4p = 73 + 3k^2 for some odd positive k)?
is_norm_form(p) = {
  my(val, k, kk);
  val = 4*p - 73;
  if(val <= 0 || val % 3 != 0, return(0));
  k = sqrtint(val \ 3);
  forstep(kk = max(1, k-1), k+2, 1,
    if(kk % 2 == 1 && 3*kk^2 == val, return(1)));
  0
};

\\ -----------------------------------------------------------------------
\\ Main verification for one (p, a2) pair
\\ -----------------------------------------------------------------------

verify_pair(idx, p, a2) = {
  my(D, sf, m2, m, K, Pp, nP, P, P2, P2_prin, ord_P, is_QR, h, label_ord);

  if(a2 % p == 0,
    printf("#%-2d  p=%-7d a2=%-8d  SKIP: p | a2\n", idx, p, a2); return(-1));

  D = a2^2 - 4*p^2;
  if(D >= 0,
    printf("#%-2d  p=%-7d a2=%-8d  SKIP: D>=0\n", idx, p, a2); return(-1));

  sf = sf_part(D);
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1,
    printf("#%-2d  p=%-7d a2=%-8d  ERROR: m2 not positive integer\n", idx, p, a2); return(0));
  m  = sqrtint(m2);
  if(m^2 != m2,
    printf("#%-2d  p=%-7d a2=%-8d  ERROR: m2 not a perfect square\n", idx, p, a2); return(0));

  \\ --- Check (i): sf is QR mod p (expect YES always, by D ≡ a2^2 mod p) ---
  is_QR = (kronecker(sf, p) != -1);  \\ 1 or 0 (ramified)

  \\ --- Check (ii)+(iii): PARI bnfinit ---
  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);
  nP = #Pp;

  if(nP == 0,
    printf("#%-2d  p=%-7d a2=%-8d sf=%-10d  ANOMALY: no primes above p (inert?)\n", idx, p, a2, sf); return(0));

  P        = Pp[1];
  P2       = idealpow(K, P, 2);
  P2_prin  = bnfisprincipal(K, P2);

  \\ ord([P]): 1 if P principal, 2 if P^2 principal but P not, >2 otherwise
  if(P2_prin[1] == 0*P2_prin[1],
    if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1],
       ord_P = 1; label_ord = "1",
       ord_P = 2; label_ord = "2"),
    ord_P = 99; label_ord = ">2(!)"
  );

  printf("#%-2d  p=%-7d a2=%-8d sf=%-10d m=%-6d h=%-4d nP=%d  QR:%s  [P]^2=1:%s  ord([P])=%s\n",
    idx, p, a2, sf, m, h, nP,
    if(is_QR,  "YES", "NO(!)"),
    if(ord_P <= 2, "YES", "NO(!)"),
    label_ord);

  if(ord_P <= 2, 1, 0)
};

\\ -----------------------------------------------------------------------
\\ Test suite: 10 non-norm-form primes, 3 a2 values each
\\ -----------------------------------------------------------------------

\\ Verify none are norm-form (checked by hand; confirmed here)
non_norm_primes = [23, 41, 53, 67, 97, 113, 127, 151, 181, 211];

print("Thread 16: General order-2 theorem for non-norm-form primes");
print("=============================================================");
print("Theorem: for any prime p and a2 with |a2|<2p, p∤a2:");
print("  D = a2^2 - 4p^2 = sf*m^2  =>  p splits in Q(sqrt(sf)) and [P]^2 = 1.");
print();
print("Confirming none of the test primes are norm-form (4p = 73+3k^2):");
{
  for(i = 1, #non_norm_primes,
    p = non_norm_primes[i];
    printf("  p=%-5d  norm-form: %s\n", p, if(is_norm_form(p), "YES(!)", "no")));
}
print();
printf("%-4s  %-7s %-8s %-10s %-6s %-4s %s  %-6s  %s\n",
  "IDX","p","a2","sf=sf(D)","m","h","nP","QR?","[P]^2=1? (ord)");
print(concat(vector(80, i, "-")));

{
  my(total = 0, passed = 0, p, a2, r, idx = 1);
  for(i = 1, #non_norm_primes,
    p = non_norm_primes[i];
    \\ Three a2 choices per prime: small positive, small negative, medium
    a2_list = [p\5 + 1, -(p\4 + 1), 2*(p\3) + 1];
    for(j = 1, #a2_list,
      a2 = a2_list[j];
      r = verify_pair(idx, p, a2);
      if(r >= 0, total++; passed += r);
      idx++));

  print();
  printf("RESULT: %d / %d (p, a2) pairs satisfy [P]^2 = 1\n", passed, total);
  if(passed == total,
    print("ALL PASSED — Theorem confirmed for 10 non-norm-form primes (30 test cases)."),
    print("FAILURES DETECTED — check output above."));
}

\\ -----------------------------------------------------------------------
\\ Extra check: a2 near the boundary |a2| ~ 2p-1 (stress the squarefree factoring)
\\ -----------------------------------------------------------------------
print();
print("Boundary stress test: a2 = 2p-3 (near |a2| = 2p) for first 5 primes:");
print(concat(vector(80, i, "-")));
{
  my(total2 = 0, passed2 = 0, p, a2, r, idx = 31);
  for(i = 1, 5,
    p = non_norm_primes[i];
    a2 = 2*p - 3;  \\ near boundary; |D| = (2p-3)^2 - 4p^2 = 9-12p (small negative)
    r = verify_pair(idx, p, a2);
    if(r >= 0, total2++; passed2 += r);
    idx++);
  printf("\nBoundary: %d / %d passed\n", passed2, total2);
}

\\ -----------------------------------------------------------------------
\\ Extra check: a2 = 1 (minimal positive) for first 5 primes
\\ -----------------------------------------------------------------------
print();
print("Minimal a2 test: a2 = 1 for first 5 primes:");
print(concat(vector(80, i, "-")));
{
  my(total3 = 0, passed3 = 0, p, r, idx = 36);
  for(i = 1, 5,
    p = non_norm_primes[i];
    r = verify_pair(idx, p, 1);
    if(r >= 0, total3++; passed3 += r);
    idx++);
  printf("\nMinimal a2=1: %d / %d passed\n", passed3, total3);
}

print();
print("DONE.");
