\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify that the order-2 Frobenius theorem generalises beyond
\\ the secp256k1 norm-form family 4p = 73 + 3k^2.
\\
\\ THEOREM (Thread 15, generalised statement):
\\   Let p be any prime and t any integer with 0 < |t| < 2*sqrt(p) and p ∤ t.
\\   Form the product Jacobian J = E × E^t (E: elliptic curve over F_p with trace t,
\\   E^t its quadratic twist with trace -t).  The Weil polynomial of J is
\\     F(T) = T^4 + a2*T^2 + p^2,  a2 = 2p - t^2.
\\   Let D = a2^2 - 4p^2 = t^2*(t^2 - 4p) < 0,
\\       sf = squarefree part of D  (so D = sf * m^2, m > 0),
\\       K  = Q(sqrt(sf))          (an imaginary quadratic field).
\\   Let P be any prime of O_K above p.  Then [P]^2 = 1 in Cl(K).
\\
\\ PROOF SKETCH (see thread15 for detailed (A)-(E)):
\\   beta = (-a2 + m*sqrt(sf)) / 2 satisfies x^2 + a2*x + p^2 = 0, so beta in O_K.
\\   N_{K/Q}(beta) = p^2, so (beta) has norm p^2 in O_K.
\\   Since p ∤ a2 (as |t| < 2*sqrt(p) and p ∤ t imply p ∤ t^2, hence p ∤ a2 = 2p-t^2 mod p),
\\   we have (beta) ≠ (p), so (beta) = P^2 or Pbar^2.  Hence [P]^2 = 1. QED.
\\
\\ EXPERIMENT: verify numerically for 10 non-norm-form primes.
\\ For each test pair (p, t) chosen so that p is NOT of the form (73+3k^2)/4:
\\   - Confirm F(T) = T^4 + a2*T^2 + p^2 with a2 = 2p - t^2
\\   - Compute D, sf, m; verify D = sf*m^2
\\   - Verify p ∤ a2
\\   - Verify [P]^2 = 1 in Cl(Q(sqrt(sf)))
\\
\\ Run: gp --stacksize 64000000 -q thread16_general_order2.gp

default(parisize, 64000000);
default(timer, 0);

\\ ---------------------------------------------------------------
\\ Utility: signed squarefree part
\\ ---------------------------------------------------------------
sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ---------------------------------------------------------------
\\ Check norm-form: is p = (73 + 3*k^2)/4 for some odd k?
\\ ---------------------------------------------------------------
is_normform(p) = {
  my(k);
  \\ 4p - 73 must be divisible by 3 and a perfect odd square
  if((4*p - 73) % 3 != 0, return(0));
  k = sqrtint((4*p - 73) \ 3);
  if(k*k == (4*p - 73) \ 3 && k % 2 == 1 && (4*p - 73) == 3*k*k, 1, 0)
};

\\ ---------------------------------------------------------------
\\ Core verification for a single (p, t) pair
\\ ---------------------------------------------------------------
\\ Returns 1 if [P]^2 = 1 (theorem holds), 0 otherwise.
\\ Prints a result line.
verify_general(p, t, label) = {
  my(a2, D, sf, m2, m, K, Pp, P, P2, h, beta_elt, I_beta,
     norm_beta, p_div_a2, P2_prin, ord_str, status);

  \\ Basic sanity
  if(!isprime(p),
    printf("  %-12s: p=%d not prime\n", label, p); return(-1));
  if(t == 0 || abs(t) >= 2*sqrt(p),
    printf("  %-12s: p=%d t=%d outside Hasse bound\n", label, p, t); return(-1));
  if(t % p == 0,
    printf("  %-12s: p=%d t=%d divisible by p — excluded\n", label, p, t); return(-1));

  \\ Biquadratic Weil polynomial setup
  a2 = 2*p - t^2;
  D  = a2^2 - 4*p^2;               \\ = t^2*(t^2 - 4p) < 0
  sf = sf_part(D);

  \\ Check D = sf*m^2 with m > 0
  m2 = D / sf;
  if(denominator(m2) != 1 || m2 <= 0,
    printf("  %-12s: p=%d D/sf=%Ps not a positive integer\n", label, p, m2);
    return(-1));
  m = sqrtint(m2);
  if(m*m != m2,
    printf("  %-12s: p=%d D/sf=%d not a perfect square\n", label, p, m2);
    return(-1));

  \\ Check (iii): p ∤ a2
  p_div_a2 = (a2 % p == 0);

  \\ Norm of beta
  norm_beta = (a2^2 - m2*sf) / 4;  \\ should equal p^2

  \\ PARI ideal computation in K = Q(sqrt(sf))
  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);

  if(#Pp == 0,
    printf("  %-12s: p=%d p inert in Q(sqrt(%d)) — unexpected\n", label, p, sf);
    return(-1));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);

  \\ Is P^2 principal?
  P2_prin = bnfisprincipal(K, P2);
  \\ P2_prin[1] is the exponent vector; all-zero <=> principal
  status = (matsize(P2_prin[1])[2] == 0 || norml2(P2_prin[1]) == 0);

  ord_str = if(status, "1",
    if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1], "1 (P prin)", ">2 ERROR"));

  printf("  %-12s p=%-7d t=%-5d a2=%-9d sf=%-8d m=%-6d h=%-4d | N(β)=p²:%s  p|a2:%s  [P]²=1:%s\n",
    label, p, t, a2, sf, m, h,
    if(norm_beta == p^2, "YES", "NO "),
    if(p_div_a2, "YES(!)", "no "),
    if(status, "YES", "NO "));

  if(status, 1, 0)
};

\\ ---------------------------------------------------------------
\\ Test suite: 10 non-norm-form primes with varied traces
\\ ---------------------------------------------------------------
\\ Primes chosen to be definitively NOT of the form (73+3k^2)/4.
\\ Traces chosen in the interior of the Hasse bound to avoid degenerate cases.
\\
\\ p=101:  4*101=404; 404-73=331; 331/3=110.33... not an integer => not norm-form
\\ p=127:  4*127=508; 508-73=435; 435/3=145; sqrt(145)~12.04 not integer => not norm-form
\\ p=157:  4*157=628; 628-73=555; 555/3=185; sqrt(185)~13.6 not integer => not norm-form
\\ p=199:  4*199=796; 796-73=723; 723/3=241; sqrt(241)~15.5 not integer => not norm-form
\\ p=251:  4*251=1004; 1004-73=931; 931/3=310.3... not integer => not norm-form
\\ p=307:  4*307=1228; 1228-73=1155; 1155/3=385; sqrt(385)~19.6 not integer => not norm-form
\\ p=401:  4*401=1604; 1604-73=1531; 1531/3=510.3... not integer => not norm-form
\\ p=503:  4*503=2012; 2012-73=1939; 1939/3=646.3... not integer => not norm-form
\\ p=601:  4*601=2404; 2404-73=2331; 2331/3=777; sqrt(777)~27.9 not integer => not norm-form
\\ p=701:  4*701=2804; 2804-73=2731; 2731/3=910.3... not integer => not norm-form

test_cases() = {
  my(cases, passed, total, ok);
  \\       [prime, trace, label]
  cases = [
    [101,  3,  "p=101,t=3"],
    [127,  5,  "p=127,t=5"],
    [157,  7,  "p=157,t=7"],
    [199,  9,  "p=199,t=9"],
    [251,  4,  "p=251,t=4"],
    [307, 11,  "p=307,t=11"],
    [401,  6,  "p=401,t=6"],
    [503, 13,  "p=503,t=13"],
    [601, 10,  "p=601,t=10"],
    [701, 15,  "p=701,t=15"]
  ];
  passed = 0; total = 0;
  for(i = 1, #cases,
    my(row = cases[i], p = row[1], t = row[2], lbl = row[3]);
    ok = verify_general(p, t, lbl);
    if(ok >= 0, total++; if(ok == 1, passed++)));
  printf("\nResult: %d / %d passed (theorem [P]^2=1 holds).\n", passed, total);
  [passed, total]
};

\\ ---------------------------------------------------------------
\\ Bonus: spot-check norm-form exclusion (sanity)
\\ ---------------------------------------------------------------
check_normform_exclusion() = {
  my(normform_count);
  normform_count = 0;
  print("\nNorm-form exclusion check (all 10 test primes should be non-norm-form):");
  for(i = 1, 10,
    my(p_list = [101, 127, 157, 199, 251, 307, 401, 503, 601, 701],
       p = p_list[i], nf = is_normform(p));
    if(nf, normform_count++; printf("  p=%d: NORM-FORM (!!)\n", p),
            printf("  p=%d: not norm-form ✓\n", p)));
  printf("  => %d / 10 accidentally norm-form.\n", normform_count);
  normform_count
};

\\ ---------------------------------------------------------------
\\ Statement of the general theorem (printed header)
\\ ---------------------------------------------------------------
print("Thread 16: General order-2 Frobenius theorem (non-norm-form verification)");
print("==========================================================================");
print();
print("THEOREM (generalised): For any prime p and trace t with 0 < |t| < 2*sqrt(p)");
print("and p ∤ t, the prime ideal P of O_{Q(sqrt(sf))} above p (where sf = sf-part(D),");
print("D = t^2*(t^2-4p)) satisfies [P]^2 = 1 in Cl(Q(sqrt(sf))).");
print();
print("The algebraic proof (A)-(E) from thread15 uses only:");
print("  (A) beta = (-a2 + m*sqrt(sf))/2 satisfies x^2 + a2*x + p^2 = 0");
print("  (B) beta lies in K = Q(sqrt(sf))");
print("  (C) N(beta) = p^2");
print("  (D) p ∤ a2  (since a2 = 2p-t^2, so a2 ≡ -t^2 mod p, and p ∤ t)");
print("  (E) Therefore (beta) = P^2 or Pbar^2 => [P]^2 = 1.");
print("None of these steps use the norm-form condition 4p = 73+3k^2.");
print();
print("Numerical verification for 10 non-norm-form primes:");
print("--------------------------------------------------------------------");

\\ Run tests
check_normform_exclusion();
print();
print("Main verification:");
test_cases();

print();
print("DONE.");
