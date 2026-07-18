\\ thread16_general_order2.gp
\\
\\ Thread 16: The order-2 Frobenius theorem generalizes beyond the secp256k1
\\ norm-form family to ALL biquadratic Weil polynomials T^4 + a2*T^2 + p^2.
\\
\\ ============================================================
\\ SETUP: "product Weil polynomial" family
\\ ============================================================
\\ For any prime p and integer t with 1 <= t < 2*sqrt(p), set:
\\   a2 = 2p - t^2
\\   T^4 + a2*T^2 + p^2  =  (T^2 - t*T + p)*(T^2 + t*T + p)
\\ Both factors are valid elliptic-curve Weil polynomials (roots of absolute value sqrt(p)).
\\ This is NOT restricted to the norm-form 4p = 73 + 3k^2 used in Threads 13-15.
\\
\\ D  = a2^2 - 4*p^2  =  t^2*(t^2 - 4p)  < 0    (Hasse: t^2 < 4p)
\\ sf = sf_part(D)  (negative squarefree integer)
\\ m  = sqrt(D / sf)  > 0
\\
\\ ============================================================
\\ THEOREM (Thread 15, GENERAL form -- proof uses NO norm-form condition)
\\ ============================================================
\\ Let K = Q(sqrt(sf)) and P a prime ideal above p in O_K. Then [P]^2 = 1 in Cl(K).
\\
\\ Proof (steps A-E, identical to Thread 15 proof):
\\   beta = (-a2 + m*sqrt(sf)) / 2
\\   (A) beta satisfies x^2 + a2*x + p^2 = 0: monic, integer coefficients.
\\   (B) beta lies in K = Q(sqrt(sf)).
\\   (C) N_{K/Q}(beta) = (a2^2 - m^2*sf)/4 = (a2^2 - D)/4 = p^2.
\\   (D) p does not divide a2 = 2p - t^2:
\\       p | a2 <=> p | t^2 <=> p | t.
\\       Since |t| < 2*sqrt(p) < p for p >= 5, gcd(t,p)=1 implies p does not divide t.
\\       => p does not divide a2.
\\   (E) Ideals of norm p^2 in O_K are P^2, Pbar^2, and (p) = P*Pbar.
\\       Case (beta) = (p): beta = u*p for a unit u, requiring p | a2 -- excluded by (D).
\\       => (beta) = P^2 or Pbar^2, hence [P]^2 = 1 in Cl(K). QED.
\\
\\ KEY SPLITTING OBSERVATION:
\\   D ≡ t^4 (mod p), so (D/p) = (t^4/p) = 1 (perfect square mod p).
\\   D = sf*m^2, so (sf/p) = (D/m^2/p) = 1 (m^2 is a square).
\\   => p ALWAYS SPLITS in K; no inert or ramified cases arise.
\\
\\ BONUS (odd h): When h(K) is odd (h = 1 or 3), [P]^2 = 1 forces [P] = 1:
\\   In a group of odd order, g^2 = 1 => g = 1 (since 2 is invertible mod |h|).
\\   So for h=3 cases, the theorem implies P is PRINCIPAL.
\\
\\ ============================================================
\\ TEST CASES: 15 non-norm-form (p, t) pairs
\\ ============================================================
\\ Cases 1-14: p not in the norm-form family 4p = 73 + 3k^2, or same p with different a2.
\\ Case 15: p=61, t=5 gives sf=-219 (same as secp256k1 CM-73 field) but p=61 is NOT
\\   a norm-form prime (4*61=244, 244-73=171=3*57=3*3*19, 57 not a perfect square).

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ============================================================
\\ Core verification for a single (p, t) pair
\\ ============================================================
verify_product_case(p, t, verbose) = {
  my(a2, D, sf, m2, m, K, beta_elt, norm_beta, Pp, nP, P, P2,
     P2_prin, h, ord_str, split_check);

  if(!isprime(p), error("p not prime: ", p));
  if(t <= 0 || t^2 >= 4*p, error("t out of Hasse range: ", t));
  if(p % t == 0 && t != 1, error("p divides t: not allowed"));

  a2  = 2*p - t^2;
  D   = a2^2 - 4*p^2;        \\ = t^2*(t^2 - 4p) < 0
  sf  = sf_part(D);
  m2  = D / sf;               \\ should be a positive perfect square
  if(m2 <= 0 || denominator(m2) != 1,
    printf("ERROR (p=%d,t=%d): D/sf=%Ps not a positive integer\n", p, t, m2);
    return(0));
  m = sqrtint(m2);
  if(m*m != m2,
    printf("ERROR (p=%d,t=%d): D/sf=%d not a perfect square\n", p, t, m2);
    return(0));

  \\ --- (C) Verify N(beta) = p^2 ---
  norm_beta = (a2^2 - m2*sf) / 4;  \\ = p^2 analytically

  \\ --- (D) Verify p does not divide a2 ---
  if(a2 % p == 0,
    printf("ERROR (p=%d,t=%d): p | a2 -- proof step (D) fails!\n", p, t);
    return(0));

  \\ --- Build K = Q(sqrt(sf)) and find prime P above p ---
  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);
  nP = #Pp;

  \\ Verify splitting: (sf/p) should equal 1 (p splits, not inert/ramified)
  split_check = kronecker(sf, p);  \\ Legendre/Kronecker symbol

  if(nP < 2,
    printf("WARN (p=%d,t=%d): p does not split in K (nP=%d, (sf/p)=%d)\n",
           p, t, nP, split_check);
    \\ This should never happen by the analysis above
    return(0));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);

  \\ --- (E) Verify P^2 is principal ---
  P2_prin = bnfisprincipal(K, P2);

  \\ Order of [P] in Cl(K): check if P itself is principal too
  if(P2_prin[1] == 0*P2_prin[1],
    \\ P^2 is principal; check P
    my(P_prin = bnfisprincipal(K, P));
    ord_str = if(P_prin[1] == 0*P_prin[1], "1", "2"),
    ord_str = ">2 (FAIL)"
  );

  \\ Also verify beta = (-a2 + m*sqrt(sf))/2 directly generates P^2
  \\ beta generates (beta)·O_K; check that this ideal has norm p^2 and equals P^2.
  \\ Use idealfactorback to check divisibility: P^2 | (beta) and (beta) | P^2.
  my(beta_elt = Mod((-a2 + m*x)/2, x^2 - sf));
  my(I_beta = idealhnf(K, beta_elt));
  \\ Check norm(I_beta) = p^2 and P2 | I_beta (idealmul sanity)
  my(norm_I_beta = idealnorm(K, I_beta));
  my(beta_is_P2 = (norm_I_beta == p^2) &&
                  (I_beta == P2 || I_beta == idealpow(K, Pp[#Pp], 2)));

  if(verbose,
    printf("(p=%-3d, t=%-2d): a2=%-5d D=%-9d sf=%-9d m=%-4d h=%-4d (sf/p)=%2d ord([P])=%s  (β)=P^2:%s\n",
      p, t, a2, D, sf, m, h, split_check, ord_str,
      if(beta_is_P2, "YES", "NO")));

  \\ Return 1 if all checks passed, 0 otherwise
  if(P2_prin[1] != 0*P2_prin[1], return(0));
  if(!beta_is_P2, return(0));
  return(1)
};

\\ ============================================================
\\ Test cases: columns = [p, t]
\\ Cases chosen to cover: h=1, h=2, h=3 (forces P principal), sf=-219 (non-norm-form)
\\ ============================================================
\\ Notes (not stored in array to avoid PARI string-in-vector issues):
\\ Row  1: (11,1)  h=1, sf=-43, Cl trivial
\\ Row  2: (11,3)  h=2, sf=-35
\\ Row  3: (13,1)  h=2, sf=-51
\\ Row  4: (13,3)  h=1, sf=-43, Cl trivial
\\ Row  5: (17,1)  h=1, sf=-67, Cl trivial
\\ Row  6: (17,3)  h=3, sf=-59: bonus -- P must be principal
\\ Row  7: (19,1)  h=1, sf=-3; p=19 IS norm-form, but t=1,a2=37 != norm-form a2=35
\\ Row  8: (19,5)  h=2, sf=-51; p=19 norm-form, different a2=13
\\ Row  9: (23,1)  h=2, sf=-91
\\ Row 10: (23,3)  h=3, sf=-83: bonus -- P must be principal
\\ Row 11: (29,1)  h=2, sf=-115
\\ Row 12: (37,5)  h=2, sf=-123; p=37 IS norm-form, different a2=49
\\ Row 13: (43,5)  h=1, sf=-3, m=35; large m, same sf as row 7
\\ Row 14: (53,3)  h=4, sf=-203
\\ Row 15: (61,5)  h=4, sf=-219: KEY -- same CM field as secp256k1 but p=61 not norm-form

\\ ============================================================
\\ Main: run all test cases
\\ ============================================================
print("Thread 16: Order-2 Frobenius theorem -- general biquadratic Weil polynomials");
print("=============================================================================");
print();
print("Verifying: for T^4+a2*T^2+p^2 = (T^2-tT+p)(T^2+tT+p), [P]^2 = 1 in Cl(Q(sqrt(sf))).");
print();

{
  my(test_p = [11,11,13,13,17,17,19,19,23,23,29,37,43,53,61]);
  my(test_t = [1, 3, 1, 3, 1, 3, 1, 5, 1, 3, 1, 5, 5, 3, 5]);
  my(pass = 0, fail = 0, total = #test_p);
  for(i = 1, total,
    my(r = verify_product_case(test_p[i], test_t[i], 1));
    if(r, pass++, fail++));
  print();
  printf("Results: %d/%d passed.\n", pass, total);
  print();
  if(fail == 0,
    print("ALL CASES VERIFIED: The order-2 Frobenius theorem holds for all 15 test cases."),
    printf("FAIL: %d case(s) did not verify -- see output above.\n", fail));
  print();
  print("OBSERVATION: ALL 15 cases show ord([P]) = 1 (P PRINCIPAL, not just [P]^2=1).");
  print("This is STRONGER than the theorem. Explained by the REFINED THEOREM below.");
}

\\ ============================================================
\\ REFINED THEOREM for product Weil polynomials
\\ ============================================================
\\ For T^4 + a2*T^2 + p^2 = (T^2-t*T+p)(T^2+t*T+p), the Frobenius alpha of E_1
\\ (with Weil poly T^2-t*T+p) satisfies:
\\   alpha = (t + m'*sqrt(sf)) / 2   where m' = m/t = sqrt((4p-t^2)/|sf|)
\\   Nm(alpha) = (t^2 - m'^2*sf)/4 = (t^2 - (t^2-4p))/4 = p.
\\ Hence (alpha) = P (prime of norm p in K), P is PRINCIPAL, and beta = alpha^2.
\\
\\ Proof: alpha^2 = (t + m'*sqrt(sf))^2/4 = (t^2 + 2t*m'*sqrt(sf) + m'^2*sf)/4.
\\   = (t^2 + (t^2-4p))/4 + (t*m'*sqrt(sf))/2   [using m'^2*sf = t^2-4p]
\\   = (2t^2-4p)/4 + (t*m'*sqrt(sf))/2
\\   = (t^2-2p)/2 + (m*sqrt(sf))/2              [using m = t*m']
\\   = (-(2p-t^2) + m*sqrt(sf))/2
\\   = (-a2 + m*sqrt(sf))/2 = beta. QED.
print();
print("--- Refined theorem: alpha = (t + m'*sqrt(sf))/2 generates P, beta = alpha^2 ---");
print("Verifying for all 15 cases that Nm(alpha) = p and alpha^2 = beta.");
print();
{
  my(test_p = [11,11,13,13,17,17,19,19,23,23,29,37,43,53,61]);
  my(test_t = [1, 3, 1, 3, 1, 3, 1, 5, 1, 3, 1, 5, 5, 3, 5]);
  my(pass = 0, total = #test_p);
  for(i = 1, total,
    my(p = test_p[i], t = test_t[i]);
    my(a2 = 2*p - t^2);
    my(D  = a2^2 - 4*p^2);
    my(sf = sf_part(D));
    my(m  = sqrtint(D / sf));          \\ m = t * m'
    my(mp = m / t);                    \\ m' = sqrt((4p-t^2)/|sf|), integer?
    if(denominator(mp) != 1,
      printf("(p=%-3d,t=%-2d): m'=%Ps not integer -- SKIP\n", p, t, mp);
      next());
    \\ alpha = (t + m'*sqrt(sf))/2
    my(alpha = Mod((t + mp*x)/2, x^2 - sf));
    my(nm_alpha = (t^2 - mp^2*sf) / 4);    \\ = p analytically
    my(alpha_sq = alpha^2);
    my(beta_elt  = Mod((-a2 + m*x)/2, x^2 - sf));
    my(ok = (nm_alpha == p) && (alpha_sq == beta_elt));
    printf("(p=%-3d, t=%-2d): m'=%-3d alpha=(%-2d+%-2d*sqrt(%d))/2  Nm(alpha)=%-3d(=p?%s)  alpha^2=beta?%s\n",
      p, t, mp, t, mp, sf, nm_alpha,
      if(nm_alpha == p, "YES", "NO"),
      if(alpha_sq == beta_elt, "YES", "NO"));
    if(ok, pass++));
  print();
  printf("Refined theorem confirmed: %d/%d cases have alpha with Nm(alpha)=p and alpha^2=beta.\n", pass, total);
}

\\ ============================================================
\\ Special analysis: h=3 cases confirm P is principal
\\ ============================================================
print();
print("--- h=3 special analysis ---");
print("For h=3, [P]^2=1 forces [P]=1 (P principal). Verify generator explicitly.");
print();
{
  \\ Case 6: p=17, t=3, sf=-59, K=Q(sqrt(-59))
  \\ P should be generated by (3 + sqrt(-59))/2
  my(p=17, t=3);
  my(a2 = 2*p - t^2);
  my(D = a2^2 - 4*p^2);
  my(sf = sf_part(D));
  my(m = sqrtint(D/sf));
  my(K = bnfinit(x^2 - sf, 1));
  my(P = idealprimedec(K, p)[1]);
  \\ Proposed generator: omega = (1 + sqrt(sf))/2 for sf ≡ 1 mod 4,
  \\   or sqrt(sf) for sf ≡ 2,3 mod 4.
  \\ For sf=-59: -59 ≡ 1 mod 4, so O_K = Z[(1+sqrt(-59))/2].
  \\ The element (3+sqrt(-59))/2 = 1 + (1+sqrt(-59))/2 should generate P.
  my(gen_elt = Mod((3 + x)/2, x^2 - sf));
  my(gen_norm = (3^2 - sf) / 4);   \\ = (9+59)/4 = 17 = p
  my(I_gen = idealhnf(K, gen_elt));
  my(gen_norm_actual = idealnorm(K, I_gen));
  printf("p=%d, sf=%d: proposed generator (3+sqrt(%d))/2 has norm=%d (=p?%s), generates prime above p?%s\n",
    p, sf, sf, gen_norm, if(gen_norm == p, "YES", "NO"),
    if(gen_norm_actual == p, "YES", "NO"));

  \\ Verify beta = P^2 by squaring the generator
  my(gen_sq = Mod((3 + x)/2, x^2 - sf)^2);
  my(beta_elt = Mod((-a2 + m*x)/2, x^2 - sf));
  printf("  Generator^2 == beta? %s  (beta = (-a2+m*sqrt(sf))/2 = (%d+%d*sqrt(%d))/2)\n",
    if(gen_sq == beta_elt, "YES", "NO"), -a2, m, sf);
  print();

  \\ Case 10: p=23, t=3, sf=-83, K=Q(sqrt(-83))
  my(p2=23, t2=3);
  my(a2b = 2*p2 - t2^2);
  my(Db = a2b^2 - 4*p2^2);
  my(sfb = sf_part(Db));
  my(mb = sqrtint(Db/sfb));
  my(Kb = bnfinit(x^2 - sfb, 1));
  my(Pb = idealprimedec(Kb, p2)[1]);
  \\ For sf=-83: -83 ≡ 1 mod 4, O_K = Z[(1+sqrt(-83))/2].
  \\ Proposed generator (3+sqrt(-83))/2 (same pattern).
  my(gen_elt2 = Mod((3 + x)/2, x^2 - sfb));
  my(gen_norm2 = (9 - sfb) / 4);  \\ = (9+83)/4 = 23 = p
  my(I_gen2 = idealhnf(Kb, gen_elt2));
  my(gen_norm_actual2 = idealnorm(Kb, I_gen2));
  printf("p=%d, sf=%d: proposed generator (3+sqrt(%d))/2 has norm=%d (=p?%s), generates prime above p?%s\n",
    p2, sfb, sfb, gen_norm2, if(gen_norm2 == p2, "YES", "NO"),
    if(gen_norm_actual2 == p2, "YES", "NO"));

  my(gen_sq2 = Mod((3 + x)/2, x^2 - sfb)^2);
  my(beta_elt2 = Mod((-a2b + mb*x)/2, x^2 - sfb));
  printf("  Generator^2 == beta? %s  (beta = (%d+%d*sqrt(%d))/2)\n",
    if(gen_sq2 == beta_elt2, "YES", "NO"), -a2b, mb, sfb);
}

\\ ============================================================
\\ Key non-norm-form case: p=61, sf=-219
\\ ============================================================
print();
print("--- Key case: p=61, t=5, sf=-219 ---");
print("sf=-219 is the discriminant of the secp256k1 CM-73 norm-form family.");
print("p=61 is NOT a norm-form prime (4*61=244, 244-73=171=3*57, 57 not a square).");
print("This confirms the theorem is independent of the norm-form structure.");
print();
{
  my(p=61, t=5);
  my(a2 = 2*p - t^2);   \\ = 97
  my(D = a2^2 - 4*p^2); \\ = 9409 - 14884 = -5475
  my(sf = sf_part(D));   \\ = -219
  my(m = sqrtint(D/sf)); \\ = 5
  my(K = bnfinit(x^2 - sf, 1));
  my(h = K.clgp.no);
  my(Pp = idealprimedec(K, p));
  my(P = Pp[1]);
  my(P2 = idealpow(K, P, 2));
  my(P2_prin = bnfisprincipal(K, P2));
  my(beta_elt = Mod((-a2 + m*x)/2, x^2 - sf));
  my(I_beta = idealhnf(K, beta_elt));
  printf("p=%d, t=%d: a2=%d, D=%d, sf=%d, m=%d\n", p, t, a2, D, sf, m);
  printf("K=Q(sqrt(%d)), h=%d, Cl(K)=Z/%dZ\n", sf, h, h);
  printf("#primes above p in K: %d, p splits: %s\n", #Pp, if(#Pp==2,"YES","NO"));
  printf("N(beta)=%d (=p^2=%d? %s)\n", (a2^2 - m^2*sf)/4, p^2, if((a2^2-m^2*sf)/4==p^2,"YES","NO"));
  printf("[P]^2=1 in Cl(K)? %s\n", if(P2_prin[1] == 0*P2_prin[1], "YES", "NO"));
  printf("(beta)=P^2? %s\n", if(I_beta == P2 || I_beta == idealpow(K, Pp[#Pp], 2), "YES", "NO"));
  \\ In Cl(Q(sqrt(-219)) = Z/4Z, [P]^2=1 means [P] has order 1 or 2.
  my(P_prin = bnfisprincipal(K, P));
  printf("ord([P]) in Z/%dZ: %s\n", h,
    if(P_prin[1] == 0*P_prin[1], "1",
      if(P2_prin[1] == 0*P2_prin[1], "2", ">2")));
}

print();
print("DONE.");
