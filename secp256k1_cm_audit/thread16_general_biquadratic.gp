\\ ================================================================
\\ Thread 16: Universal Order-2 Frobenius — General Biquadratic Case
\\ ================================================================
\\
\\ Context: Thread 15 proved that for norm-form primes p (4p=73+3k^2),
\\   the Frobenius ideal [P] in Cl(Q(sqrt(sf))) satisfies [P]^2 = 1.
\\   The proof (steps A–E) used only:
\\     (A) beta = (-a2 + m*sqrt(sf))/2 is an alg. int. (minpoly x^2+a2*x+p^2)
\\     (B) beta in K = Q(sqrt(sf))
\\     (C) N_{K/Q}(beta) = p^2
\\     (D) p does not divide a2
\\     (E) therefore (beta) = P^2 or Pbar^2 => [P]^2 = 1
\\
\\ The proof makes NO use of the norm-form condition 4p=73+3k^2.
\\ This script verifies the theorem for 10 non-norm-form primes.
\\
\\ Construction: for prime p and integer alpha with 0 < |alpha| < 2*sqrt(p)
\\   and p ∤ alpha, set:
\\     a2 = 2*p - alpha^2   (so T^4+a2*T^2+p^2 = (T^2+alpha*T+p)(T^2-alpha*T+p))
\\     D  = a2^2 - 4*p^2 = alpha^2*(alpha^2 - 4*p) < 0  (since |alpha|<2*sqrt(p))
\\     sf = squarefree kernel of D (negative)
\\     m  = sqrt(D/sf)
\\     beta = (-a2 + m*sqrt(sf))/2
\\
\\ All five conditions (A)–(E) hold automatically; we verify numerically.
\\ ================================================================

\\ Helper: check and verify one (p, alpha) instance
verify_case(p, alpha) = {
  my(a2, D, sf, m2, m, beta_elt, norm_beta, K, h, Pp, P, P2, P2_prin,
     P_prin, p2_is_princ, p_is_princ, ord_str);

  a2 = 2*p - alpha^2;
  D  = a2^2 - 4*p^2;   \\ = alpha^2*(alpha^2 - 4*p) < 0

  if(D >= 0,
    printf("  p=%-9d alpha=%-4d: D >= 0 (%d), skipping\n", p, alpha, D);
    return(-1));

  \\ Squarefree part of D (negative since D<0)
  sf = core(D);
  m2 = D \ sf;
  if(m2 <= 0,
    printf("  p=%-9d alpha=%-4d: D/sf <= 0 (%d), skipping\n", p, alpha, m2);
    return(-1));
  if(!issquare(m2),
    printf("  p=%-9d alpha=%-4d: D/sf=%d not a perfect square\n", p, alpha, m2);
    return(-1));
  m = sqrtint(m2);

  \\ (D): p must not divide a2
  if(a2 % p == 0,
    printf("  p=%-9d alpha=%-4d: p|a2 (a2=%d) — degenerate case, skip\n", p, alpha, a2);
    return(-1));

  \\ (A)+(B): beta = (-a2 + m*sqrt(sf))/2 in K=Q(sqrt(sf))
  beta_elt = Mod((-a2 + m*x)/2, x^2 - sf);

  \\ (C): Norm of beta = p^2
  norm_beta = (a2^2 - m2*sf) \ 4;

  \\ Verify norm numerically
  if(norm_beta != p^2,
    printf("  p=%-9d alpha=%-4d: NORM FAILURE (norm=%d != p^2=%d)\n", p, alpha, norm_beta, p^2);
    return(0));

  \\ (E): verify in Cl(K)
  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  Pp = idealprimedec(K, p);

  if(#Pp == 0,
    printf("  p=%-9d alpha=%-4d sf=%-10d: p inert in K, no prime above p\n", p, alpha, sf);
    return(-1));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2, 1);
  P_prin  = bnfisprincipal(K, P,  1);

  p2_is_princ = (matsize(P2_prin[1]) == [0,0] || P2_prin[1] == 0*P2_prin[1]);
  p_is_princ  = (matsize(P_prin[1])  == [0,0] || P_prin[1]  == 0*P_prin[1]);

  ord_str = if(p_is_princ, "1", if(p2_is_princ, "2", ">2(ERR)"));

  printf("  p=%-9d a=%-4d a2=%-10d sf=%-11d m=%-7d h=%-5d N=p^2:%s p∤a2:YES [P]^2=1:%-8s ord(%s)\n",
    p, alpha, a2, sf, m, h,
    if(norm_beta == p^2, "YES","NO "),
    if(p2_is_princ, "YES","NO(ERR)"),
    ord_str);

  if(p2_is_princ, 1, 0)
};

\\ ================================================================
\\ Check that p is NOT of the norm-form 4p=73+3k^2
\\ ================================================================
is_norm_form(p) = {
  my(r);
  r = 4*p - 73;
  if(r <= 0 || r % 3 != 0, return(0));
  r = r \ 3;
  issquare(r)
};

\\ ================================================================
\\ Main: select 10 diverse non-norm-form primes and test alpha values
\\ ================================================================

print("Thread 16: General Biquadratic Weil Polynomial — Order-2 Frobenius Theorem");
print("==========================================================================");
print("Testing T^4 + a2*T^2 + p^2 with a2 = 2p - alpha^2 (|alpha| < 2*sqrt(p)).");
print("Theorem: if p ∤ a2 then [P]^2 = 1 in Cl(Q(sqrt(sf(D)))), D=a2^2-4p^2.");
print("Goal: verify for 10 non-norm-form primes with multiple alpha values each.");
print();

{
  my(test_cases, total, pass, r);
  \\ Format: [p, [alpha values]]
  \\ Choose primes NOT of the form 4p=73+3k^2 (verified by is_norm_form above)
  \\ and diverse alpha values to exercise different CM fields.
  \\
  \\ Small primes:  p=23,29,31,43,53,67,71,83,97,101
  \\ Check norm-form status inline.
  test_cases = [
    [23,   [1, 2, 3, 4]],
    [29,   [1, 3, 5]],
    [31,   [1, 4, 5]],
    [43,   [1, 3, 5, 7]],
    [53,   [1, 4, 7]],
    [67,   [1, 5, 8]],
    [71,   [1, 6, 8]],
    [97,   [1, 7, 9, 12]],
    [101,  [1, 8, 13]],
    [127,  [1, 9, 15, 21]],
    [983,  [1, 15, 31]],
    [9973, [1, 30, 77, 150]]
  ];

  total = 0; pass = 0;

  for(i = 1, #test_cases,
    my(p = test_cases[i][1], alphas = test_cases[i][2]);
    if(!isprime(p), next);
    if(is_norm_form(p),
      printf("p=%d IS norm-form — skipping (shouldn't happen)\n", p);
      next);
    printf("--- p = %d (NOT norm-form) ---\n", p);
    for(j = 1, #alphas,
      my(alpha = alphas[j]);
      if(alpha >= 2*sqrtint(p), next);   \\ need |alpha| < 2*sqrt(p)
      r = verify_case(p, alpha);
      if(r == 1, total++; pass++);
      if(r == 0, total++));
    print());

  printf("Summary: %d/%d non-norm-form (p,alpha) cases passed [P]^2=1.\n", pass, total);
  printf("All-pass: %s\n", if(pass == total, "YES — theorem confirmed for non-norm-form primes", "NO — FAILURES PRESENT"));
}

print();
print("THEOREM STATEMENT (Thread 16 generalisation):");
print("  Let p be any prime, alpha any integer with 0 < |alpha| < 2*sqrt(p) and p ∤ alpha.");
print("  Set a2 = 2p - alpha^2, D = a2^2 - 4p^2, sf = squarefree(D), m = sqrt(D/sf).");
print("  In K = Q(sqrt(sf)), write (p) = P*Pbar. Then [P]^2 = 1 in Cl(K).");
print("  The same proof (A)-(E) applies regardless of whether p is a norm-form prime.");
print();
print("DONE.");
