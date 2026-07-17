\\ thread16_general_primes.gp
\\
\\ Thread 16: Verify the "universal order-2 Frobenius" theorem for non-norm-form primes.
\\
\\ THEOREM (Thread 15, GENERAL form -- proof is completely algebraic, no norm-form needed):
\\   Let p be ANY prime, a2 ANY integer with 0 < |a2| < 2p and p ∤ a2.
\\   Let D = a2^2 - 4p^2 < 0, sf = squarefree_part(D) (negative), m = sqrt(D/sf) ∈ Z>0.
\\   Then:
\\   (SPLIT) p ALWAYS splits in K = Q(sqrt(sf)): sf ≡ (a2*m^{-1})^2 (mod p), a nonzero square.
\\   (ORD2)  The prime ideal P of norm p in K satisfies [P]^2 = 1 in Cl(K).
\\
\\ KEY SPLITTING FACT: D = sf*m^2, so D ≡ a2^2 (mod p) [since p^2 ≡ 0 mod p... wait, D=a2^2-4p^2 ≡ a2^2 mod p].
\\   Thus sf*m^2 ≡ a2^2 (mod p).  Since p ∤ m (because p ∤ D), m is invertible mod p,
\\   so sf ≡ (a2/m)^2 (mod p) is a nonzero quadratic residue. Hence (sf/p) = 1 (Legendre).
\\   In particular, p cannot be inert in K when p ∤ a2. ✓
\\
\\ EXPERIMENT:
\\   10 non-norm-form primes p ∈ {101,103,107,113,127,131,137,139,149,151},
\\   each tested with 5 values of a2: {1, -1, p-1, p+1, -(2p-3)}.
\\   Additionally, a "sharpness check" with a2 = p (violates p ∤ a2 hypothesis)
\\   to test whether [P]^2 = 1 still holds (it should fail sometimes when p splits
\\   in Q(sqrt(-3)) and h(Q(sqrt(-3))) is divisible by 4, but for small p may still hold).
\\
\\ Run: gp -q thread16_general_primes.gp

default(parisize, 64000000);
default(timer, 0);

\\ ==============================================================
\\ Utility: squarefree part (with sign preserved)
\\ ==============================================================

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ==============================================================
\\ Verify the theorem for a single (p, a2) pair.
\\ Returns 1 if all checks pass, 0 if the pair is degenerate/excluded.
\\   mode=0: standard (p ∤ a2 required)
\\   mode=1: sharpness check (a2=p, p | a2 case; just report what [P]^ord is)
\\ ==============================================================

verify_pair(p, a2, mode) = {
  my(D, sf, m2, m, K, Pp, P, P2, P2_prin, ord_P, leg_sf, beta_elt,
     splitting_ok, pass_str, note);

  \\ --- Basic conditions ---
  if(abs(a2) == 0 || abs(a2) >= 2*p,
    printf("  [SKIP] a2=%d out of range for p=%d\n", a2, p);
    return(0));

  if(mode == 0 && (a2 % p == 0),
    printf("  [SKIP] p | a2 (a2=%d, p=%d)\n", a2, p);
    return(0));

  \\ --- Compute D, sf, m ---
  D  = a2^2 - 4*p^2;
  if(D >= 0, printf("  [ERROR] D=%d >= 0 for a2=%d, p=%d\n", D, a2, p); return(0));

  sf = sf_part(D);
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1,
    printf("  [ERROR] D/sf=%Ps not a positive integer\n", m2); return(0));
  m  = sqrtint(m2);
  if(m*m != m2, printf("  [ERROR] D/sf=%d not a perfect square\n", m2); return(0));

  \\ --- Verify (SPLIT): check Legendre symbol (sf/p) ---
  \\ Should be +1 when p ∤ a2 (squarefree sf ≡ (a2/m)^2 mod p)
  leg_sf = kronecker(sf, p);
  splitting_ok = (leg_sf == 1);

  \\ --- Compute ideal class group of K = Q(sqrt(sf)) ---
  K    = bnfinit(x^2 - sf, 1);
  Pp   = idealprimedec(K, p);

  \\ When p splits, #Pp = 2. If inert, #Pp = 1 (the inert prime has norm p^2).
  if(#Pp == 0,
    printf("  [WARN] no prime above p=%d in Q(sqrt(%d))?\n", p, sf);
    return(0));

  \\ For a split prime, pick the first factor of norm p.
  \\ For inert prime (norm p^2), [P] is principal => [P]^2 = 1 trivially.
  P = Pp[1];

  \\ Check order of [P] in Cl(K) via bnfisprincipal on P and P^2
  P2       = idealpow(K, P, 2);
  P2_prin  = bnfisprincipal(K, P2);

  \\ Determine order of [P]: 1, 2, or higher
  if(P2_prin[1] == 0*P2_prin[1],
    \\ P^2 is principal
    my(P1_prin = bnfisprincipal(K, P));
    if(P1_prin[1] == 0*P1_prin[1],
      ord_P = 1,   \\ P itself is principal
      ord_P = 2),  \\ P^2 is principal but P is not => ord = 2
    ord_P = -1    \\ P^2 is not principal => ord > 2
  );

  \\ Note for sharpness check (mode=1)
  note = if(mode == 1, " [SHARPNESS: a2=p]", "");

  pass_str = if(ord_P == 1 || ord_P == 2, "PASS", "FAIL");

  printf("  p=%-5d a2=%-7d sf=%-8d m=%-6d h=%-4d leg(sf/p)=%-2d #Pp=%d  ord[P]=%-3s  [P]^2=1:%s%s\n",
    p, a2, sf, m, K.clgp.no, leg_sf, #Pp,
    if(ord_P > 0, Str(ord_P), ">2"),
    if(ord_P == 1 || ord_P == 2, "YES", "NO"),
    note);

  if(ord_P == 1 || ord_P == 2, 1, 0)
};

\\ ==============================================================
\\ Main experiment
\\ ==============================================================

print("Thread 16: Universal order-2 Frobenius -- extension to non-norm-form primes");
print("============================================================================");
print();
print("THEOREM: For prime p, integer a2 with 0 < |a2| < 2p and p ∤ a2:");
print("  (SPLIT) p splits in K = Q(sqrt(sf)) where sf = sqf(a2^2 - 4p^2).");
print("  (ORD2)  [P]^2 = 1 in Cl(K) for the prime P of norm p in K.");
print();

{
  my(primes10, a2_sets, total, pass, r, sharpness_pass, sharpness_total);

  \\ 10 non-norm-form primes (none satisfy 73 + 3k^2 = 4p for odd k)
  primes10 = [101, 103, 107, 113, 127, 131, 137, 139, 149, 151];

  total = 0; pass = 0;
  sharpness_total = 0; sharpness_pass = 0;

  for(i = 1, #primes10,
    my(p = primes10[i]);
    printf("--- p = %d ---\n", p);

    \\ 5 standard a2 values (all with p ∤ a2 by construction)
    a2_sets = [1, -1, p-1, p+1, -(2*p-3)];
    for(j = 1, #a2_sets,
      r = verify_pair(p, a2_sets[j], 0);
      total++; pass += r);

    \\ Sharpness check: a2 = p (violates p ∤ a2 hypothesis)
    r = verify_pair(p, p, 1);
    sharpness_total++;
    sharpness_pass += r;

    print());

  printf("Standard cases: %d/%d passed ([P]^2=1 with p ∤ a2 hypothesis)\n",
    pass, total);
  printf("Sharpness (a2=p): %d/%d had [P]^2=1 (hypothesis violated but may still hold trivially)\n",
    sharpness_pass, sharpness_total);

  if(pass == total,
    print("\nALL STANDARD CASES CONFIRMED. The general theorem holds for 10 non-norm-form primes × 5 a2 values."),
    printf("\nWARNING: %d cases FAILED.\n", total - pass));
}

print();
print("REMARK (splitting): When p ∤ a2, leg(sf/p) should always equal 1.");
print("  Proof: D = a2^2 - 4p^2 ≡ a2^2 (mod p), and D = sf*m^2.");
print("  Since p ∤ D (p ∤ a2 => p ∤ a2^2 => p ∤ D), p is invertible mod D/sf = m^2,");
print("  so sf ≡ D*(m^{-1})^2 ≡ (a2*m^{-1})^2 (mod p), a nonzero square.");

{
  \\ Quick sweep: verify leg(sf/p) = 1 for ALL tested pairs with p ∤ a2
  my(primes10, a2_sets, leg_all_1);
  primes10 = [101, 103, 107, 113, 127, 131, 137, 139, 149, 151];
  leg_all_1 = 1;
  for(i = 1, #primes10,
    my(p = primes10[i]);
    a2_sets = [1, -1, p-1, p+1, -(2*p-3)];
    for(j = 1, #a2_sets,
      my(a2 = a2_sets[j], D, sf, m2);
      D  = a2^2 - 4*p^2;
      sf = sf_part(D);
      if(kronecker(sf, p) != 1, leg_all_1 = 0;
        printf("  SPLIT FAIL: p=%d a2=%d leg=%d\n", p, a2, kronecker(sf,p)))));
  if(leg_all_1, print("  Verified: leg(sf/p) = 1 for ALL 50 standard pairs. ✓"),
                print("  ERROR: some split check failed."));
}

print();
print("DONE.");
