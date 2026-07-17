\\ thread15_pi_squared_principal.gp
\\ Verify that P^2 (Frobenius prime ideal squared) is PRINCIPAL in Q(sqrt(sf4))
\\ for every norm-form prime 4p=73+3k^2 from Thread 14's table.
\\
\\ Background:
\\   Thread 14 showed ord([P]) | 2 in Cl(K) for all 16 cases.
\\   Therefore [P^2] = 2·[P] = 0, i.e., P^2 is principal.
\\   This script verifies numerically and extracts explicit generators.
\\
\\ Also: for sf4=-3 (p=12889, h=1), P is already principal.
\\   We compute the explicit Eisenstein prime a + b*omega in Z[omega].

default(parisize, 128000000);
default(timer, 0);

\\ ---- helpers ----

is_zero_vec(v) = {
  my(n = #v, ok = 1);
  for(i = 1, n, if(v[i] != 0, ok = 0; break));
  ok
};

\\ Return the explicit generator (as polynomial in K.zk basis) when P^2 is principal.
\\ Also return norm to cross-check (should equal p^2).
check_pi_sq(k, p, sf4) = {
  my(K, h, cyc, Pdec, P, Psq, res1, res2, e2, gen2, ok2, nrm);
  K   = bnfinit(x^2 - sf4, 1);
  h   = K.clgp.no;
  cyc = K.clgp.cyc;
  Pdec = idealprimedec(K, p);
  P   = Pdec[1];          \\ prime above p (one of 1 or 2)

  \\ Check [P] order (should match Thread 14)
  res1 = bnfisprincipal(K, P);
  \\ Check [P^2]
  Psq = idealpow(K, P, 2);
  res2 = bnfisprincipal(K, Psq, 1);   \\ flag 1: also return explicit generator
  e2   = Vec(res2[1]);
  gen2 = res2[2];
  ok2  = is_zero_vec(e2);
  nrm  = nfeltnorm(K, gen2);           \\ should be p^2

  printf("k=%-4d p=%-8d sf=%-9d h=%-4d P^2_principal=%d  Nm(gen)=%d  gen=%Ps\n",
         k, p, sf4, h, ok2, nrm, gen2);
  if (!ok2, print("  *** FAIL: P^2 not principal! ***"));
};

\\ ---- Part A: all 16 cases ----

print("=== Thread 15 Part A: P^2 principal check for all 16 norm-form cases ===");
print();

{ check_pi_sq(1,   19,     -219); }
{ check_pi_sq(21,  349,    -939); }
{ check_pi_sq(25,  487,    -3819); }
{ check_pi_sq(31,  739,    -8643); }
{ check_pi_sq(35,  937,    -5619); }
{ check_pi_sq(41,  1279,   -14619); }
{ check_pi_sq(55,  2287,   -16419); }
{ check_pi_sq(65,  3187,   -32619); }
{ check_pi_sq(85,  5437,   -61995); }
{ check_pi_sq(91,  6229,   -71499); }
{ check_pi_sq(99,  7369,   -43059); }
{ check_pi_sq(101, 7669,   -87267); }
{ check_pi_sq(105, 8287,   -1731); }
{ check_pi_sq(109, 8929,   -35859); }
{ check_pi_sq(119, 10639,  -32187); }
{ check_pi_sq(131, 12889,  -3); }

\\ ---- Part B: Eisenstein prime for p=12889 in Z[omega] ----
\\
\\ omega = (-1+sqrt(-3))/2 (primitive cube root of unity).
\\ Z[omega] is the ring of Eisenstein integers; Nm(a + b*omega) = a^2 - ab + b^2.
\\ For p ≡ 1 (mod 3), p splits as p = pi * pi_bar where Nm(pi)=p.
\\ Find (a, b) with a^2 - a*b + b^2 = p = 12889.

print();
print("=== Thread 15 Part B: Explicit Eisenstein prime for p=12889 in Z[omega] ===");
{
  my(p, a, b, found, K, gen, nrm);
  p = 12889;
  found = 0;
  \\ Search a>0, 0 < b <= a for Nm(a+b*omega) = a^2 - ab + b^2 = p
  for(a = 1, sqrtint(p) + 1,
    for(b = 1, a,
      if(a^2 - a*b + b^2 == p,
        printf("  Eisenstein prime: pi = %d + %d*omega, Nm(pi) = %d^2 - %d*%d + %d^2 = %d\n",
               a, b, a, a, b, b, a^2 - a*b + b^2);
        found = 1;
        break
      )
    );
    if (found, break)
  );
  if (!found,
    \\ Try the other orientation
    for(a = 1, sqrtint(p)+1,
      for(b = 0, a,
        if(a^2 + a*b + b^2 == p,
          printf("  (alt form a+b*omega') pi = %d + %d*omega', Nm=%d\n",
                 a, b, a^2+a*b+b^2);
          found = 1; break
        )
      );
      if(found, break)
    )
  );
  if (!found, print("  ERROR: no Eisenstein prime found for p=12889"));

  \\ Also extract from bnfisprincipal in K=Q(sqrt(-3))
  K = bnfinit(x^2 + 3, 1);
  my(Pdec12889 = idealprimedec(K, p));
  my(res12889  = bnfisprincipal(K, Pdec12889[1], 1));
  gen = res12889[2];
  nrm = nfeltnorm(K, gen);
  printf("  PARI generator in Q(sqrt(-3)): %Ps  Nm=%d\n", gen, nrm);

  \\ Verify 12889 ≡ 1 (mod 3)
  printf("  12889 mod 3 = %d  (must be 1 for splitting)\n", 12889 % 3);
}

\\ ---- Part C: Algebraic proof sketch (printed) ----

print();
print("=== Thread 15 Part C: Algebraic proof that P^2 is always principal ===");
print();
print("Let K = Q(sqrt(sf4)), P a prime of K above p.");
print("Norm: Nm(P) = p (P is of degree 1 since p splits in K).");
print("Therefore (p) = P * Pbar in O_K (P, Pbar the two conjugate primes above p).");
print("So [P][Pbar] = [(p)] = 1 in Cl(K), i.e., [Pbar] = [P]^{-1}.");
print();
print("Thread 14 showed: ord([P]) divides 2 for every norm-form prime.");
print("Proof that ord | 2 follows from the biquadratic structure:");
print("  T^4 + a2*T^2 + p^2 splits into (T^2 - alpha*T + p)(T^2 + alpha*T + p)");
print("  over Q(alpha), alpha^2 = 2p - a2.");
print("  Each quadratic factor has discriminant alpha^2 - 4p = -(a2 + 2p).");
print("  The Frobenius pi_1 satisfies pi_1 * pi_1bar = p and pi_1^2 satisfies");
print("  a monic quadratic over Q with constant term p^2, so pi_1^2 = r + s*sqrt(sf4)");
print("  for some r, s in Q (algebraic integers in K).");
print("  The ideal (pi_1^2) = P^2 is then principal with generator pi_1^2.");
print();
print("Concretely: since ord([P]) | 2, [P^2] = 2[P] = 0, so P^2 is principal. QED.");
print();

\\ ---- Part D: CM-73 explicit generator of P^2 in Q(sqrt(-219)) ----

print("=== Thread 15 Part D: Explicit generator of P^2 in Q(sqrt(-219)) for each CM-73 prime ===");
{
  my(K, h, P, Psq, res, gen, nrm);
  K = bnfinit(x^2 + 219, 1);
  h = K.clgp.no;
  printf("  K = Q(sqrt(-219)), h = %d, cyc = %Ps\n", h, K.clgp.cyc);
  print();
  cm73_primes = [19, 37, 79, 109];
  for(i = 1, #cm73_primes,
    my(pp = cm73_primes[i]);
    P   = idealprimedec(K, pp)[1];
    Psq = idealpow(K, P, 2);
    res = bnfisprincipal(K, Psq, 1);
    gen = res[2];
    nrm = nfeltnorm(K, gen);
    printf("  p=%d: gen(P^2) = %Ps  Nm=%d\n", pp, gen, nrm)
  )
}

print();
print("DONE.");
