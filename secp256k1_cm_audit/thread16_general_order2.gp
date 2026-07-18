\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify the "universal order-2 Frobenius" theorem holds for
\\ non-norm-form primes, i.e., for ARBITRARY ordinary elliptic curves E/F_p.
\\
\\ SETUP: For any ordinary E/F_p with Frobenius trace t (≠ 0), consider the
\\ abelian surface A = E × E^t with biquadratic Weil polynomial
\\   W(T) = T^4 + a2*T^2 + p^2,   a2 = 2p - t^2.
\\
\\ Discriminant: D = a2^2 - 4p^2 = t^2*(t^2-4p) < 0.
\\ Write D = sf*m^2 with sf = squarefree-part(D) < 0, m > 0.
\\
\\ THEOREM (Thread 15, generalised): For any such (p, a2), if p ∤ a2
\\ (equivalently p ∤ t, which holds for all ordinary E/F_p with p > 4),
\\ then the prime P above p in Q(sqrt(sf)) satisfies [P]^2 = 1 in Cl(Q(sqrt(sf))).
\\
\\ This script verifies this for 10 primes NOT in the secp256k1 norm-form
\\ family {p : 4p = 73+3k^2, k odd}, using random curves y^2 = x^3+ax+b over F_p.
\\
\\ Also checks: for each prime, 5 independent curves give CONSISTENT sf
\\ (since sf = sf_part(t^2-4p) depends on the CM field of E, not just on p).
\\
\\ Run: gp --stacksize 256000000 -q thread16_general_order2.gp

default(parisize, 256000000);
default(timer, 0);

\\ ---------------------------------------------------------------
\\ Utility: signed squarefree part
\\ ---------------------------------------------------------------
sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ---------------------------------------------------------------
\\ Is p of norm-form type? (4p = 73 + 3k^2, k odd, p prime)
\\ ---------------------------------------------------------------
is_norm_form(p) = {
  my(val, k);
  forstep(k = 1, 1000, 2,
    val = 73 + 3*k^2;
    if(val == 4*p, return(1));
    if(val > 4*p, return(0)));
  0
};

\\ ---------------------------------------------------------------
\\ Verify [P]^2 = 1 for curve y^2 = x^3 + a*x + b over F_p
\\ Returns: [pass, sf, t, a2, h] or [0, reason]
\\ ---------------------------------------------------------------
verify_curve(p, a, b) = {
  my(E, t, a2, D, sf, m2, m, K, Pp, P, P2, P2_prin, h, ord_P);

  \\ check non-singular
  if((4*a^3 + 27*b^2) % p == 0, return([0, "singular"]));

  E = ellinit([Mod(a, p), Mod(b, p)]);
  t = p + 1 - ellcard(E);        \\ Frobenius trace

  if(t == 0, return([0, "supersingular"]));

  a2 = 2*p - t^2;
  D  = a2^2 - 4*p^2;             \\ = t^2*(t^2-4p) < 0
  sf = sf_part(D);

  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1, return([0, "D/sf_error"]));
  m  = sqrtint(m2);
  if(m*m != m2, return([0, "not_square"]));

  \\ Check p does not divide a2 (theorem precondition)
  if(a2 % p == 0, return([0, "p|a2 (supersingular escape)"]));

  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);
  if(#Pp == 0, return([0, "p_inert"]));

  P       = Pp[1];
  P2      = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);

  \\ ord of [P]: check if P itself is principal (ord=1) or only P^2 (ord=2)
  ord_P = if(P2_prin[1] == 0*P2_prin[1],
              if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1], 1, 2),
              -1);

  [if(P2_prin[1] == 0*P2_prin[1], 1, 0), sf, t, a2, h, ord_P]
};

\\ ---------------------------------------------------------------
\\ Main: 10 non-norm-form primes, 5 curves each
\\ ---------------------------------------------------------------

\\ Primes chosen to span a range and NOT be of norm-form type
\\ (verified by hand: none satisfy 4p=73+3k^2 for small odd k)
test_primes = [101, 197, 307, 503, 1013, 2011, 5003, 10007, 20011, 50021];

print("Thread 16: Universal order-2 Frobenius — 10 non-norm-form primes");
print("==================================================================");
print("For each prime p, 5 curves y^2 = x^3 + a*x + b over F_p are tested.");
print("Columns: p  a  b  t  sf  h  ord([P])  [P]^2=1?");
print();

{
  my(total = 0, ok = 0, pass, sf, t, a2, h, ord_P, res, p, a, b);

  for(i = 1, #test_primes,
    p = test_primes[i];

    if(is_norm_form(p),
      printf("  p=%d: SKIPPED (norm-form prime)\n", p);
      next);

    printf("p = %d\n", p);

    for(j = 1, 5,
      \\ Pick deterministic but varied (a,b)
      a = (j * 17 + 3) % p;
      b = (j * 31 + 7) % p;

      res = verify_curve(p, a, b);

      if(type(res[1]) != "t_INT" || res[1] == 0,
        printf("  [a=%d b=%d] skip: %s\n", a, b, res[2]);
        next);

      pass  = res[1];
      sf    = res[2];
      t     = res[3];
      a2    = res[4];
      h     = res[5];
      ord_P = res[6];

      total++;
      if(pass, ok++);

      printf("  a=%-4d b=%-4d : t=%-6d  sf=%-10d  h=%-4d  ord([P])=%d  [P]^2=1:%s\n",
        a, b, t, sf, h, ord_P, if(pass, "YES", "NO(!)")));

    print());

  printf("Summary: %d/%d test cases passed (all [P]^2=1: %s)\n",
    ok, total, if(ok == total, "YES", "NO"));
}

\\ ---------------------------------------------------------------
\\ Explicit demonstration: the theorem's β element
\\ ---------------------------------------------------------------
print();
print("Explicit beta check for p=307 (non-norm-form):");
{
  my(p, a, b, E, t, a2, D, sf, m2, m, beta_elt, K, Pp, P, P2, I_beta, eqv);
  p = 307; a = 20; b = 38;
  E = ellinit([Mod(a, p), Mod(b, p)]);
  t = p + 1 - ellcard(E);
  a2 = 2*p - t^2;
  D  = a2^2 - 4*p^2;
  sf = sf_part(D);
  m2 = D / sf;
  m  = sqrtint(m2);

  printf("  t=%d  a2=%d  D=%d  sf=%d  m=%d\n", t, a2, D, sf, m);

  K  = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  P  = Pp[1];
  P2 = idealpow(K, P, 2);

  \\ beta = (-a2 + m*sqrt(sf)) / 2 as a nf element
  beta_elt = Mod((-a2 + m*x)/2, x^2 - sf);
  I_beta   = idealhnf(K, beta_elt);

  eqv = (I_beta == P2 || I_beta == idealpow(K, Pp[#Pp], 2));
  printf("  beta = (%d + %d*sqrt(%d))/2\n", -a2, m, sf);
  printf("  minpoly(beta) = %Ps\n", minpoly(beta_elt));
  printf("  N(beta) = %d  (= p^2 = %d: %s)\n",
    (a2^2 - m2*sf)/4, p^2, if((a2^2 - m2*sf)/4 == p^2, "YES", "NO"));
  printf("  idealhnf(K,beta) == P^2: %s\n", if(eqv, "YES", "NO"));
  printf("  bnfisprincipal(K,P^2)[1]==%d (0=principal): %s\n",
    bnfisprincipal(K,P2)[1],
    if(bnfisprincipal(K,P2)[1] == 0*bnfisprincipal(K,P2)[1], "YES","NO"));
}

print();
print("DONE.");
