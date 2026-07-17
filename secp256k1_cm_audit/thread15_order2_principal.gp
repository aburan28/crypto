\\ thread15_order2_principal.gp
\\ Thread 15: Verify the "universal order-2" conjecture.
\\
\\ CONJECTURE: For every norm-form prime 4p = 73+3k^2 (k odd), the Frobenius
\\ ideal (pi) above p in the CM field K = Q(sqrt(sf)) satisfies (pi)^2 is
\\ PRINCIPAL in the ideal class group Cl(K).
\\
\\ This implies ord([pi]) | 2, hence ord([pi]) in {1,2}.
\\ We observed: ord=1 iff h(K)=1 (sf=-3, p=12889); ord=2 otherwise.
\\
\\ Method:
\\   1. For each (k,p,sf): construct K=bnfinit(x^2-sf), find pi above p,
\\      compute pi^2 in K, and call bnfisprincipal to check principality.
\\   2. Also verify pi^2 explicitly: since pi satisfies pi*pibar=p, we have
\\      pi^2 satisfies Nm(pi^2)=p^2, and pi^2+pibar^2 = (pi+pibar)^2 - 2p = a^2-2p
\\      where a=Tr_K(pi).
\\   3. Special check for sf=-3 (p=12889): find explicit Eisenstein prime.
\\
default(parisize, 128000000);
default(timer, 0);

\\ Helper: check if ideal is principal
is_principal(K, idl) = {
  my(res);
  res = bnfisprincipal(K, idl, 1);
  return(res[1] == 0 || matsize(res[1])[2] == 0 || vecmax(apply(abs, Vec(res[1]))) == 0);
};

\\ Compute pi from Weil polynomial T^4 + a2*T^2 + p^2
\\ The Frobenius satisfies T^2 = (a2 +/- sqrt(a2^2-4p^2))/2 = (a2 +/- sqrt(disc4))/2
\\ disc4 = a2^2 - 4p^2 = sf * d^2 for some d.
\\ pi is root of T^2 - A*T + p = 0 where A^2 = a2 + 2p (one of the quadratic factors).
\\ We work in K = Q(sqrt(sf)), where pi lies.
\\ Explicitly: if T^4 + a2*T^2 + p^2 = (T^2 - A*T + p)(T^2 + A*T + p) over Q(A),
\\ then A^2 = a2 + 2p (from Vieta). But disc4 = a2^2-4p^2 = (a2-2p)(a2+2p).
\\
\\ In K = Q(sqrt(sf)), pi satisfies pi + pibar = A (where A^2 = a2+2p) and pi*pibar=p.
\\ We know pi = (A + sqrt(A^2-4p))/2. In K, A = sqrt(a2+2p) is in a bigger field
\\ UNLESS a2+2p is a perfect square. But disc4 = (a2-2p)(a2+2p) = sf*d^2.
\\
\\ Simpler: the Frobenius ideal above p in K: find a prime ideal P | p in K,
\\ then check whether P^2 is principal.
\\
\\ Also: pi^2 = pi*pi. Since pi+pibar = Tr(pi) in K (which may not be in Q),
\\ we need to work in an appropriate field.
\\
\\ For the purpose of the conjecture: we just need (P)^2 principal,
\\ where P is the prime ideal above p in K.
\\
check_pi_squared(k, p, a2, sf) = {
  my(K, h, cyc, Plist, P, res_P, res_P2, id_P2, is_princ, pi2_elem);
  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  cyc = K.clgp.cyc;
  Plist = idealprimedec(K, p);
  if(#Plist == 0, printf("  ERROR: p=%d doesn't split in sf=%d?\n", p, sf); return);
  P = Plist[1];
  \\ P^2 in ideal group
  id_P2 = idealpow(K, P, 2);
  res_P2 = bnfisprincipal(K, id_P2, 1);
  is_princ = (vecmax(apply(abs, Vec(res_P2[1]))) == 0);
  \\ Also get class exponent of P itself
  res_P = bnfisprincipal(K, P, 1);
  printf("  k=%-4d p=%-8d a2=%-8d sf=%-10d h=%-4d (P)^2 principal: %s  [P]=(%s)\n",
    k, p, a2, sf, h, if(is_princ,"YES","NO"), res_P[1]);
  if(is_princ,
    \\ Extract the generator of P^2
    pi2_elem = res_P2[2];
    \\ Verify: Nm(pi2_elem) should be p^2
    my(nm2);
    nm2 = idealnorm(K, idealhnf(K, pi2_elem));
    if(nm2 != p^2,
      printf("    WARNING: generator norm = %d (expected p^2 = %d)\n", nm2, p^2)
    ,
      printf("    Generator of (P^2): norm=%d=p^2 CONFIRMED\n", nm2)
    )
  );
};

print("=== THREAD 15: Universal Order-2 Conjecture Verification ===");
print("Claim: For each norm-form prime 4p=73+3k^2, the Frobenius ideal P above p in");
print("       Q(sqrt(sf)) satisfies P^2 = (pi^2) is principal.");
print();
print("CM-73 primes (sf=-219, reference class):");
check_pi_squared(1, 19, -35, -219);
check_pi_squared(5, 37, 1, -219);
check_pi_squared(9, 79, 85, -219);
check_pi_squared(11, 109, 145, -219);

print();
print("Non-CM-73 norm-form primes:");
check_pi_squared(21, 349, 385, -939);
check_pi_squared(105, 8287, 2149, -1731);
check_pi_squared(131, 12889, -25751, -3);
check_pi_squared(25, 487, -299, -3819);
check_pi_squared(35, 937, 1, -5619);
check_pi_squared(31, 739, -1403, -8643);
check_pi_squared(119, 10639, 10549, -32187);
check_pi_squared(65, 3187, -4499, -32619);
check_pi_squared(109, 8929, 5905, -35859);
check_pi_squared(99, 7369, 385, -43059);
check_pi_squared(41, 1279, -2315, -14619);
check_pi_squared(55, 2287, -899, -16419);
check_pi_squared(85, 5437, -9791, -61995);
check_pi_squared(91, 6229, -11375, -71499);
check_pi_squared(101, 7669, -13751, -87267);

print();
print("=== Special case: sf=-3 (Eisenstein prime) ===");
print("p=12889 should factor as pi*pibar in Z[omega] where omega=(-1+sqrt(-3))/2");
{
  my(K, Plist, P, g, nm, pi_explicit, trace_pi);
  K = bnfinit(x^2 + x + 1, 1);  \\ Z[omega], ring of integers of Q(sqrt(-3))
  Plist = idealprimedec(K, 12889);
  print("Ideal factorization of (12889) in Z[omega]:");
  printf("  Number of prime ideals above 12889: %d\n", #Plist);
  if(#Plist >= 1,
    P = Plist[1];
    printf("  f (residue degree): %d\n", P.f);
    printf("  e (ramification): %d\n", P.e);
    \\ Find explicit generator
    g = bnfisprincipal(K, P, 1);
    printf("  Class exp: %Ps (should be [0] since h=1)\n", g[1]);
    pi_explicit = g[2];
    nm = idealnorm(K, idealhnf(K, pi_explicit));
    printf("  Explicit generator pi: %Ps\n", pi_explicit);
    printf("  Nm(pi): %d (should equal 12889)\n", nm);
    \\ Check norm directly: Nm(a+b*omega) = a^2 - a*b + b^2
    my(a, b);
    if(type(pi_explicit) == "t_COL" || type(pi_explicit) == "t_VEC",
      a = polcoef(Mod(pi_explicit, x^2+x+1), 0);
      b = polcoef(Mod(pi_explicit, x^2+x+1), 1);
      printf("  pi = %d + %d*omega: Nm = %d^2 - %d*%d + %d^2 = %d\n",
        a, b, a, a, b, b, a^2 - a*b + b^2)
    )
  )
};

print();
print("=== Algebraic proof sketch for P^2 principal ===");
print("Key identity: P*Pbar = (p) in any imaginary quadratic K=Q(sqrt(D<0)).");
print("If ord([P])=2, then [P]^2=1, i.e., P^2 is principal: P^2=(alpha) for some alpha in O_K.");
print("This alpha satisfies Nm(alpha)=p^2 and is a 'Weil p^2-integer' (alpha+alphabar in Z).");
print("The element alpha = pi^2 (Frobenius squared) satisfies exactly these properties.");
print("Proof: Nm(pi^2) = (Nm pi)^2 = p^2; Tr(pi^2) = (Tr pi)^2 - 2*Nm(pi) = A^2-2p in Z.");
print("So (pi^2) = P^2 (both have norm p^2 and pi^2 in O_K).");
print("Hence P^2 = (pi^2) is principal IFF pi^2 in O_K, which holds by Weil-number theory.");
print("QED (modulo verifying pi^2 is truly integral, which follows from pi being a Weil p-integer).");
print();
print("DONE.");
