\\ thread16_general_weil.gp
\\
\\ Thread 16: Order-2 theorem generalises beyond the secp256k1 norm-form family.
\\
\\ CONSTRUCTION: For prime p > 4 and t with 0 < t^2 < 4p, the abelian surface
\\ E x E^t (E has Frobenius trace t; E^t, its quadratic twist, has trace -t)
\\ has biquadratic Weil polynomial
\\
\\   (T^2 - t*T + p)(T^2 + t*T + p) = T^4 + a2*T^2 + p^2,  a2 = 2p - t^2.
\\
\\ KEY LEMMA (Thread 16, new): p always SPLITS in K = Q(sqrt(sf)) where
\\   sf = sf(t^2*(t^2-4p)) = -core(4p-t^2).
\\
\\   Proof: Let n = 4p-t^2 > 0.
\\   (1) p does not divide n: n = -t^2 (mod p), and p does not divide t^2
\\       (since |t|<2*sqrt(p)<p for p>4 and t!=0). So p does not divide n.
\\   (2) (core(n)/p) = (n/p) [Legendre]: write n=core(n)*k^2; p does not
\\       divide k (else p^2|n, but p does not divide n), so (k^2/p)=1.
\\   (3) (n/p) = (-t^2/p) = (-1/p)*(t^2/p) = (-1/p)*1 = (-1/p).
\\   (4) (sf/p) = (-core(n)/p) = (-1/p)*(core(n)/p) = (-1/p)*(-1/p) = 1.
\\   => p splits in K. QED.
\\
\\ Combined with:
\\   (c) p does not divide a2 (since a2 = -t^2 mod p != 0),
\\ the Thread 15 proof (A)-(E) gives [P]^2 = 1 in Cl(K) for all such (p,t).
\\
\\ Run: gp --stacksize 128000000 -q thread16_general_weil.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_norm_form(p) = {
  my(val, k2);
  val = 4*p - 73;
  if(val <= 0 || val % 3 != 0, return(0));
  k2 = val \ 3;
  issquare(k2)
};

\\ Verify Thread-15 proof conditions for one (p, t) pair (product-of-twists).
verify_general(p, t) = {
  my(a2, n, D, sf, m2, m, K, Pp, P, P2, P2_prin, norm_beta, kron_sfp, h, ok);
  if(!isprime(p), return(-2));
  if(is_norm_form(p), return(-3));
  if(t == 0 || t^2 >= 4*p, return(-4));

  a2 = 2*p - t^2;
  n  = 4*p - t^2;
  D  = a2^2 - 4*p^2;        \\ = -t^2 * n < 0
  sf = sf_part(D);           \\ = -core(n)
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1, return(-5));
  m = sqrtint(m2);
  if(m*m != m2, return(-6));

  norm_beta = (a2^2 - m2*sf) / 4;

  kron_sfp = kronecker(sf, p);

  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  Pp = idealprimedec(K, p);
  if(#Pp == 0, return(-7));   \\ p inert: cannot happen by Key Lemma

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);
  ok = (P2_prin[1] == 0 * P2_prin[1]);

  printf("p=%-6d t=%-4d a2=%-8d sf=%-11d m=%-7d h=%-4d  N=p^2:%s  p|a2:%s  (sf/p)=%2d  [P]^2=1:%s\n",
    p, t, a2, sf, m, h,
    if(norm_beta == p^2, "YES", "NO!"),
    if(a2 % p == 0, "YES(!)", "no"),
    kron_sfp,
    if(ok, "YES", "NO!"));
  ok
};

\\ ============================================================
\\ Part A: 15 specific non-norm-form (p, t) pairs
\\ ============================================================

{
  my(pts, np, nt, i, p, t, r, total, passed);

  print("Thread 16: Order-2 theorem — non-norm-form primes");
  print("===================================================");
  print("Columns: p, t, a2, sf, m, h=|Cl(K)|, N(beta)=p^2?, p|a2?, (sf/p), [P]^2=1?");
  print();

  \\ Test cases: 15 (p,t) pairs with non-norm-form primes
  \\ Norm-form primes to avoid: {19,37,79,109,349,487,739,937,1279,...}
  np = [101, 103, 107, 113, 127, 131, 137, 139, 149, 151, 157, 163, 503, 1009, 10007];
  nt = [  5,   7,   3,   9,   5,  11,   7,   3,   5,   9,   7,  11,  13,  21,    7];
  total = 0; passed = 0;
  for(i = 1, #np,
    p = np[i]; t = nt[i];
    r = verify_general(p, t);
    if(r == -3, printf("  SKIP p=%d: norm-form prime\n", p); next());
    if(r < 0,   printf("  SKIP p=%d t=%d: code %d\n", p, t, r); next());
    total++;
    if(r, passed++));

  printf("\nPassed: %d / %d\n", passed, total);
}

\\ ============================================================
\\ Part B: Kronecker sweep — (sf/p) = 1 for ALL valid (p,t)?
\\ ============================================================

{
  my(kron_ok, kron_total, tmax, n, sf, kr, p, t);
  print();
  print("Part B: Kronecker sweep over all non-norm-form p in [5,200], 1<=t<2*sqrt(p)");
  print("Expected: (sf/p) = 1 for every pair.");

  kron_ok = 1; kron_total = 0;
  forprime(p = 5, 200,
    if(is_norm_form(p), next());
    tmax = floor(2*sqrt(p)) - 1;
    for(t = 1, tmax,
      n = 4*p - t^2;
      if(n <= 0, next());
      sf = sf_part(-n);
      kr = kronecker(sf, p);
      kron_total++;
      if(kr != 1,
        kron_ok = 0;
        printf("  FAIL: p=%d t=%d sf=%d (sf/p)=%d\n", p, t, sf, kr))));

  printf("Checked %d (p,t) pairs. All (sf/p)=1: %s\n", kron_total, if(kron_ok, "YES", "NO!"));
}

\\ ============================================================
\\ Part C: Symbolic confirmation of (-1/p)^2 = 1 = (sf/p)
\\ for sample primes (smallest t=1 case)
\\ ============================================================

{
  my(n, sf, kr_sf, kr_m1, p);
  print();
  print("Part C: Symbolic check (-1/p)^2 = (sf/p) for primes 5..53, t=1");

  forprime(p = 5, 53,
    n = 4*p - 1;
    sf = sf_part(-n);
    kr_sf = kronecker(sf, p);
    kr_m1 = kronecker(-1, p);
    printf("  p=%-4d (-1/p)=%-2d (-1/p)^2=%d (sf/p)=%d  match:%s\n",
      p, kr_m1, kr_m1^2, kr_sf, if(kr_m1^2 == kr_sf, "YES", "NO!")));
}

print();
print("DONE.");
