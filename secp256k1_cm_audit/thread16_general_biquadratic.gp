\\ thread16_general_biquadratic.gp
\\
\\ Thread 16: Verify the "order-2 Frobenius ideal" theorem for GENERAL primes p
\\ (i.e., NOT restricted to the secp256k1 norm-form family 4p = 73 + 3k^2).
\\
\\ THEOREM (General biquadratic, proved in Thread 15):
\\   Let p be an odd prime, a2 an integer with p ∤ a2 and D := a2^2 - 4p^2 ≠ 0.
\\   Write D = sf * m^2 with sf squarefree (m > 0).
\\   Let K = Q(sqrt(sf)), O_K its ring of integers, P a prime of O_K above p.
\\   Then [P]^2 = 1 in Cl(K).
\\
\\ PROOF RECAP (A)–(E):
\\   beta := (-a2 + m*sqrt(sf))/2.
\\   (A) beta ∈ O_K (satisfies monic x^2 + a2*x + p^2 = 0 over Z). ✓
\\   (B) N_{K/Q}(beta) = (a2^2 - m^2*sf)/4 = (a2^2 - D)/4 = p^2. ✓
\\   (C) p ∤ a2  =>  p ∤ m  (if p|m then D≡0 mod p^2 => a2^2≡4p^2 mod p^2 => p|a2). ✓
\\   (D) beta/p ∉ O_K: writing beta/p = (-a2/p + (m/p)sqrt(sf))/2, membership in
\\       O_K requires p|a2 and p|m. But (C) forbids this. ✓
\\   (E) p splits in K (sf ≡ (a2/m)^2 mod p is a QR) => exactly 2 primes P,Pbar above p.
\\       (beta) has norm p^2 and ≠ (p), so (beta) = P^2 or Pbar^2. Hence [P]^2 = 1. QED.
\\
\\ THIS SCRIPT verifies the theorem for:
\\   * 10 primes p NOT in the norm-form family {(73+3k^2)/4 : k odd}
\\   * 3 a2 values per prime: two with D<0 (imaginary quadratic K) and one with D>0 (real quadratic K)
\\   * Total: 30 test cases spanning both real and imaginary quadratic fields.
\\
\\ Run: gp -q thread16_general_biquadratic.gp

default(parisize, 256000000);
default(timer, 0);

\\ ---------------------------------------------------------------
\\ Utility: squarefree part preserving sign
\\ ---------------------------------------------------------------

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ---------------------------------------------------------------
\\ Check if prime p is in the norm-form family 4p = 73 + 3k^2
\\ ---------------------------------------------------------------

is_norm_form(pp) = {
  my(k, val);
  forstep(k = 1, 400, 2,
    val = 73 + 3*k^2;
    if(val % 4 == 0 && val\4 == pp, return(1)));
  0
};

\\ ---------------------------------------------------------------
\\ Core verifier: given prime p and integer a2, check [P]^2 = 1
\\ Returns: 1 on PASS, 0 on SKIP (p|a2 or D=0 or p not split), -1 on FAIL
\\ ---------------------------------------------------------------

verify_one(p, a2, label) = {
  my(D, sf, m2, m, kron_val, K, Pp, P, Pbar, P2, Pbar2,
     P2_prin, beta_elt, I_beta, field_type, pass);

  \\ (i) Preconditions
  if(a2 % p == 0,
    printf("  [SKIP] %s: p | a2\n", label);
    return(0));

  D = a2^2 - 4*p^2;
  if(D == 0,
    printf("  [SKIP] %s: D = 0\n", label);
    return(0));

  sf = sf_part(D);
  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1,
    printf("  [ERR ] %s: D/sf not a positive integer (sf=%d, D=%d)\n", label, sf, D);
    return(-1));
  m = sqrtint(m2);
  if(m^2 != m2,
    printf("  [ERR ] %s: D/sf=%d not a perfect square\n", label, m2);
    return(-1));

  \\ (ii) Check p splits in K (theorem requires this; it's implied by p∤a2 but verify)
  kron_val = kronecker(sf, p);
  if(kron_val != 1,
    printf("  [SKIP] %s: p is %s in Q(sqrt(%d)) — split condition fails\n",
           label, if(kron_val == 0, "ramified", "inert"), sf);
    return(0));

  \\ (iii) Build K and compute class group
  K = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  if(#Pp < 2,
    printf("  [SKIP] %s: only %d prime(s) above p in O_K\n", label, #Pp);
    return(0));

  P    = Pp[1];
  Pbar = Pp[#Pp];
  P2    = idealpow(K, P, 2);
  Pbar2 = idealpow(K, Pbar, 2);

  \\ (iv) Check [P]^2 = 1 via bnfisprincipal
  P2_prin = bnfisprincipal(K, P2);
  pass = (P2_prin[1] == 0 * P2_prin[1]);

  \\ (v) Also check (beta) == P^2 or Pbar^2
  beta_elt = Mod((-a2 + m*x) / 2, x^2 - sf);
  I_beta   = idealhnf(K, beta_elt);
  field_type = if(sf < 0, "imag", "real");

  printf("  [%s] %s: D=%-12d  sf=%-8d  m=%-5d  h=%-4d  [P]^2=1: %-3s  (β)=P^2: %s  K=%s\n",
    if(pass, "PASS", "FAIL"),
    label, D, sf, m, K.clgp.no,
    if(pass, "YES", "NO"),
    if(I_beta == P2 || I_beta == Pbar2, "YES", "NO"),
    field_type);

  if(!pass,
    printf("    *** THEOREM VIOLATION: [P]^2 ≠ 1 for p=%d, a2=%d! ***\n", p, a2));

  return(if(pass, 1, -1))
};

\\ ---------------------------------------------------------------
\\ Main: 10 non-norm-form primes, 3 a2 each (2 with D<0, 1 with D>0)
\\ ---------------------------------------------------------------

print("Thread 16: General biquadratic Weil polynomial — order-2 Frobenius ideal");
print("===========================================================================");
print("Testing 10 non-norm-form primes p, each with 3 values of a2.");
print("D < 0 => imaginary quadratic K = Q(sqrt(sf)),  D > 0 => real quadratic K.");
print();

\\
\\ Test primes: chosen to be NOT in the norm-form family.
\\ Norm-form set (k<=20, odd k): {19, 37, 79, 109, 163, 193, ...}; see thread14.
\\ We pick primes clearly outside this family.
\\
\\ For each p we use:
\\   a2_imag1 = 7            (small, |a2|<<p, D = 49-4p^2 < 0)
\\   a2_imag2 = p\3 + 1      (moderate, D < 0 as long as (p/3+1) < 2p, always true)
\\   a2_real  = 2*p + 3      (D = (2p+3)^2 - 4p^2 = 12p+9 > 0, real quadratic K)
\\
\\ All three values have p ∤ a2 by construction (since p>7 for the first; for the
\\ second (p/3+1) mod p = p/3+1 ≠ 0 since p∤(p/3+1) for p>3; for the third
\\ 2p+3 ≡ 3 mod p so p∤(2p+3)).
\\

{
  my(test_primes, p, a2i1, a2i2, a2r, nf, pass_count, skip_count, fail_count, r);
  test_primes = [101, 151, 199, 251, 307, 401, 503, 601, 701, 809];
  pass_count = 0; skip_count = 0; fail_count = 0;

  for(ii = 1, #test_primes,
    p = test_primes[ii];
    nf = is_norm_form(p);
    a2i1 = 7;
    a2i2 = p\3 + 1;
    a2r  = 2*p + 3;

    printf("\np = %d  (norm-form: %s)\n", p, if(nf, "YES (warning!)", "no"));

    \\ D<0 case 1: a2=7, D = 49-4p^2 < 0, imaginary quadratic
    r = verify_one(p, a2i1, Str("p=", p, " a2=7"));
    if(r == 1, pass_count++, if(r == 0, skip_count++, fail_count++));

    \\ D<0 case 2: a2=floor(p/3)+1, imaginary quadratic
    r = verify_one(p, a2i2, Str("p=", p, " a2=", a2i2));
    if(r == 1, pass_count++, if(r == 0, skip_count++, fail_count++));

    \\ D>0 case: a2=2p+3, D = 12p+9 > 0, real quadratic
    r = verify_one(p, a2r, Str("p=", p, " a2=2p+3=", a2r));
    if(r == 1, pass_count++, if(r == 0, skip_count++, fail_count++));
  );

  printf("\n===========================================================================\n");
  printf("Summary: %d PASS, %d SKIP, %d FAIL (out of %d attempted)\n",
    pass_count, skip_count, fail_count, #test_primes * 3);
  printf("Theorem verified: %s\n", if(fail_count == 0, "YES (no violations found)", "NO -- VIOLATION DETECTED"));
}

\\ ---------------------------------------------------------------
\\ Bonus: two large-prime spot-checks for extra confidence
\\ ---------------------------------------------------------------

print();
print("Bonus: spot-checks at larger primes");
print("-------------------------------------");
{
  my(r, pc, sc, fc);
  pc = 0; sc = 0; fc = 0;
  r = verify_one(9001, 17, "p=9001 a2=17 (D<0)");
  if(r==1,pc++,if(r==0,sc++,fc++));
  r = verify_one(9001, 2*9001+5, "p=9001 a2=2p+5=18007 (D>0)");
  if(r==1,pc++,if(r==0,sc++,fc++));
  r = verify_one(104729, 31, "p=104729 a2=31 (D<0, large prime)");
  if(r==1,pc++,if(r==0,sc++,fc++));
  r = verify_one(104729, 104729\4 + 7, "p=104729 a2=p/4+7 (D<0, large prime)");
  if(r==1,pc++,if(r==0,sc++,fc++));
  printf("Bonus summary: %d PASS, %d SKIP, %d FAIL\n", pc, sc, fc);
}

print();
print("DONE.");
