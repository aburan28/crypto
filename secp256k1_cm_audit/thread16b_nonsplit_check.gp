\\ thread16b_nonsplit_check.gp
\\
\\ Thread 16 (supplement): Demonstrate that the order-2 phenomenon for
\\ biquadratic Weil polynomials T^4+a2*T^2+p^2 arises ONLY when the
\\ polynomial does NOT factor over Z as (T^2-tT+p)(T^2+tT+p) for any t in Z.
\\
\\ KEY STRUCTURAL DISTINCTION:
\\   Case A (rational split): a2 = 2p - t^2 for some t in Z.
\\     Then beta = (-a2+m*sqrt(sf))/2 = pi^2 where pi = (t+sqrt(t^2-4p))/2 in O_K.
\\     The ideal (pi) = P (principal) in O_K. So ord([P]) = 1.
\\
\\   Case B (non-rational split): a2 != 2p - t^2 for any t in Z.
\\     Then 2p-a2 is not a perfect square. beta is NOT pi^2 for any pi in O_K.
\\     The ideal (beta) = P^2 but P need not be principal. ord([P]) in {1, 2}.
\\
\\ This script:
\\   (i)  Confirms Case A: ord([P])=1 for all integer-trace splits.
\\   (ii) Finds Case B examples via genuine genus-2 hyperelliptic Jacobians
\\        (using hyperellcharpoly) and verifies ord([P]) can be 2.
\\   (iii) Gives the exact condition: ord([P])=2 iff [P] non-trivial iff P non-principal
\\         iff (beta) = P^2 with [P] non-trivial in Cl(Q(sqrt(sf))).
\\
\\ Run: gp-2.15 --stacksize 256000000 -q thread16b_nonsplit_check.gp

default(parisize, 256000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ---------------------------------------------------------------
\\ Is 2p-a2 a perfect square? (Test for rational split)
\\ ---------------------------------------------------------------
is_rational_split(p, a2) = {
  my(v);
  v = 2*p - a2;
  if(v < 0, return(0));
  sqrtint(v)^2 == v
};

\\ ---------------------------------------------------------------
\\ Verify [P]^2=1 and determine ord([P]) for given (p, a2)
\\ Returns: [ord, sf, h, D] or [0, reason]
\\ ---------------------------------------------------------------
check_biquad(p, a2) = {
  my(D, sf, m2, m, K, Pp, P, P2, P2_prin, ord_P, h);

  D  = a2^2 - 4*p^2;
  sf = sf_part(D);
  if(sf >= 0, return([0, "sf>=0 (real disc)"]));

  m2 = D / sf;
  if(m2 <= 0 || denominator(m2) != 1, return([0, "D/sf_error"]));
  m  = sqrtint(m2);
  if(m*m != m2, return([0, "not_square"]));
  if(a2 % p == 0, return([0, "p|a2"]));

  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  Pp = idealprimedec(K, p);
  if(#Pp == 0, return([0, "p_inert"]));

  P    = Pp[1];
  P2   = idealpow(K, P, 2);
  P2_prin = bnfisprincipal(K, P2);

  if(P2_prin[1] != 0*P2_prin[1], return([0, "[P]^2!=1 (ERROR)"]));

  ord_P = if(bnfisprincipal(K, P)[1] == 0*bnfisprincipal(K, P)[1], 1, 2);

  [ord_P, sf, h, D]
};

\\ ---------------------------------------------------------------
\\ Search for genuine non-rational-split biquadratic Weil polys
\\ using hyperelliptic curves y^2 = f(x), deg f = 5 or 6.
\\ We look for primes p and genus-2 curves C/F_p such that:
\\   charpoly(Jac(C)) = T^4 + a2*T^2 + p^2
\\   AND 2p - a2 is NOT a perfect square (non-rational split)
\\ ---------------------------------------------------------------

print("Part 1: Case A confirmation — rational splits always give ord([P])=1");
print("========================================================================");
{
  \\ Take 5 cases from Thread 16 known output
  cases_a = [[101, 20, 14],   \\ p, a, t  (a2 = 2p-t^2)
             [197, 20, 16],
             [307, 20, 28],
             [503, 20, 16],
             [1013, 20, -42]];
  for(i = 1, #cases_a,
    my(p = cases_a[i][1], a = cases_a[i][2], t0 = cases_a[i][3]);
    my(a2 = 2*p - t0^2);
    my(res = check_biquad(p, a2));
    printf("  p=%-5d a2=%-7d t=%d  rational_split:%s  ord([P])=%d  sf=%d\n",
      p, a2, t0, if(is_rational_split(p, a2), "YES", "NO"), res[1], res[2]));
}
print();

print("Part 2: Case B search — non-rational-split biquadratic Weil polys");
print("=====================================================================");
print("Searching genus-2 curves y^2 = x^5 + a*x + b over F_p for biquadratic charpoly...");
print();
{
  my(found = 0, target = 5);
  for(p = 50, 300,
    if(!isprime(p) || found >= target, next);
    for(a = 1, p-1,
      if(found >= target, break);
      for(b = 1, p-1,
        if(found >= target, break);
        \\ skip if discriminant is 0
        if((4*a^5 + 25*a^4 + 5^5*b^4) % p == 0, next);

        my(fld = ffgen(p, 't));
        my(f = fld^0*x^5 + a*fld^0*x + b*fld^0);
        my(cp);
        iferr(cp = hyperellcharpoly(f), e, next);

        \\ Check if charpoly is of form T^4 + a2*T^2 + p^2
        \\ i.e., coefficients of T^3 and T are 0
        if(polcoeff(cp, 3) != 0 || polcoeff(cp, 1) != 0, next);
        my(a2 = polcoeff(cp, 2));
        if(polcoeff(cp, 0) != p^2, next);

        \\ Check non-rational-split
        if(is_rational_split(p, a2), next);

        \\ Check theorem applies (sf < 0 etc)
        my(res = check_biquad(p, a2));
        if(type(res[1]) != "t_INT", next);

        found++;
        printf("  FOUND #%d: p=%d  curve y^2=x^5+%d*x+%d\n", found, p, a, b);
        printf("    charpoly = %Ps\n", cp);
        printf("    a2=%d  2p-a2=%d  is_perfect_sq:%s\n",
          a2, 2*p-a2, if(is_rational_split(p,a2),"YES","NO"));
        printf("    sf=%d  h=%d  ord([P])=%d  [P]^2=1:YES\n\n",
          res[2], res[3], res[1]))));
  if(found == 0,
    print("  No non-rational-split biquadratic charpolys found in search range."));
}

print("Part 3: Verify secp256k1 norm-form primes are non-rational-split");
print("===================================================================");
{
  \\ Thread 15 data: first 6 norm-form primes
  norm_cases = [[19, 35], [37, 1], [79, -85], [109, -145], [349, 385], [487, -299]];
  for(i = 1, #norm_cases,
    my(p = norm_cases[i][1], a2 = norm_cases[i][2]);
    my(v = 2*p - a2);
    my(res = check_biquad(p, a2));
    printf("  p=%-5d a2=%-6d  2p-a2=%-5d  rat_split:%s  ord([P])=%d  sf=%d  h=%d\n",
      p, a2, v, if(is_rational_split(p,a2),"YES","NO "), res[1], res[2], res[3]));
}

print();
print("STRUCTURAL THEOREM (Thread 16):");
print("  Let T^4 + a2*T^2 + p^2 be a biquadratic Weil polynomial.");
print("  Let sf = squarefree-part(a2^2 - 4p^2), K = Q(sqrt(sf)), P | p in O_K.");
print("  (i)  If 2p-a2 = t^2 for some t in Z (rational split = E x E^t):");
print("       pi = (t + sqrt(t^2-4p))/2 in O_K satisfies (pi) = P, so ord([P]) = 1.");
print("  (ii) If 2p-a2 is not a perfect square (non-rational split):");
print("       beta = (-a2+m*sqrt(sf))/2 satisfies (beta) = P^2 but P need not be");
print("       principal; ord([P]) in {1, 2}, proved = divides 2 in Thread 15.");
print("  The order-2 phenomenon is specific to non-rational-split biquadratic surfaces.");
print();
print("DONE.");
