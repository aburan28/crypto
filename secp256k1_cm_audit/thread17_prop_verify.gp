\\ thread17_prop_verify.gp
\\
\\ Thread 17: Verify Proposition prop:order2-frobenius from structural_completeness.tex.
\\ For odd prime p and integer a2 with p!|a2, a2 != 0:
\\   D = a2^2 - 4p^2 = sf*m^2  (sf squarefree, K=Q(sqrt(sf)))
\\   (i)  m^2*sf - a2^2 + 4p^2 = 0            [N(beta)=p^2]
\\   (ii) [P]^2 = 1 in Cl(K)
\\   (iii) p splits in K (not inert, not ramified)
\\
\\ Following thread16 PARI idioms exactly.
\\ Run: gp -q < secp256k1_cm_audit/thread17_prop_verify.gp

default(parisize, 128000000);
default(timer, 0);

\\ Squarefree part via built-in core()
sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ Verify the proposition for (p, a2).
\\ Returns 1=pass, 0=fail, -1=skip.
verify_prop(p, a2) = {
  my(D, sf, m2, m, K, Pp, P2, cls, result);
  if (a2 == 0 || a2 % p == 0, return(-1));
  D = a2^2 - 4*p^2;
  if (D == 0, return(-1));
  sf = sf_part(D);
  if (sf == 1, return(-1));  \\ D a perfect square => K=Q, not quadratic
  m2 = D / sf;
  if (denominator(m2) != 1 || m2 <= 0, return(-1));
  m = sqrtint(m2);
  if (m^2 != m2, return(-1));
  \\ Part (i): algebraic identity m^2*sf - a2^2 + 4p^2 = D - D = 0
  if (m2*sf - a2^2 + 4*p^2 != 0,
    printf("FAIL(i) p=%d a2=%d\n", p, a2); return(0));
  \\ Build K = Q(sqrt(sf))
  K = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  \\ Part (iii): p must split => exactly 2 primes above p, each with f=1
  if (#Pp != 2 || Pp[1].f != 1 || Pp[2].f != 1,
    printf("FAIL(iii) p=%d a2=%d nPp=%d\n", p, a2, #Pp); return(0));
  \\ Part (ii): P^2 is principal
  P2 = idealpow(K, Pp[1], 2);
  cls = bnfisprincipal(K, P2);
  result = (cls[1] == 0*cls[1]);
  if (!result,
    printf("FAIL(ii) [P]^2 not principal p=%d a2=%d\n", p, a2); return(0));
  return(1);
}

\\====================================================================
\\ Battery I: 20 pairs D<0 (imaginary quadratic K)
\\====================================================================
print("Battery I: D < 0 (imaginary quadratic K)");
{
  my(total=0, ok=0, r);
  my(cases = [
    [3,1], [5,2], [7,3], [11,4], [13,5],
    [17,6], [19,7], [23,8], [29,9], [31,10],
    [37,11], [41,12], [43,13], [47,14], [53,15],
    [59,16], [61,17], [67,18], [71,19], [73,20]
  ]);
  for(i=1, #cases,
    my(p=cases[i][1], a2=cases[i][2]);
    r = verify_prop(p, a2);
    if (r >= 0, total++; ok += r);
  );
  printf("Battery I: %d/%d passed\n\n", ok, total);
}

\\====================================================================
\\ Battery II: 10 pairs D>0 (real quadratic K, |a2| > 2p)
\\====================================================================
print("Battery II: D > 0 (real quadratic K)");
{
  my(total=0, ok=0, r);
  my(cases = [
    [3,7], [5,11], [7,15], [11,23], [13,27],
    [17,35], [19,40], [23,48], [29,60], [31,64]
  ]);
  for(i=1, #cases,
    my(p=cases[i][1], a2=cases[i][2]);
    r = verify_prop(p, a2);
    if (r >= 0, total++; ok += r);
  );
  printf("Battery II: %d/%d passed\n\n", ok, total);
}

\\====================================================================
\\ Battery III: mass sweep — all (p, a2) with p prime <= 53, 1 <= a2 < p
\\====================================================================
print("Battery III: mass sweep p<=53, 1<=a2<p");
{
  my(total=0, ok=0, r);
  forprime(p=3, 53,
    for(a2=1, p-1,
      r = verify_prop(p, a2);
      if (r >= 0, total++; ok += r);
    );
  );
  printf("Battery III: %d/%d passed\n\n", ok, total);
}

\\====================================================================
\\ Summary
\\====================================================================
print("=== SUMMARY ===");
print("Proposition prop:order2-frobenius: all batteries completed.");
print("Check output above for any FAIL lines.");
quit
