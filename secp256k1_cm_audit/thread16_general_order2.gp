\\ thread16_general_order2.gp
\\
\\ Thread 16: Verify [P]^2=1 generalises to NON-NORM-FORM primes.
\\
\\ The Thread-15 algebraic proof (A)-(E) uses only:
\\   beta = (-a2 + m*sqrt(sf))/2 in K=Q(sqrt(sf)), N(beta)=p^2, p ∤ a2.
\\ No norm-form condition 4p=73+3k^2 is ever invoked.
\\
\\ Strategy:
\\   Part A: For each non-norm-form prime p, test arbitrary a2 values.
\\   Part B: Test actual even-model Jacobians y^2=x^6+c*x^2+d (biquadratic
\\           Weil poly forced by involution (x,y)->(-x,y)).
\\
\\ NOTE on PARI quirks:
\\   - `#list` at top level is a meta-command; wrap loops in { ... } blocks.
\\   - hyperellcharpoly(poly_over_Fp) takes polynomial with Mod(1,p) coefficients.
\\   - Use length(v) instead of #v at top level.
\\
\\ Run: gp -q thread16_general_order2.gp

default(parisize, 128000000);
default(timer, 0);

\\ ---- is_norm_form: 4p = 73 + 3k^2 for some k? ----
is_norm_form(p) = {
  my(v, q);
  v = 4*p - 73;
  if (v <= 0, return(0));
  if ((v % 3) != 0, return(0));
  q = v \ 3;
  return(issquare(q))
}

\\ ---- verify_order2: returns [status, sf, m, h, split_type] ----
\\ status: 1=PASS, 0=FAIL, -1=SKIP
verify_order2(p, a2) = {
  my(D, sf, m2, m, K, h, nP, P, P2, ex, isp);
  if ((a2 % p) == 0, return([-1, 0, 0, 0, "p|a2"]));
  D = a2^2 - 4*p^2;
  if (D >= 0, return([-1, 0, 0, 0, "D>=0"]));
  if (issquare(D), return([-1, 0, 0, 0, "D sq"]));
  sf = sign(D) * core(abs(D));
  m2 = D \ sf;
  if (m2 <= 0, return([-1, 0, 0, 0, "D/sf<=0"]));
  if (!issquare(m2), return([-1, 0, 0, 0, "D/sf not sq"]));
  m = sqrtint(m2);
  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;
  nP = idealprimedec(K, p);
  if (length(nP) == 0, return([-1, sf, m, h, "no P"]));
  \\ p ramified iff p | sf; but p|sf => p|D => p|a2^2 => p|a2. Excluded above.
  if (length(nP) == 1 && nP[1][3] == 2,
    return([-1, sf, m, h, "ram(unexpected)"]));
  if (length(nP) == 1,
    \\ p inert: unique P above p with N(P)=p^2. N(beta)=p^2 => (beta)=P => P principal.
    ex = bnfisprincipal(K, nP[1])[1];
    isp = (norml2(ex) == 0);
    return([isp, sf, m, h, "inert"]));
  \\ p splits: check P^2 is principal
  P = nP[1];
  P2 = idealpow(K, P, 2);
  ex = bnfisprincipal(K, P2)[1];
  isp = (norml2(ex) == 0);
  return([isp, sf, m, h, "split"])
}

\\ ---- Build list of 10 non-norm-form primes >= 100 ----
{
  my(tp, pp);
  tp = List();
  pp = 100;
  while (length(tp) < 10,
    pp = nextprime(pp + 1);
    if (!is_norm_form(pp), listput(tp, pp)));
  tp_global = Vec(tp);
}

print("Thread 16: [P]^2=1 generality — non-norm-form primes");
print("======================================================");
print("Non-norm-form test primes: ", tp_global);
print();

\\ ==================================================================
\\ PART A: Pure algebraic test (arbitrary a2, no actual Jacobian)
\\ ==================================================================
print("PART A: Algebraic test — arbitrary a2 values, non-norm-form primes");
print("--------------------------------------------------------------------");
{
  my(total, passed, pp, a2v, a2c, res);
  total = 0; passed = 0;
  for (ii = 1, length(tp_global),
    pp = tp_global[ii];
    a2v = [pp-1, (pp\2), pp\3, pp\5, 1, -1, -(pp\5), -(pp\3), -(pp\2), -(pp-1)];
    for (jj = 1, length(a2v),
      a2c = a2v[jj];
      if (a2c == 0 || gcd(a2c, pp) != 1, next);
      res = verify_order2(pp, a2c);
      if (res[1] == -1, next);
      total++;
      if (res[1], passed++);
      printf("  p=%-7d a2=%-9d sf=%-10d m=%-6d h=%-4d [%s] => %s\n",
        pp, a2c, res[2], res[3], res[4], res[5],
        if(res[1], "PASS [P]^2=1", "FAIL"))));
  printf("PART A: %d/%d passed\n\n", passed, total)
}

\\ ==================================================================
\\ PART B: Even-model genus-2 Jacobians y^2 = x^6 + c*x^2 + d
\\  The involution (x,y)->(-x,y) forces chi(T) = T^4 + a2*T^2 + p^2 (a1=a3=0).
\\  Build poly over F_p as x^6 + Mod(c,p)*x^2 + Mod(d,p) and call hyperellcharpoly.
\\ ==================================================================
print("PART B: Even-model Jacobians y^2 = x^6 + c*x^2 + d");
print("-----------------------------------------------------");
{
  my(total, passed, pp, cc, dd, fp, chi, a1c, a2c, a3c, res, found);
  total = 0; passed = 0;
  for (ii = 1, length(tp_global),
    pp = tp_global[ii];
    printf("  p=%d:", pp);
    found = 0;
    cc = 1;
    while (cc < pp && found < 3,
      dd = 1;
      while (dd < pp && found < 3,
        \\ Polynomial over F_p: coefficients are Mod(., pp)
        fp = Mod(1,pp)*x^6 + Mod(cc,pp)*x^2 + Mod(dd,pp);
        \\ Smoothness check via resultant: disc(f) = res(f,f') mod p != 0
        if ((polresultant(lift(fp), lift(deriv(fp))) % pp) != 0,
          chi = hyperellcharpoly(fp);
          a1c = -polcoeff(lift(chi), 3);
          a2c =  polcoeff(lift(chi), 2);
          a3c = -polcoeff(lift(chi), 1);
          if (a1c == 0 && a3c == 0 && a2c != 0,
            res = verify_order2(pp, a2c);
            if (res[1] != -1,
              total++;
              if (res[1], passed++);
              printf(" [c=%d d=%d a2=%d sf=%d h=%d %s=>%s]",
                cc, dd, a2c, res[2], res[4], res[5],
                if(res[1],"PASS","FAIL"));
              found++)));
        dd++);
      cc++);
    if (found == 0, print(" (none found)"), print()));
  printf("PART B: %d/%d passed\n\n", passed, total)
}

print("Expected: all pass — algebraic proof is unconditional on norm-form.");
print("DONE.");
