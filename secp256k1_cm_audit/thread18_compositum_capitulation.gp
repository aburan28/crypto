/* thread18_compositum_capitulation.gp  (v2 — fixed comparison logic)
 *
 * Thread 18 — Frobenius ideal capitulation in biquadratic composita
 *
 * Question: For biquadratic Weil poly W = f1*f2 (f_i = x^2 + t_i*x + p),
 * Prop (Thread 17) gives [P_i]^2 = 1 in Cl(K_i).
 * The Jacobian's CM field is the compositum K = Q(sqrt(sf1), sqrt(sf2)).
 * Does [P] capitalise to order 1 in Cl(K), or stay order 2?
 *
 * Part D compares with irreducible quartic palindromes (non-biquadratic
 * Weil polys), showing that order > 2 IS possible in that case,
 * confirming the order-2 property is specific to the biquadratic type.
 */

/* Helper: order of ideal I in Cl(K) using SNF invariant factors */
idealcl_order(K, I) = {
  local(v, cyc, ord, d, g);
  v = bnfisprincipal(K, I); /* returns [exponents, alpha] */
  cyc = K.clgp.cyc;
  if(#cyc == 0, return(1));   /* trivial class group */
  /* v[1] = exponent vector [e1,...,ek], cyc = [d1,...,dk] */
  /* ord([I]) = lcm_{i}( d_i / gcd(e_i, d_i) ) */
  ord = 1;
  for(i=1, min(#cyc, #v[1]),
    d = cyc[i];
    g = gcd(v[1][i], d);
    ord = lcm(ord, d/g)
  );
  return(ord)
}

{
printf("=== Thread 18: Compositum capitulation check ===\n\n");

/* -----------------------------------------------------------
 * Part A: Exhaustive biquadratic order check for p=113
 * Tests ALL irreducible x^2+ax+p: expect ord([P]) ≤ 2.
 * ----------------------------------------------------------- */
printf("Part A: Exhaustive biquadratic check (p=113)\n");
p = 113;
sq = floor(2*sqrt(p+0.0));
count_a = 0; max_ord_a = 1; bad_a = 0;
for(a = -sq, sq,
  f1 = Pol([1, a, p]);
  if(!polisirreducible(f1), next);
  K1 = bnfinit(f1, 1);
  P1 = idealprimedec(K1, p)[1];
  ord_a = idealcl_order(K1, P1);
  count_a++;
  if(ord_a > max_ord_a, max_ord_a = ord_a);
  if(ord_a > 2, bad_a++;
     printf("  COUNTEREX: a=%d, ord=%d (h=%d)\n", a, ord_a, K1.clgp.no))
);
printf("  %d quadratics checked, max ord([P]) = %d", count_a, max_ord_a);
if(bad_a, printf(" *** %d COUNTEREXAMPLES ***", bad_a),
          printf(" [CONFIRMED: all ≤ 2 as per Prop 4.X]"));
printf("\n\n");

/* -----------------------------------------------------------
 * Part B: Compositum Q(sqrt(sf1),sqrt(sf2)) — capitulation check
 * For each Howe-glueable-type pair (t1,t2) over F_{p2},
 * compute ord([P]) in Cl(compositum).
 * -----------------------------------------------------------*/
printf("Part B: Compositum Q(sqrt(sf1),sqrt(sf2)) — capitulation check (p=89)\n");
p2 = 89;
sq2 = floor(2*sqrt(p2+0.0));
comp_count = 0; order1_n = 0; order2_n = 0; orderN_n = 0;

for(t1 = -sq2, sq2,
  f1 = Pol([1, t1, p2]);
  if(!polisirreducible(f1), next);
  sf1 = core(t1^2 - 4*p2);
  for(t2 = t1+1, sq2,
    f2 = Pol([1, t2, p2]);
    if(!polisirreducible(f2), next);
    sf2 = core(t2^2 - 4*p2);
    if(sf1 == sf2, next);  /* same CM field — isogenous pair */

    clist = polcompositum(f1, f2);
    if(#clist == 0, next);
    fc = clist[1];
    Kc = bnfinit(fc, 1);
    PP = idealprimedec(Kc, p2);

    for(j = 1, #PP,
      P = PP[j];
      ord_c = idealcl_order(Kc, P);
      comp_count++;
      if(ord_c == 1, order1_n++,
         ord_c == 2, order2_n++,
         orderN_n++)
    );
    if(comp_count >= 300, break(2))
  )
);
printf("  %d prime ideals in composita checked:\n", comp_count);
printf("    ord=1 (capitulated):   %d  (%.1f%%)\n",
       order1_n, 100.0*order1_n/max(comp_count,1));
printf("    ord=2 (stays order-2): %d  (%.1f%%)\n",
       order2_n, 100.0*order2_n/max(comp_count,1));
printf("    ord>2:                 %d\n\n", orderN_n);

/* -----------------------------------------------------------
 * Part C: Explicit examples showing ord in each layer
 * For pairs where both K_i have non-trivial class group,
 * print the order in K_i and in the compositum.
 * ----------------------------------------------------------- */
printf("Part C: Explicit order table (p=71)\n");
printf("  %s | %s | %s\n",
       "t1,t2     sf1,sf2      h(K1),h(K2),h(Kc)",
       "ord_K1  ord_K2",
       "ords_Kc");
p3 = 71;
sq3 = floor(2*sqrt(p3+0.0));
shown = 0;
for(t1 = -sq3, sq3,
  if(shown >= 12, break);
  f1 = Pol([1, t1, p3]);
  if(!polisirreducible(f1), next);
  sf1 = core(t1^2 - 4*p3);
  K1 = bnfinit(f1, 1);
  h1 = K1.clgp.no;
  for(t2 = t1+1, sq3,
    if(shown >= 12, break);
    f2 = Pol([1, t2, p3]);
    if(!polisirreducible(f2), next);
    sf2 = core(t2^2 - 4*p3);
    if(sf1 == sf2, next);
    K2 = bnfinit(f2, 1);
    h2 = K2.clgp.no;

    clist = polcompositum(f1, f2);
    if(#clist == 0, next);
    fc = clist[1];
    Kc = bnfinit(fc, 1);
    hc = Kc.clgp.no;

    P1 = idealprimedec(K1, p3)[1];
    P2 = idealprimedec(K2, p3)[1];
    o1 = idealcl_order(K1, P1);
    o2 = idealcl_order(K2, P2);

    PP_c = idealprimedec(Kc, p3);
    ords_c = "";
    for(j=1,#PP_c,
      oc = idealcl_order(Kc, PP_c[j]);
      ords_c = Str(ords_c, oc, if(j<#PP_c,",",""))
    );

    printf("  t1=%3d,t2=%3d | sf=%4d,%5d | h=%d,%d,%3d | o1=%d,o2=%d | [%s]\n",
           t1, t2, sf1, sf2, h1, h2, hc, o1, o2, ords_c);
    shown++
  )
);

/* -----------------------------------------------------------
 * Part D: Irreducible quartic palindromic polys — can ord > 2?
 * f = x^4 + c3*x^3 + c2*x^2 + c3*x + 1 (palindromic, totally imaginary)
 * We expect to find class group elements of order > 2 above some prime q,
 * showing the order-2 property is NOT universal for all CM quartic fields.
 * ----------------------------------------------------------- */
printf("\nPart D: Irreducible quartic palindromes — ord > 2 examples\n");
found_d = 0;
for(c3 = -5, 5,
  for(c2 = -8, 20,
    fq = Pol([1, c3, c2, c3, 1]);
    if(!polisirreducible(fq), next);
    /* Check totally imaginary: all 4 roots complex */
    rr = polroots(fq * 1.0);
    ok = 1;
    for(j=1,4, if(abs(imag(rr[j])) < 1e-6, ok=0; break));
    if(!ok, next);
    Kq = bnfinit(fq, 1);
    hq = Kq.clgp.no;
    if(hq == 1, next);
    cyc_q = Kq.clgp.cyc;
    /* Max possible order of any element */
    max_possible = vecmax(cyc_q);
    if(max_possible <= 2, next);  /* class group has no order-4+ elements */
    /* Find a prime q where ord([P_q]) > 2 */
    for(q = 2, 40,
      if(!isprime(q), next);
      PPq = idealprimedec(Kq, q);
      for(j=1,#PPq,
        ord_q = idealcl_order(Kq, PPq[j]);
        if(ord_q > 2,
          printf("  f=x^4+%d*x^3+%d*x^2+%d*x+1, q=%d, h=%d, cyc=%s, P_%d: ord=%d\n",
                 c3, c2, c3, q, hq, Str(cyc_q), j, ord_q);
          found_d++;
          if(found_d >= 10, break(4))
        )
      )
    )
  )
);
if(!found_d, printf("  No order>2 found — widen search range\n"));

printf("\n=== Done ===\n");
}
