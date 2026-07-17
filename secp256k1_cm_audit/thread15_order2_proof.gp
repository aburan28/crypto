\\ thread15_order2_proof.gp — Thread 15
\\ Verify P^2 is principal; compute generators; test genus-theory explanation.
\\ Run: gp --stacksize 128000000 -q thread15_order2_proof.gp

default(parisize, 128000000);
default(timer, 0);

\\ ---------------------------------------------------------------
\\ cl_order(ce, cyc): order of class element in Cl(K)
\\ ---------------------------------------------------------------
cl_order(ce, cyc) = {
  my(n, lce, maxord, ord, tmp, z, j);
  n = #cyc;
  if(n == 0, return(1));
  lce = vector(n, j, lift(ce[j]));
  z = 1;
  for(j = 1, n, if(lce[j] % cyc[j] != 0, z = 0; break));
  if(z, return(1));
  maxord = lcm(Vec(cyc));
  for(ord = 2, maxord,
    tmp = vector(n, j, lift(Mod(lce[j] * ord, cyc[j])));
    z = 1;
    for(j = 1, n, if(tmp[j] != 0, z = 0; break));
    if(z, return(ord))
  );
  -1
};

\\ ---------------------------------------------------------------
\\ quad_norm(disc_K, a, b): Norm of element [a,b]~ in Q(sqrt(sf))
\\ Uses integral-basis formula depending on disc mod 4.
\\ If disc % 4 == 0: zk=[1,x], Nm = a^2 - (disc/4)*b^2
\\ If disc % 4 != 0 (disc ≡ 1 mod 4): zk=[1,(x-1)/2], Nm = a^2 - a*b + (1-disc)/4*b^2
\\ ---------------------------------------------------------------
quad_norm(disc_K, a, b) = {
  if(disc_K % 4 == 0,
    a^2 - (disc_K \ 4) * b^2
  ,
    a^2 - a*b + (1 - disc_K) \ 4 * b^2
  )
};

\\ ---------------------------------------------------------------
\\ prime_discs(sf): prime discriminants dividing disc(Q(sqrt(sf)))
\\ ---------------------------------------------------------------
prime_discs(sf) = {
  my(D, rem, res, q, pd);
  D   = quaddisc(sf);
  res = [];
  if(D % 2 == 0,
    if(D % 8 == 0,
      res = concat(res, [if(D > 0, 8, -8)])
    ,
      res = concat(res, [if(D < 0, -4, 4)])
    )
  );
  rem = abs(D);
  while(rem % 2 == 0, rem = rem \ 2);
  q = 3;
  while(q * q <= rem,
    if(rem % q == 0,
      pd = if((q - 1) % 4 == 0, q, -q);
      res = concat(res, [pd]);
      while(rem % q == 0, rem = rem \ q)
    );
    q += 2
  );
  if(rem > 1,
    pd = if((rem - 1) % 4 == 0, rem, -rem);
    res = concat(res, [pd])
  );
  res
};

\\ ---------------------------------------------------------------
\\ in_pg(sf, p): 1 if all genus chars (pd_i / p) equal 1
\\ ---------------------------------------------------------------
in_pg(sf, p) = {
  my(pds, ok, j);
  pds = prime_discs(sf);
  ok  = 1;
  for(j = 1, #pds, if(kronecker(pds[j], p) != 1, ok = 0; break));
  ok
};

\\ ---------------------------------------------------------------
\\ Cases from Thread 14  (matched k, p, sf)
\\ ---------------------------------------------------------------
K_sf = [-219,-939,-1731,-3,-3819,-5619,-8643,-32187,-32619,-35859,-43059,-14619,-16419,-61995,-71499,-87267];
K_k  = [1,21,105,131,25,35,31,119,65,109,99,41,55,85,91,101];
K_p  = [19,349,8287,12889,487,937,739,10639,3187,8929,7369,1279,2287,5437,6229,7669];
NC   = #K_sf;

print("=== Thread 15: P^2 Principality + Genus Analysis ===");
print();

\\ ---------------------------------------------------------------
\\ Part A+B: P^2 principal check + explicit generator norm
\\ ---------------------------------------------------------------
print("--- Part A+B: P^2 principal check and generator norm ---");
print("  sf            k    p         h    ord(P) P2ok  pg   gen_norm    p^2");

{
  my(i, sf, kk, p, K, h, cyc, Pl, P, P2, r1, r2, ce1, ce2, ord1, gen2, gnorm, pg, disc_K, allz, j);
  for(i = 1, NC,
    sf = K_sf[i]; kk = K_k[i]; p = K_p[i];
    K  = bnfinit(x^2 - sf, 1);
    h  = K.clgp.no;
    cyc = K.clgp.cyc;
    disc_K = K.disc;
    Pl = idealprimedec(K, p);
    if(#Pl == 0,
      printf("  %-12d  %-4d  %-8d  INERT\n", sf, kk, p); next
    );
    P   = Pl[1];
    P2  = idealpow(K, P, 2);
    r1  = bnfisprincipal(K, P);
    ce1 = r1[1];
    ord1 = cl_order(ce1, cyc);
    r2  = bnfisprincipal(K, P2, 1);
    ce2 = r2[1];
    allz = 1;
    for(j = 1, #ce2, if(lift(ce2[j]) != 0, allz = 0; break));
    gen2 = r2[2];
    gnorm = quad_norm(disc_K, gen2[1], gen2[2]);
    pg = in_pg(sf, p);
    printf("  %-12d  %-4d  %-8d  %-5d  %-6d %-5d %-4d %-11d %d\n",
           sf, kk, p, h, ord1, allz, pg, gnorm, p^2)
  )
}

print();

\\ ---------------------------------------------------------------
\\ Part C: Genus characters per case
\\ ---------------------------------------------------------------
print("--- Part C: Genus characters per case ---");
print("  (All=1 → principal genus; multiple -1 → higher genus, still ord|2 holds)");

{
  my(i, sf, kk, p, pds, chars, j);
  for(i = 1, NC,
    sf = K_sf[i]; kk = K_k[i]; p = K_p[i];
    pds   = prime_discs(sf);
    chars = vector(#pds, j, kronecker(pds[j], p));
    printf("  sf=%-12d  p=%-8d  disc=%-8d  pdiscs=%Ps  chars=%Ps\n",
           sf, p, quaddisc(sf), pds, chars)
  )
}

print();

\\ ---------------------------------------------------------------
\\ Part D: Explicit [a,b] coordinates of generator of P^2
\\ ---------------------------------------------------------------
print("--- Part D: Generator beta of P^2 (coordinates in nf basis [1,omega]) ---");

{
  my(i, sf, kk, p, K, Pl, P, P2, r2, gen2, disc_K, gnorm);
  for(i = 1, NC,
    sf = K_sf[i]; kk = K_k[i]; p = K_p[i];
    K  = bnfinit(x^2 - sf, 1);
    disc_K = K.disc;
    Pl = idealprimedec(K, p);
    if(#Pl == 0, printf("  sf=%-12d p=%-8d INERT\n", sf, p); next);
    P   = Pl[1];
    P2  = idealpow(K, P, 2);
    r2  = bnfisprincipal(K, P2, 1);
    gen2  = r2[2];
    gnorm = quad_norm(disc_K, gen2[1], gen2[2]);
    \\ Print: sf, p, [a,b], norm, expected p^2
    printf("  sf=%-12d  p=%-8d  beta=[%d,%d]  Nm(beta)=%d  p^2=%d  match=%d\n",
           sf, p, gen2[1], gen2[2], gnorm, p^2, (gnorm==p^2))
  )
}

print();
print("--- Summary ---");
print("KEY FINDING from Part C: genus chars NOT all 1 for cases:");
print("  sf=-3819,-8643,-43059,-61995,-87267  (chars contain -1 pairs).");
print("CONCLUSION: genus-theory principal-genus membership does NOT explain ord([P])|2.");
print("OPEN: Why does [P]^2=0 hold even outside the principal genus?");
print("PROPOSED Thread 16: Investigate degree-4 CM field structure of the Weil poly.");
print("  Specifically: is [pi]^2=0 in Cl(L) for the degree-4 CM field L of A/Fp?");
print("  If so, the restriction map Cl(L)->Cl(Q(sqrt(sf))) maps [pi_L]^2 to [P_sf]^2=0.");
print();
print("DONE.");
quit()
