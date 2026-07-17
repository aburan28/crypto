\\ thread15_order2_verify.gp
\\ Thread 15: Verify the "universal order-2" conjecture for norm-form Frobenius.
\\
\\ Claim: for every prime p in the norm-form family 4p = 73+3k^2,
\\   let K = Q(sqrt(sf)) where sf = squarefree-part(a2^2 - 4p^2),
\\   and let Pp be a prime of K above p.  Then
\\
\\     Pp^2 = (pi_sq)  as ideals in O_K,
\\
\\   where  pi_sq := (-a2 + m*sqrt(sf))/2  and  m = sqrt((a2^2-4p^2)/sf).
\\
\\ Key: Nm_K(pi_sq) = (a2^2-disc4)/4 = p^2, so (pi_sq) has norm p^2 = Nm(Pp^2).
\\
\\ All 25 norm-form primes from Thread 14 (k<=199) are tested.
\\ Run:  gp --stacksize 128000000 -q thread15_order2_verify.gp

default(parisize, 128000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

check_one(k, p, a2, sf) = {
  my(disc4, sf_chk, ratio, m, K, el, I_pi2, n_pi2);
  my(Pp_list, Pp, Pp2, D, n_D, res, ce, cyc, j, all_zero, ok_this, tag, h);

  disc4 = a2^2 - 4*p^2;

  sf_chk = sf_part(disc4);
  if(sf_chk != sf,
    printf("  k=%-3d WARN sf mismatch: computed %d, table %d\n",k,sf_chk,sf);
    return(0)
  );

  ratio = disc4 / sf;
  m = sqrtint(ratio);
  if(m^2 != ratio,
    printf("  k=%-3d ERROR: disc4/sf=%d not a perfect square\n",k,ratio);
    return(0)
  );

  K = bnfinit(x^2 - sf, 1);
  h = K.clgp.no;

  \\ Element pi_sq = (-a2 + m*sqrt(sf))/2 in K  (K.pol = x^2-sf, x=sqrt(sf))
  el = Mod((-a2 + m*x)/2, K.pol);

  I_pi2 = idealhnf(K, el);
  n_pi2 = idealnorm(K, I_pi2);

  Pp_list = idealprimedec(K, p);
  Pp = Pp_list[1];
  Pp2 = idealmul(K, Pp, Pp);

  \\ Equality check: idealdiv should give unit ideal (norm 1) if equal
  D = idealdiv(K, I_pi2, Pp2);
  n_D = idealnorm(K, D);
  \\ If not matching Pp, try conjugate Pp_bar (split case, 2 primes above p)
  if(n_D != 1 && #Pp_list >= 2,
    my(Pp2c);
    Pp2c = idealmul(K, Pp_list[2], Pp_list[2]);
    D = idealdiv(K, I_pi2, Pp2c);
    n_D = idealnorm(K, D)
  );

  res = bnfisprincipal(K, Pp2);
  ce = res[1];
  cyc = K.clgp.cyc;
  all_zero = 1;
  for(j = 1, #cyc, if(ce[j] != 0, all_zero = 0; break));

  ok_this = (n_pi2 == p^2) && (n_D == 1) && all_zero;
  if(ok_this, tag = "OK", tag = "FAIL");

  printf("  k=%-3d p=%-7d sf=%-9d m=%-5d h=%-4d | norm_ok=%d eq=%d princ=%d [%s]\n",
    k, p, sf, m, h, (n_pi2 == p^2), (n_D == 1), all_zero, tag);
  ok_this
};

\\ Data from Thread 14: [k, p, a2, sf]  (single-line — multi-line fails in gp 2.15)
nf_data = [[1,19,-35,-219],[5,37,1,-219],[9,79,85,-219],[11,109,145,-219],[21,349,385,-939],[25,487,-299,-3819],[31,739,-1403,-8643],[35,937,1,-5619],[41,1279,-2315,-14619],[55,2287,-899,-16419],[65,3187,-4499,-32619],[85,5437,-9791,-61995],[91,6229,-11375,-71499],[99,7369,385,-43059],[101,7669,-13751,-87267],[105,8287,2149,-1731],[109,8929,5905,-35859],[119,10639,10549,-32187],[131,12889,-25751,-3],[135,13687,5701,-65019],[141,14929,-15575,-136299],[145,15787,-4499,-108219],[159,18979,-32915,-212619],[179,24049,-32975,-243219],[195,28537,-57071,-342435]];

print("Thread 15: Verify (pi_sq)=Pp^2 for norm-form primes 4p=73+3k^2");
print("pi_sq = (-a2 + m*sqrt(sf))/2,  m^2 = (a2^2-4p^2)/sf");
print("===================================================================");

ok_count = 0;
fail_count = 0;
for(idx = 1, #nf_data, {
  my(row, res_ok);
  row = nf_data[idx];
  res_ok = check_one(row[1], row[2], row[3], row[4]);
  if(res_ok, ok_count++, fail_count++)
});

printf("\nSummary: %d OK, %d FAIL out of %d cases.\n",
       ok_count, fail_count, #nf_data);

\\ ---- Symbolic norm check (a2^2-disc4)/4 = p^2 ----
print();
print("Algebraic norm: Nm(pi_sq) = (a2^2 - disc4)/4 = p^2  [first 4 cases]");
for(idx = 1, min(4, #nf_data), {
  my(row, k, p, a2, sf, disc4, ratio, m, nm);
  row = nf_data[idx];
  k = row[1]; p = row[2]; a2 = row[3]; sf = row[4];
  disc4 = a2^2 - 4*p^2;
  ratio = disc4/sf;
  m = sqrtint(ratio);
  nm = (a2^2 - ratio*sf) / 4;
  printf("  k=%-3d p=%-5d: (a2^2-disc4)/4 = %d (expect %d) [%s]\n",
    k, p, nm, p^2, if(nm == p^2, "OK", "FAIL"))
});

\\ ---- Integrality: minpoly of pi_sq has Z coefficients ----
print();
print("Integrality: minpoly of pi_sq should lie in Z[x]  [first 4 cases]");
for(idx = 1, min(4, #nf_data), {
  my(row, k, p, a2, sf, disc4, ratio, m, K, el, mp, dn);
  row = nf_data[idx];
  k = row[1]; p = row[2]; a2 = row[3]; sf = row[4];
  disc4 = a2^2 - 4*p^2;
  ratio = disc4/sf;
  m = sqrtint(ratio);
  K = nfinit(x^2 - sf);
  el = Mod((-a2 + m*x)/2, K.pol);
  mp = minpoly(el);
  dn = denominator(mp);
  printf("  k=%-3d: minpoly = %Ps  denom=%d [%s]\n",
    k, mp, dn, if(dn == 1, "integral", "NOT integral"))
});

\\ ---- h=1 case: sf=-3, k=131, p=12889 ----
print();
print("h=1 special case: k=131, p=12889, sf=-3 [Q(sqrt(-3)), Eisenstein]");
{
  my(K3, Pp3, r1, r2);
  K3 = bnfinit(x^2 + 3, 1);
  Pp3 = idealprimedec(K3, 12889)[1];
  r1 = bnfisprincipal(K3, Pp3);
  r2 = bnfisprincipal(K3, idealmul(K3, Pp3, Pp3));
  printf("  Pp  class=[%Ps], Nm=%d\n", r1[1], idealnorm(K3, Pp3));
  printf("  Pp2 class=[%Ps], Nm=%d\n", r2[1], idealnorm(K3, idealmul(K3,Pp3,Pp3)))
}

print();
print("DONE.");
quit
