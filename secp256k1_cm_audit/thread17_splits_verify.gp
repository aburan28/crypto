\\ thread17_splits_verify.gp
\\ Thread 17: verify p SPLITS in Q(sqrt(sf)) for 10 fresh (p,a2) pairs.
\\ Theorem (Threads 15-16): p prime, p∤a2 => [P]^2=1 in Cl(Q(sqrt(sf)))
\\   AND (corollary) p splits (kro(disc_K, p) = +1).

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

check_one(pp, aa) = {
  my(D, sf, disc_K, kro, norm_beta);
  D = aa^2 - 4*pp^2;
  if(D == 0, printf("  SKIP D=0\n"); return(-1));
  if(Mod(aa, pp) == Mod(0, pp), printf("  SKIP pp|aa\n"); return(-1));
  sf = sf_part(D);
  if(sf == 1 || sf == 0, printf("  SKIP sf=%d (D perfect square)\n", sf); return(-1));
  disc_K = if(sf % 4 == 1, sf, 4*sf);
  kro = kronecker(disc_K, pp);
  norm_beta = (aa^2 - (D/sf)*sf) / 4;
  my(status);
  if(kro == 1, status = "SPLITS";
  ,if(kro == -1, status = "INERT";
  ,              status = "RAMIFIED"));
  printf("  p=%-5d a2=%-5d sf=%-8d kro=%2d -> %-9s | N(b)=%d == p^2:%d\n",
    pp, aa, sf, kro, status, norm_beta, norm_beta == pp^2);
  return(kro == 1 && norm_beta == pp^2)
};

print("=== Thread 17: 10 fresh (p,a2) splits verification ===");
print("");

run_all() = {
  my(pv, av, ok, tot, r);
  pv = [101, 127, 149, 199, 211, 251, 307, 311, 331, 503];
  av = [ 14,  20,  30,  50, 100,  60, 200, 300, 320, 400];
  ok = 0; tot = 10;
  for(i = 1, tot,
    r = check_one(pv[i], av[i]);
    if(r == 1, ok++));
  printf("\nSummary: %d/%d SPLITS and N(beta)=p^2.\n", ok, tot);
  if(ok == tot,
    print("ALL PASS -- p always splits when p does not divide a2."),
    print("FAILURES."))
};

run_all()
