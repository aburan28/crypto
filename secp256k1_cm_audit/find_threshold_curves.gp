\\
\\ Find y^2 = x^3 + 2 curves over F_p with lambda/n in [0.07, 0.34]
\\ Used for Thread 20: GLV-HNP Phase 2 lambda/n threshold bisection.
\\
\\ Output: p, n, lambda, lambda/n  (tab-separated)
\\
\\ Run: gp -q find_threshold_curves.gp
\\

print("p\tn\tlambda\tlambda/n");

{
  my(E, n, disc, sq, inv2, lam1, lam2, ratio, lam);
  forprime(p=500, 60000,
    if(p%3!=1, next);
    E = ellinit([0,0,0,0,2], p);
    n = ellcard(E);
    if(!isprime(n), next);
    \\ Roots of x^2+x+1 mod n: need -3 to be QR mod n
    disc = Mod(-3, n);
    if(!issquare(disc), next);
    sq = sqrt(disc);
    inv2 = Mod(1,n)/2;
    lam1 = lift((-1+sq)*inv2);
    lam2 = lift((-1-sq)*inv2);
    \\ Pick root in (0, n/2) to use as the "small" eigenvalue
    if(lam1 > 0 && lam1 < n/2,
      lam = lam1; ratio = lam1*1.0/n,
      if(lam2 > 0 && lam2 < n/2,
        lam = lam2; ratio = lam2*1.0/n,
        next));
    if(ratio > 0.07 && ratio < 0.35,
      printf("%d\t%d\t%d\t%.5f\n", p, n, lam, ratio));
  );
}
