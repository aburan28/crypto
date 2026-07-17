\\ thread16_chlrs_rosenhain_fix.gp — exhaustive Möbius search (PARI-correct)
\\ Run: gp -q secp256k1_cm_audit/thread16_chlrs_rosenhain_fix.gp
default(parisize, 64000000);
default(timer, 0);
print("=== Thread 16: CHLRS Rosenhain normalization fix ===");

\\ Parameters
p = 1009; b0 = Mod(11,p); t_E1 = 43;
a2_target = 2*p - t_E1^2;
printf("p=%d b=%d target_a2=%d\n\n", p, lift(b0), a2_target);

\\ omega: cube root of unity in F_p
om = polrootsmod(x^2+x+1,p)[1];

\\ d: first non-square mod p
d_ns = 0;
fordiv(1,p, for(dd=2,p-1, if(kronecker(dd,p)==-1, d_ns=dd; break); break));
\\ simpler:
d_ns=0; for(dd=2,p, if(kronecker(dd,p)==-1, d_ns=dd; break));
printf("omega=%Ps d_ns=%d\n", om, d_ns);

\\ F_{p^3} via ffgen(p^3, 'a)
F3a = ffgen(p^3, 'a);
printf("F_%d^3 constructed, generator a\n", p);

\\ Roots of x^3 + b0 in F_{p^3}:
\\ polrootsff(polynomial-over-Fp, p, generator-of-Fpn)
roots1 = polrootsff(x^3 + lift(b0), p, F3a);
printf("roots of x^3+%d in F_{p^3}: %d roots\n", lift(b0), #roots1);
if(#roots1 != 3, print("ERROR: expected 3 roots in F_{p^3}"));

\\ E2 roots: multiply each E1 root by d_ns
roots2 = vector(3, i, Mod(d_ns,p)*roots1[i]);
printf("E2 (x^3+%d) roots: %d roots\n\n", lift(Mod(d_ns^3,p)*b0), #roots2);

\\ Concatenate: pts[1..3]=E1 roots, pts[4..6]=E2 roots
pts = vector(6, i, if(i<=3, roots1[i], roots2[i-3]));

\\ Is elt in F_p? <=> elt^p == elt
in_Fp(elt) = (elt^p == elt);

\\ Rosenhain Frob a2: C: y^2=x(x-1)(x-l1)(x-l2)(x-l3) over F_p
\\ l1,l2,l3 are integers (already reduced mod p)
rosen_a2(l1,l2,l3) = {
  my(fld, h5, fp);
  fld = ffgen(p,'u);
  h5 = fld^0*x*(x-Mod(1,p))*(x-Mod(l1,p))*(x-Mod(l2,p))*(x-Mod(l3,p));
  polcoeff(hyperellcharpoly(h5), 2)
};

\\ Möbius image: T s.t. P1->0, P2->1, P3->inf
\\ T(Q) = (Q-P1)(P2-P3) / ((Q-P3)(P2-P1))
\\ returns -1 on degenerate (den=0)
mob(Q, P1, P2, P3) = {
  my(num, den);
  num = (Q-P1)*(P2-P3);
  den = (Q-P3)*(P2-P1);
  if(den==0, return(-1));
  num/den
};

\\ Main search
print("--- Exhaustive Möbius search (120 ordered triples, 6 pts) ---");
n_valid=0; n_match=0;

run_search(show_all) = {
  my(P1,P2,P3,lv,lam,l1,l2,l3,a2v,ok,ri);
  my(a2_dist_set);
  a2_dist_set = Set([]);
  for(i1=1,6, for(i2=1,6,
    if(i2==i1, next());
    for(i3=1,6,
      if(i3==i1||i3==i2, next());
      P1=pts[i1]; P2=pts[i2]; P3=pts[i3];
      \\ remaining 3 indices
      ri=vector(3); ok=1;
      {
        my(cnt=0);
        for(jj=1,6,
          if(jj==i1||jj==i2||jj==i3, next());
          cnt++; ri[cnt]=jj
        )
      };
      \\ compute lambdas
      lv=vector(3); ok=1;
      for(jj=1,3,
        lam=mob(pts[ri[jj]],P1,P2,P3);
        if(lam==-1, ok=0; break);
        if(!in_Fp(lam), ok=0; break);
        lv[jj]=lam
      );
      if(!ok, next());
      l1=lift(lv[1]); l2=lift(lv[2]); l3=lift(lv[3]);
      if(l1==0||l1==1||l2==0||l2==1||l3==0||l3==1, next());
      if(l1==l2||l1==l3||l2==l3, next());
      n_valid++;
      a2v=rosen_a2(l1,l2,l3);
      a2_dist_set=setunion(a2_dist_set,Set([a2v]));
      if(a2v==a2_target,
        n_match++;
        printf("  MATCH: (%d->0, %d->1, %d->inf) lam=(%d,%d,%d) a2=%d\n",
          i1,i2,i3,l1,l2,l3,a2v);
      , if(show_all,
        printf("  triple=(%d,%d,%d) lam=(%d,%d,%d) a2=%d\n",
          i1,i2,i3,l1,l2,l3,a2v)
      ))
    )
  ));
  printf("Valid triples: %d  Matches: %d\n", n_valid, n_match);
  printf("Distinct a2 values seen: %Ps\n", a2_dist_set);
  if(setsearch(a2_dist_set,[a2_target]),
    print("target a2 IS present"),
    print("target a2 ABSENT from all valid Rosenhain models"))
};

run_search(0);
print("");
print("=== DONE: thread16_chlrs_rosenhain_fix.gp ===");
