\\ howe_richelot_v3.gp
\\ Compute ALL THREE Galois-equivariant Richelot images of
\\ y^2=(x^3+7)(x^3+13) over F_43 (the three sigma_k matchings).
\\ Also check: these images lie in WRONG isogeny class,
\\ confirming the Richelot of the naive cover is not the Howe-glued curve.
\\ THEN: do the INVERSE problem — which of the 7 canonical Z/3Z classes
\\ has ONE of its Richelot images equal to E1 x E2 (the split surface)?
\\ That class is the Howe-glued curve.
\\ Run: gp -q howe_richelot_v3.gp

default(parisize, 64000000); default(timer, 0);
p = 43; z3 = 6; z3sq = 36;

\\ F_{43^3} = F_43[a]/(a^3 - 36). Elements: [c0,c1,c2].
fadd(u,v) = [(u[1]+v[1])%p,(u[2]+v[2])%p,(u[3]+v[3])%p];
fneg(u) = [(-u[1])%p,(-u[2])%p,(-u[3])%p];
fscl(c,u) = [(c*u[1])%p,(c*u[2])%p,(c*u[3])%p];
fcon(c) = [c%p,0,0];
fmul(u,v) = [(u[1]*v[1]+36*(u[2]*v[3]+u[3]*v[2]))%p,(u[1]*v[2]+u[2]*v[1]+36*u[3]*v[3])%p,(u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%p];

a_ = [0,1,0]; a2_ = [0,0,1];

\\ Richelot with Z/3Z-equivariant pairing: s=S*alpha, q=Q*alpha^2
\\ Returns [a_monic, b_monic] of the image y^2 = x^6 + a*x^3 + b,
\\ or [-1,-1] if result is not Z/3Z form.
richelot_zt3(S, Q) = {
  my(s,q,G1c,G1x,G1x2,G2c,G2x,G2x2,G3c,G3x,G3x2);
  my(G1pc,G1px,G2pc,G2px,G3pc,G3px,D0,D0sc,D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r);
  my(lc,lcinv,aa,bb);
  s = fscl(S, a_); q = fscl(Q, a2_);
  G1c=q; G1x=fneg(s); G1x2=fcon(1);
  G2c=fscl(z3sq,q); G2x=fneg(fscl(z3,s)); G2x2=fcon(1);
  G3c=fscl(z3,q); G3x=fneg(fscl(z3sq,s)); G3x2=fcon(1);
  G1px=fcon(2); G1pc=G1x;
  G2px=fcon(2); G2pc=G2x;
  G3px=fcon(2); G3pc=G3x;
  D0=fadd(fadd(fmul(G1x2,fadd(fmul(G2x,G3c),fneg(fmul(G3x,G2c)))),fneg(fmul(G1x,fadd(fmul(G2x2,G3c),fneg(fmul(G3x2,G2c)))))),fmul(G1c,fadd(fmul(G2x2,G3x),fneg(fmul(G3x2,G2x)))));
  if (D0[2]!=0||D0[3]!=0, print("Delta not scalar for S=",S," Q=",Q); return([-1,-1]));
  D0sc = D0[1];
  if (D0sc==0, print("Delta=0 for S=",S," Q=",Q); return([-1,-1]));
  D0inv = fcon(lift(Mod(D0sc,p)^(-1)));
  H1n2=fadd(fadd(fmul(G2px,G3x),fmul(G2pc,G3x2)),fneg(fadd(fmul(G2x2,G3pc),fmul(G2x,G3px))));
  H1n1=fadd(fadd(fmul(G2px,G3c),fmul(G2pc,G3x)),fneg(fadd(fmul(G2x,G3pc),fmul(G2c,G3px))));
  H1n0=fadd(fmul(G2pc,G3c),fneg(fmul(G2c,G3pc)));
  H2n2=fadd(fadd(fmul(G3px,G1x),fmul(G3pc,G1x2)),fneg(fadd(fmul(G3x2,G1pc),fmul(G3x,G1px))));
  H2n1=fadd(fadd(fmul(G3px,G1c),fmul(G3pc,G1x)),fneg(fadd(fmul(G3x,G1pc),fmul(G3c,G1px))));
  H2n0=fadd(fmul(G3pc,G1c),fneg(fmul(G3c,G1pc)));
  H3n2=fadd(fadd(fmul(G1px,G2x),fmul(G1pc,G2x2)),fneg(fadd(fmul(G1x2,G2pc),fmul(G1x,G2px))));
  H3n1=fadd(fadd(fmul(G1px,G2c),fmul(G1pc,G2x)),fneg(fadd(fmul(G1x,G2pc),fmul(G1c,G2px))));
  H3n0=fadd(fmul(G1pc,G2c),fneg(fmul(G1c,G2pc)));
  H1x2=fmul(H1n2,D0inv); H1x1=fmul(H1n1,D0inv); H1x0=fmul(H1n0,D0inv);
  H2x2=fmul(H2n2,D0inv); H2x1=fmul(H2n1,D0inv); H2x0=fmul(H2n0,D0inv);
  H3x2=fmul(H3n2,D0inv); H3x1=fmul(H3n1,D0inv); H3x0=fmul(H3n0,D0inv);
  P0=fmul(H1x0,H2x0); P1=fadd(fmul(H1x0,H2x1),fmul(H1x1,H2x0));
  P2=fadd(fadd(fmul(H1x0,H2x2),fmul(H1x1,H2x1)),fmul(H1x2,H2x0));
  P3=fadd(fmul(H1x1,H2x2),fmul(H1x2,H2x1)); P4=fmul(H1x2,H2x2);
  Q0r=fmul(P0,H3x0); Q1r=fadd(fmul(P0,H3x1),fmul(P1,H3x0));
  Q2r=fadd(fadd(fmul(P0,H3x2),fmul(P1,H3x1)),fmul(P2,H3x0));
  Q3r=fadd(fadd(fmul(P1,H3x2),fmul(P2,H3x1)),fmul(P3,H3x0));
  Q4r=fadd(fadd(fmul(P2,H3x2),fmul(P3,H3x1)),fmul(P4,H3x0));
  Q5r=fadd(fmul(P3,H3x2),fmul(P4,H3x1)); Q6r=fmul(P4,H3x2);
  if (Q0r[2]!=0||Q0r[3]!=0||Q1r[2]!=0||Q1r[3]!=0||Q2r[2]!=0||Q2r[3]!=0||Q3r[2]!=0||Q3r[3]!=0||Q4r[2]!=0||Q4r[3]!=0||Q5r[2]!=0||Q5r[3]!=0||Q6r[2]!=0||Q6r[3]!=0, print("Non-F43 coeffs for S=",S," Q=",Q); return([-1,-1]));
  lc=Q6r[1];
  if (lc==0, return([-1,-1]));
  lcinv=lift(Mod(lc,p)^(-1));
  aa=(Q3r[1]*lcinv)%p; bb=(Q0r[1]*lcinv)%p;
  [aa,bb]
};

\\ The three sigma_k matchings for naive cover y^2=(x^3+7)(x^3+13):
\\ alpha^3=36, beta=2*alpha, z3=6
\\ sigma_0: s=3*alpha, q=2*alpha^2
\\ sigma_1: s=(1+2*z3)*alpha=13*alpha, q=2*z3*alpha^2=12*alpha^2
\\ sigma_2: s=(1+2*z3^2)*alpha=(1+72)*alpha=30*alpha, q=2*z3^2*alpha^2=2*36*alpha^2=72*alpha^2=29*alpha^2
print("=== Naive cover: 3 Richelot images ===");
r0 = richelot_zt3(3, 2); print("sigma_0: a=",r0[1],"  b=",r0[2]);
r1 = richelot_zt3(13, 12); print("sigma_1: a=",r1[1],"  b=",r1[2]);
r2 = richelot_zt3(30, 29); print("sigma_2: a=",r2[1],"  b=",r2[2]);
print("");

\\ Verify char poly of each Richelot image
for (i=1,3,
  res = [r0,r1,r2][i];
  aa=res[1]; bb=res[2];
  if (aa>=0,
    h = Mod(aa,p)*x^3*x^3 + Mod(aa,p)*x^3 + Mod(bb,p);
    h = (Mod(1,p)*x^6) + Mod(aa,p)*x^3 + Mod(bb,p);
    cp = hyperellcharpoly(h);
    nj = subst(cp,variable(cp),1);
    print("sigma_",i-1," char poly: ",cp,"  #Jac=",nj)
  )
);
print("");

\\ === Main question: for each of the 7 Z/3Z classes, compute its 3 Richelot images.
\\ The Howe-glued class will have one image that "looks like" E1 x E2.
\\ Over F_43, E1 x E2 is a DEGENERATE ppav: the naive cover y^2=(x^3+7)(x^3+13)
\\ has its Jacobian in a DIFFERENT isogeny class from E1 x E2.
\\ So "the Richelot image IS E1 x E2" means the image H1H2H3 factors as (x^3-r1)(x^3-r2)
\\ where y^2=x^3-r1 has the same char poly as E1 and y^2=x^3-r2 as E2.
\\ Equivalently: the image sextic should be y^2=(x^3+7)(x^3+13) = x^6+20x^3+5 in Z/3Z form... BUT
\\ the "degenerate" case gives Δ=0, not a normal Richelot.
\\
\\ ALTERNATIVE APPROACH: instead of Richelot, check which class has
\\ the Frobenius char poly equal to target T^4-83T^2+1849 AND
\\ its Igusa invariants lie on the Howe-glued locus.
\\
\\ For now: compute Richelot images of ALL 7 classes and check char polys.
\\
print("=== Richelot images of the 7 canonical Z/3Z classes ===");
print("(Correct Howe-glued class should produce image with char poly T^4-83T^2+1849)");
print("");

ca=[1,2,3,4,6,8,16]; cb=[2,8,42,32,39,42,39];

for (i=1,7,
  aa=ca[i]; bb=cb[i];
  print("Class ",i,": a=",aa,"  b=",bb);
  \\ Need to find alpha_i, beta_i as cube roots of -(roots of T^2+aT+b=0) in F_{43^3}
  \\ T^2 + aa*T + bb = 0. Roots: r1,r2 = (-aa ± sqrt(aa^2-4bb)) / 2
  disc_ab = (aa^2 - 4*bb) % p;
  \\ Check if disc_ab is a square mod p
  if (disc_ab == 0,
    r1r2 = [(-aa * lift(Mod(2,p)^(-1))) % p];
    print("  discriminant=0: r1=r2=",r1r2[1])
  ,
    sq = lift(Mod(disc_ab, p)^((p-1)/2));
    if (sq == 1,
      \\ disc is a square, roots in F_43
      sqrtd = lift(Mod(disc_ab,p)^((p+1)/4));
      two_inv = lift(Mod(2,p)^(-1));
      r1 = ((-aa + sqrtd)*two_inv) % p;
      r2 = ((-aa - sqrtd + p)*two_inv) % p;
      print("  disc=",disc_ab,"  sqrtd=",sqrtd,"  r1=",r1,"  r2=",r2);
      \\ E1: y^2=x^3-r1=x^3+(-r1 mod p), trace check
      nr1 = (-r1) % p; nr2 = (-r2) % p;
      E_r1 = ellinit([0,nr1],p); tr1 = p+1-ellcard(E_r1);
      E_r2 = ellinit([0,nr2],p); tr2 = p+1-ellcard(E_r2);
      print("  E: y^2=x^3+",nr1," trace=",tr1,"  E: y^2=x^3+",nr2," trace=",tr2);
      \\ Check if one of these matches E1 (trace=13) and E2 (trace=-13)
      if ((tr1==13&&tr2==-13)||(tr1==-13&&tr2==13), print("  *** HOWE CANDIDATE (product of E1,E2) ***"))
    ,
      print("  disc=",disc_ab," not a square => r1,r2 in F_{43^2}, not handled here")
    )
  );
  \\ Also compute 3 Richelot images of this class
  \\ To do this, we need the cube roots of r1 and r2 in F_{43^3}
  \\ and the corresponding (S,Q) parameters.
  \\ This requires working over F_{43^3} — skipping for now, noted as BLOCKED.
  print("")
);

print("=== Summary ===");
print("A class with disc(T^2+aT+b) a square and whose roots give traces ±13");
print("would be the Howe-glued curve via the 'direct product' structure.");
print("Other classes require F_{43^2} roots — need quadratic extension work.");
