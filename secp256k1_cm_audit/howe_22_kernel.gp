\\ howe_22_kernel.gp
\\ Identify the Howe-glued class among the 7 Z/3Z canonical classes over F_43.
\\
\\ Strategy:
\\ For each class (a,b), the sextic y^2=(x^3-r1)(x^3-r2) has two natural
\\ elliptic sub-factors E_{r1}: y^2=x^3-r1 and E_{r2}: y^2=x^3-r2.
\\ The Howe-glued class for (E1=y^2=x^3+7, E2=y^2=x^3+13) is the one where
\\ {trace(E_{r1}), trace(E_{r2})} = {+13, -13}.
\\
\\ Also: enumerate all 15 quadratic factorizations of the sextic f = product of
\\ 6 Weierstrass points and check each for the degenerate Richelot to E1xE2.
\\
\\ Run: gp -q howe_22_kernel.gp

default(parisize, 256000000);
default(timer, 0);

p = 43;
ca = [1,2,3,4,6,8,16];
cb = [2,8,42,32,39,42,39];
t_target = 13;

\\ Cube root of r in F_43 if it exists, else -1
cube_root_fp(r) = {
  my(t);
  if(r == 0, return(0));
  t = lift(Mod(r,p)^((p-1)/3));
  if(t != 1, return(-1));
  lift(Mod(r,p)^((2*(p-1)/3+1)/3))
};

\\ All cube roots of r in F_43 (0, 1, or 3 roots)
all_cube_roots_fp(r) = {
  my(alpha0, z3, v);
  alpha0 = cube_root_fp(r);
  if(alpha0 < 0, return([]));
  z3 = lift(polrootsmod(Mod(1,p)*x^2+x+1, p)[1]);
  v = List();
  listput(v, alpha0);
  listput(v, lift(Mod(alpha0*z3, p)));
  listput(v, lift(Mod(alpha0*z3^2, p)));
  Vec(v)
};

print("================================================================");
print("Step 1: Sub-factors E_{r1} and E_{r2} for each Z/3Z class");
print("================================================================");
print("");
print("For class (a,b): roots r1,r2 of T^2+aT+b=0 over F_43.");
print("E_{ri}: y^2=x^3-ri  (= y^2=x^3+(-ri mod 43))");
print("Check if traces match {+13,-13} (Howe pair).");
print("");

inv2 = lift(Mod(2,p)^(-1));

for(i=1,7,
  my(aa,bb,disc,sq,r1,r2,mr1,mr2,Er1,Er2,t1,t2,flag);
  aa=ca[i]; bb=cb[i];
  disc = lift(Mod(aa^2-4*bb,p));
  sq = lift(Mod(disc,p)^((p+1)/4));
  r1 = lift(Mod((-aa+sq)*inv2,p));
  r2 = lift(Mod((-aa-sq)*inv2,p));
  mr1 = lift(Mod(-r1,p));
  mr2 = lift(Mod(-r2,p));
  Er1 = ellinit([0,mr1],p);
  Er2 = ellinit([0,mr2],p);
  t1 = p+1-ellcard(Er1);
  t2 = p+1-ellcard(Er2);
  flag = if((abs(t1)==t_target && abs(t2)==t_target && t1*t2<0), " <-- HOWE MATCH!", "");
  print("Class ",i," (a=",aa," b=",bb,"): r1=",r1," r2=",r2,
    "  E_{r1}=y^2=x^3+",mr1," t=",t1,
    "  E_{r2}=y^2=x^3+",mr2," t=",t2,flag)
);

print("");
print("================================================================");
print("Step 2: Exhaustive 15-factorization Richelot over F_{43^3}");
print("For each class, find all 15 pairings of Weierstrass points");
print("and check which (if any) degenerate Richelot gives E1 x E2.");
print("================================================================");
print("");

\\ F_{43^3} arithmetic: represent as [a,b,c] = a + b*w + c*w^2, w^3=36
fadd3(u,v) = [(u[1]+v[1])%p,(u[2]+v[2])%p,(u[3]+v[3])%p];
fneg3(u) = [(-u[1]+p)%p,(-u[2]+p)%p,(-u[3]+p)%p];
fscl3(c,u) = [(c*u[1])%p,(c*u[2])%p,(c*u[3])%p];
fmul3(u,v) = {
  [(u[1]*v[1]+36*(u[2]*v[3]+u[3]*v[2]))%p,
   (u[1]*v[2]+u[2]*v[1]+36*u[3]*v[3])%p,
   (u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%p]
};

\\ Compute cube root of r in F_{43^3} using the PARI extension field
\\ Returns [alpha0, alpha1, alpha2] as F_{43^3} elements, plus zeta3
compute_cube_roots_fp3(r) = {
  my(ww, z3, alpha0, alpha1, alpha2);
  \\ F_{43^3} = F_43[w]/(w^3-r)
  \\ alpha0 = w (in this representation, w^3=r so [w]=[0,1,0] in F_43[t]/(t^3-r))
  \\ But our arithmetic uses w^3=36, so we need to work in F_43[t]/(t^3-r) separately.
  \\ Instead: use cube_root_fp3 logic from prior scripts.
  \\ tier-0: r is a cube in F_43 -> alpha0 in F_43
  \\ tier-1: r/36 is a cube -> alpha0 = c*w (w^3=36)
  \\ tier-2: r/6 is a cube -> alpha0 = c*w^2 (w^3=36, so (c*w^2)^3=c^3*36^2=r -> c^3=r/36^2=r/6 mod p)
  my(t0, t1, t2, c);
  t0 = lift(Mod(r,p)^((p-1)/3));
  if(t0 == 1,
    c = lift(Mod(r,p)^((2*(p-1)/3+1)/3));
    alpha0 = [c,0,0];
    z3 = lift(polrootsmod(Mod(1,p)*x^2+x+1,p)[1]);
    alpha1 = [lift(Mod(c*z3,p)),0,0];
    alpha2 = [lift(Mod(c*z3^2,p)),0,0];
    return([alpha0,alpha1,alpha2])
  );
  my(r36 = lift(Mod(r * Mod(36,p)^(-1), p)));
  t1 = lift(Mod(r36,p)^((p-1)/3));
  if(t1 == 1,
    c = lift(Mod(r36,p)^((2*(p-1)/3+1)/3));
    z3 = lift(polrootsmod(Mod(1,p)*x^2+x+1,p)[1]);
    alpha0 = [0,c,0];
    alpha1 = [0,lift(Mod(c*z3,p)),0];
    alpha2 = [0,lift(Mod(c*z3^2,p)),0];
    return([alpha0,alpha1,alpha2])
  );
  my(r6 = lift(Mod(r * Mod(6,p)^(-1), p)));
  t2 = lift(Mod(r6,p)^((p-1)/3));
  if(t2 == 1,
    c = lift(Mod(r6,p)^((2*(p-1)/3+1)/3));
    z3 = lift(polrootsmod(Mod(1,p)*x^2+x+1,p)[1]);
    alpha0 = [0,0,c];
    alpha1 = [0,0,lift(Mod(c*z3,p))];
    alpha2 = [0,0,lift(Mod(c*z3^2,p))];
    return([alpha0,alpha1,alpha2])
  );
  error("no cube root for r=",r)
};

\\ For a pairing of 6 points into 3 pairs, compute the "degenerate Richelot" target.
\\ Input: three pairs {[P1,P2],[P3,P4],[P5,P6]} where each Pi is in F_{43^3}.
\\ G_i = (x-Pi)(x-Pj) has:
\\   G_i_x = -(Pi+Pj)  (linear coefficient)
\\   G_i_c = Pi*Pj     (constant term)
\\ If G1_c = G2_c = G3_c (Delta=0), compute the degenerate target.
\\ Returns [trace1, trace2] of the two elliptic factor curves, or [-999,-999] if not degenerate.
degenerate_richelot_traces(W) = {
  my(G1s, G1c, G2s, G2c, G3s, G3c, qval, diff1, diff2, diff3);
  my(s1, s2, s3, e1, e2, e3, E1c, E2c, tE1, tE2);
  G1s = fadd3(W[1], W[2]);  G1c = fmul3(W[1], W[2]);
  G2s = fadd3(W[3], W[4]);  G2c = fmul3(W[3], W[4]);
  G3s = fadd3(W[5], W[6]);  G3c = fmul3(W[5], W[6]);
  diff1 = fadd3(G1c, fneg3(G2c));
  diff2 = fadd3(G2c, fneg3(G3c));
  diff3 = fadd3(G1c, fneg3(G3c));
  \\ Check if all constant terms equal (degenerate case)
  if(diff1[1]!=0||diff1[2]!=0||diff1[3]!=0||diff2[1]!=0||diff2[2]!=0||diff2[3]!=0,
    return([-999,-999]));
  \\ All G_ic equal. This is a degenerate Richelot.
  \\ The constant term q = G1_c must be in F_43 (not F_{43^3}) for E1,E2 to be over F_43.
  qval = G1c;
  if(qval[2]!=0||qval[3]!=0, return([-998,-998]));
  \\ q is in F_43. Now find the two elliptic curves.
  \\ The three pairs {Wi, Wj} have constant terms all equal to q.
  \\ The sextic factors as (x^3 - r_alpha)(x^3 - r_beta) where r_alpha = product of one group's roots^3.
  \\ Alternatively: compute the traces directly from the group structure.
  \\ The 6 Weierstrass points split into 2 natural elliptic groups.
  \\ But with a general pairing, the groups may interleave.
  \\ Key: the two elliptic curves have 2-torsion given by the COLUMN vs ROW grouping.
  \\ For degenerate Richelot with equal constant terms:
  \\ The two factors correspond to the two HALVES of the factorization.
  \\ Specifically: {W[1],W[3],W[5]} form one elliptic curve and {W[2],W[4],W[6]} form another.
  \\   E1': y^2 = (x-W[1])(x-W[3])(x-W[5]) but these are 2-torsion points, so:
  \\   x^3 - r_alpha = (x-W[1])(x-W[3])(x-W[5])
  \\   r_alpha = W[1]*W[3]*W[5]  (product of the 3 x-coords)
  \\ This product must be in F_43 for E_alpha to be defined over F_43.
  my(prod135, prod246, r_a, r_b);
  prod135 = fmul3(fmul3(W[1],W[3]),W[5]);
  prod246 = fmul3(fmul3(W[2],W[4]),W[6]);
  if(prod135[2]!=0||prod135[3]!=0, return([-997,-997]));
  if(prod246[2]!=0||prod246[3]!=0, return([-997,-997]));
  r_a = prod135[1];  \\ r_alpha: the "a" group product
  r_b = prod246[1];  \\ r_beta: the "b" group product
  \\ E_alpha: y^2 = x^3 - r_a  =>  y^2 = x^3 + (-r_a mod p)
  \\ E_beta:  y^2 = x^3 - r_b  =>  y^2 = x^3 + (-r_b mod p)
  E1c = lift(Mod(-r_a, p));
  E2c = lift(Mod(-r_b, p));
  tE1 = p+1-ellcard(ellinit([0,E1c],p));
  tE2 = p+1-ellcard(ellinit([0,E2c],p));
  [tE1, tE2]
};

\\ All 15 ways to partition 6 elements into 3 pairs
\\ Elements indexed 1..6. Returns list of [pair1, pair2, pair3] as index-pairs.
all_15_pairings() = {
  my(pairings, idx, p1, p2, p3);
  pairings = List();
  for(a=2,6,
    p1 = [1,a];
    for(b=1,5,
      if(b==1 || b==a, next);
      for(c=b+1,6,
        if(c==a, next);
        p2 = [b,c];
        p3 = [];
        for(d=1,6,
          if(d!=1 && d!=a && d!=b && d!=c, p3 = concat(p3,[d]))
        );
        if(#p3==2 && p3[1]<p3[2],
          listput(pairings, [p1, p2, p3])
        )
      )
    )
  );
  Vec(pairings)
};

pairings15 = all_15_pairings();
print("Total pairings enumerated: ", #pairings15, " (expected 15)");
print("");

for(i=1,7,
  my(aa,bb,disc,sq,r1,r2,alpha_roots,beta_roots,W6,found_howe);
  aa=ca[i]; bb=cb[i];
  disc = lift(Mod(aa^2-4*bb,p));
  sq = lift(Mod(disc,p)^((p+1)/4));
  r1 = lift(Mod((-aa+sq)*inv2,p));
  r2 = lift(Mod((-aa-sq)*inv2,p));
  \\ alpha_i = cube roots of r1, beta_i = cube roots of r2
  alpha_roots = compute_cube_roots_fp3(r1);
  beta_roots  = compute_cube_roots_fp3(r2);
  W6 = [alpha_roots[1],alpha_roots[2],alpha_roots[3],
        beta_roots[1], beta_roots[2], beta_roots[3]];
  print("Class ",i," (a=",aa," b=",bb,")  r1=",r1,"  r2=",r2);
  found_howe = 0;
  for(pi=1,15,
    my(pr,idx1,idx2,idx3,idx4,idx5,idx6,W,res);
    pr = pairings15[pi];
    idx1=pr[1][1]; idx2=pr[1][2];
    idx3=pr[2][1]; idx4=pr[2][2];
    idx5=pr[3][1]; idx6=pr[3][2];
    W = [W6[idx1],W6[idx2],W6[idx3],W6[idx4],W6[idx5],W6[idx6]];
    res = degenerate_richelot_traces(W);
    if(res[1]==-999, next);  \\ not degenerate
    if(res[1]<-900,
      print("  pairing ",pi," [",idx1,",",idx2,"|",idx3,",",idx4,"|",idx5,",",idx6,"]: degenerate but E not over F_43 (code=",res[1],")");
      next
    );
    print("  pairing ",pi," [",idx1,",",idx2,"|",idx3,",",idx4,"|",idx5,",",idx6,"]: ",
      "E_a trace=",res[1],"  E_b trace=",res[2],
      if(abs(res[1])==t_target && abs(res[2])==t_target && res[1]*res[2]<0,
        " <== HOWE PAIR!", ""));
    if(abs(res[1])==t_target && abs(res[2])==t_target && res[1]*res[2]<0, found_howe=1)
  );
  if(found_howe==0, print("  (no degenerate Richelot gives Howe pair {+-13} for this class)"));
  print("")
);

print("================================================================");
print("Done.");
