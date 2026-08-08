\\ howe_22_kernel_v2.gp
\\
\\ Thread 3 (CHLRS Igusa forward map), 2026-08-08 autolab session.
\\
\\ Fixes the two known PARI 2.15.4 bugs in howe_22_kernel.gp (2026-06-09 log):
\\ (a) multiline for-loop / array-literal bodies outside a single {...} block
\\ are parsed line-by-line by gp script mode, silently corrupting loop vars
\\ into free polynomial indeterminates; (b) all_15_pairings() returned 30.
\\ Fix: every multi-line statement lives inside exactly one non-nested {...}
\\ block (gp 2.15.4 does not support nested braces), and the 15 perfect
\\ matchings of {1..6} are hardcoded as one single-line literal.
\\
\\ Goal: identify which of the 7 canonical Z/3Z classes over F_43
\\ (char poly T^4-83T^2+1849, tabulated 2026-06-03) is the actual Howe-glued
\\ curve for (E1: y^2=x^3+7, E2: y^2=x^3+13).
\\
\\ Run: gp -q howe_22_kernel_v2.gp

default(parisize, 256000000);

p = 43;
ca = [1,2,3,4,6,8,16];
cb = [2,8,42,32,39,42,39];
t_target = 13;
inv2 = lift(Mod(2,p)^(-1));

M15 = [[[1,2],[3,4],[5,6]],[[1,2],[3,5],[4,6]],[[1,2],[3,6],[4,5]],[[1,3],[2,4],[5,6]],[[1,3],[2,5],[4,6]],[[1,3],[2,6],[4,5]],[[1,4],[2,3],[5,6]],[[1,4],[2,5],[3,6]],[[1,4],[2,6],[3,5]],[[1,5],[2,3],[4,6]],[[1,5],[2,4],[3,6]],[[1,5],[2,6],[3,4]],[[1,6],[2,3],[4,5]],[[1,6],[2,4],[3,5]],[[1,6],[2,5],[3,4]]];
print("Hardcoded matchings: ", #M15, " (expect 15)");

fadd3(u,v) = [(u[1]+v[1])%p,(u[2]+v[2])%p,(u[3]+v[3])%p];
fneg3(u)   = [(-u[1])%p,(-u[2])%p,(-u[3])%p];
fmul3(u,v) = [(u[1]*v[1]+36*(u[2]*v[3]+u[3]*v[2]))%p, (u[1]*v[2]+u[2]*v[1]+36*u[3]*v[3])%p, (u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%p];

cube_roots_fp3(r) = {
  my(z3, t0, t1, t2, c, r36, r6);
  z3 = lift(polrootsmod(Mod(1,p)*x^2+x+1,p)[1]);
  t0 = lift(Mod(r,p)^((p-1)/3));
  if(t0 == 1, c = lift(Mod(r,p)^((2*(p-1)/3+1)/3)); return([[c,0,0],[lift(Mod(c*z3,p)),0,0],[lift(Mod(c*z3^2,p)),0,0]]));
  r36 = lift(Mod(r,p)/Mod(36,p));
  t1 = lift(Mod(r36,p)^((p-1)/3));
  if(t1 == 1, c = lift(Mod(r36,p)^((2*(p-1)/3+1)/3)); return([[0,c,0],[0,lift(Mod(c*z3,p)),0],[0,lift(Mod(c*z3^2,p)),0]]));
  r6 = lift(Mod(r,p)/Mod(6,p));
  t2 = lift(Mod(r6,p)^((p-1)/3));
  if(t2 == 1, c = lift(Mod(r6,p)^((2*(p-1)/3+1)/3)); return([[0,0,c],[0,0,lift(Mod(c*z3,p))],[0,0,lift(Mod(c*z3^2,p))]]));
  error("cube_roots_fp3: no cube root for r=", r)
};

richelot_deg_traces(W) = {
  my(G1c,G2c,G3c,d12,d23,prodA,prodB,rA,rB,tA,tB);
  G1c = fmul3(W[1],W[2]);
  G2c = fmul3(W[3],W[4]);
  G3c = fmul3(W[5],W[6]);
  d12 = fadd3(G1c, fneg3(G2c));
  d23 = fadd3(G2c, fneg3(G3c));
  if(d12!=[0,0,0] || d23!=[0,0,0], return([-999,-999,-999]));
  prodA = fmul3(fmul3(W[1],W[3]),W[5]);
  prodB = fmul3(fmul3(W[2],W[4]),W[6]);
  if(prodA[2]!=0 || prodA[3]!=0, return([-998,-998,-998]));
  if(prodB[2]!=0 || prodB[3]!=0, return([-998,-998,-998]));
  rA = prodA[1]; rB = prodB[1];
  tA = p+1-ellcard(ellinit([0, lift(Mod(-rA,p))], p));
  tB = p+1-ellcard(ellinit([0, lift(Mod(-rB,p))], p));
  [tA, tB, 1]
};

print("================================================================");
print("15-pairing sweep for all 7 canonical classes over F_43");
print("Target: some pairing gives degenerate traces {+13,-13}");
print("================================================================");

hitcount = 0;
for(i=1,7, {
  my(aa,bb,disc,sq,r1,r2,aroots,broots,W6,hits,pi,m,W,res,tag);
  aa=ca[i]; bb=cb[i];
  disc = lift(Mod(aa^2-4*bb,p));
  sq = lift(Mod(disc,p)^((p+1)/4));
  r1 = lift(Mod((-aa+sq)*inv2,p));
  r2 = lift(Mod((-aa-sq)*inv2,p));
  aroots = cube_roots_fp3(r1);
  broots = cube_roots_fp3(r2);
  W6 = [aroots[1],aroots[2],aroots[3],broots[1],broots[2],broots[3]];
  print("Class ",i," (a=",aa," b=",bb,")  r1=",r1,"  r2=",r2);
  hits = 0;
  for(pi=1,15,
    m = M15[pi];
    W = [W6[m[1][1]],W6[m[1][2]],W6[m[2][1]],W6[m[2][2]],W6[m[3][1]],W6[m[3][2]]];
    res = richelot_deg_traces(W);
    if(res[3]==1,
      tag = "";
      if(abs(res[1])==t_target && abs(res[2])==t_target && res[1]*res[2]<0, tag = " <== HOWE PAIR (traces +-13) !!!"; hits++; hitcount++);
      print("  matching ",pi," ",m,": trace_A=",res[1]," trace_B=",res[2],tag)
    )
  );
  if(hits==0, print("  (no degenerate matching hits {+13,-13} for class ",i,")"));
  print("")
});

print("================================================================");
print("SUMMARY: ", hitcount, " total (class,matching) hits");
print("Done.");
print("================================================================");
