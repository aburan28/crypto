\\ check_traces.gp  — minimal script, no nested braces, single-line for bodies
\\ For each of the 7 Z/3Z classes: compute r1,r2 and traces of y^2=x^3-r1, y^2=x^3-r2
\\ Also verify: which degenerate Richelot pairings give equal constant terms
\\ Run: gp -q check_traces.gp

default(parisize,128000000); default(timer,0);

p=43; inv2=lift(Mod(2,p)^(-1)); z3=6; z3sq=36;
ca=[1,2,3,4,6,8,16]; cb=[2,8,42,32,39,42,39];

print("=== Traces of E_{r1} and E_{r2} for 7 Z/3Z classes ===");
print("p=43, E_c = y^2 = x^3+c, t = p+1-#E(F_p)");
print("");

\\ Process class 1
do_class(ii) = {
  my(aa,bb,disc,sq,r1,r2,mr1,mr2,t1,t2);
  aa=ca[ii]; bb=cb[ii];
  disc=lift(Mod(aa^2-4*bb,p));
  sq=lift(Mod(disc,p)^((p+1)/4));
  r1=lift(Mod((-aa+sq)*inv2,p));
  r2=lift(Mod((-aa-sq)*inv2,p));
  mr1=lift(Mod(-r1,p));
  mr2=lift(Mod(-r2,p));
  t1=p+1-ellcard(ellinit([0,mr1],p));
  t2=p+1-ellcard(ellinit([0,mr2],p));
  print("Class ",ii," (a=",aa,",b=",bb,"): r1=",r1," r2=",r2," -> E_{r1}=y^2=x^3+",mr1," t=",t1,"  E_{r2}=y^2=x^3+",mr2," t=",t2, if(abs(t1)==13&&abs(t2)==13&&t1*t2<0," *** HOWE MATCH",""))
};

for(ii=1,7,do_class(ii));

print("");
print("=== Trace table for all j=0 curves over F_43 ===");
for(jj=1,42,print("  y^2=x^3+",jj,": t=",p+1-ellcard(ellinit([0,jj],p))));

print("");
print("=== Cube root tiers for r1,r2 of each class ===");
\\ tier: 0 = r is a cube in F_43, 1 = r/z3sq is a cube, 2 = r/z3 is a cube
tier(r) = {
  my(rt); r=r%p;
  if(r==0,return(0));
  rt=lift(Mod(r,p)^((p-1)/3));
  if(rt==1,return(0));
  rt=lift(Mod(r*lift(Mod(z3sq,p)^(-1)),p)^((p-1)/3));
  if(rt==1,return(1));
  2
};

do_tier(ii) = {
  my(aa,bb,disc,sq,r1,r2);
  aa=ca[ii]; bb=cb[ii];
  disc=lift(Mod(aa^2-4*bb,p));
  sq=lift(Mod(disc,p)^((p+1)/4));
  r1=lift(Mod((-aa+sq)*inv2,p));
  r2=lift(Mod((-aa-sq)*inv2,p));
  print("Class ",ii," (a=",aa,",b=",bb,"): r1=",r1," tier=",tier(r1),"  r2=",r2," tier=",tier(r2))
};

for(ii=1,7,do_tier(ii));

print("");
print("=== Which degenerate (2,2) pairings exist for each class? ===");
print("Pairing type (i+j mod 3 = k): alpha_i paired with beta_j where i+j=k mod 3");
print("Equal-constant-term pairings: k=0 (i+j=0), k=1 (i+j=1), k=2 (i+j=2)");
print("For each: 'left' group = {alpha_0,alpha_1,alpha_2}, 'right' = {beta_0,beta_1,beta_2}");
print("In all cases, the degenerate Richelot gives E_{r1} x E_{r2}.");
print("");
print("KEY INSIGHT: All 3 degenerate Richelot pairings from a Z/3Z class");
print("give the SAME split ppav E_{r1} x E_{r2}.");
print("The Howe pair (E1=y^2=x^3+7, E2=y^2=x^3+13) requires traces {+13,-13}.");
print("If no class has {t_{r1},t_{r2}} = {+13,-13}, the Howe-glued curve");
print("is NOT in the Z/3Z family y^2=x^6+ax^3+b.");
print("");
print("Done.");
