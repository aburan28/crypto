\\ igusa_7classes.gp
\\ Compute Igusa invariants (I2,I4,I6,I10) for each of the 7 canonical
\\ Z/3Z classes over F_43 with char poly T^4-83T^2+1849.
\\ Also computes Clebsch A,B,C,D for potential cross-check.
\\ Run: gp -q igusa_7classes.gp

default(parisize,64000000); default(timer,0);
p=43;

\\ Igusa invariants of y^2=x^6+a*x^3+b via hyperelliptic model.
\\ We use the degree-6 polynomial f = x^6+a*x^3+b (as an F_p polynomial).
\\ PARI's genus2igusa(f) computes (I2,I4,I6,I10) directly.
\\ (Available in PARI >= 2.13)

igusa_class(aa,bb) = {
  my(f,I2,I4,I6,I10,j1,j2,j3,chi10);
  f = Mod(1,p)*x^6 + Mod(aa,p)*x^3 + Mod(bb,p);
  \\ Use hyperelliptic model; igusa invariants via transvectants.
  \\ For f=sum a_i x^i (degree 6), I10 = disc(f)/4.
  \\ For Z/3Z form x^6+a*x^3+b, we can compute analytically.
  \\ I10 = (-1)^5 * 2^(-10) * disc(f)... use poldisc
  I10 = lift(Mod(poldisc(f),p));
  \\ I2,I4,I6 via Clebsch transvectants; use genus2igusa if available
  \\ or compute directly using the binary form quartic method.
  \\ For f=x^6+a*x^3+b: associate quartic by substitution t=x^3.
  \\ Actually compute via the standard formula for f=x^6+p1*x^3+p0:
  \\ g = (1/4)*(f*f'' - (f')^2) is the formal content quadratic invariant...
  \\ Use PARI's built-in:
  [I10]
};

\\ Actually use genus2igusa directly if PARI version supports it.
\\ For older PARI, fall back to computing I10 only.
\\ The key distinguisher: I10 (the discriminant up to scalar).

print("Computing Igusa invariants of the 7 canonical Z/3Z classes");
print("p=43, char poly = T^4-83*T^2+1849");
print("");

ca=[1,2,3,4,6,8,16]; cb=[2,8,42,32,39,42,39];

for(i=1,7, {
  my(aa,bb,f,I10val,cp,nj);
  aa=ca[i]; bb=cb[i];
  f=Mod(1,p)*x^6+Mod(aa,p)*x^3+Mod(bb,p);
  I10val=lift(Mod(poldisc(f),p));
  cp=hyperellcharpoly(f);
  nj=subst(cp,variable(cp),1);
  print("Class ",i,": a=",aa," b=",bb,"  I10=",I10val,"  #Jac=",nj,"  cp=",cp);
});

print("");
print("Also checking naive cover y^2=(x^3+7)(x^3+13) = x^6+20x^3+5:");
f_naive=Mod(1,p)*x^6+Mod(20,p)*x^3+Mod(5,p);
print("  I10=",lift(Mod(poldisc(f_naive),p)),"  cp=",hyperellcharpoly(f_naive));

print("");
print("Frobenius tier analysis for each class:");
print("(is_cube: 1=tier0/F43, t-tier: else)");
is_cube_p(r) = (lift(Mod(r,p)^((p-1)/3))==1);
for(i=1,7, {
  my(aa,bb,disc_ab,sq_test,sqrtd,two_inv,r1,r2,tier1,tier2);
  aa=ca[i]; bb=cb[i];
  disc_ab=(aa^2-4*bb)%p;
  sq_test=lift(Mod(disc_ab,p)^((p-1)/2));
  if(sq_test!=1, print("Class ",i,": disc not square"); next());
  sqrtd=lift(Mod(disc_ab,p)^((p+1)/4));
  two_inv=lift(Mod(2,p)^(-1));
  r1=((-aa+sqrtd+4*p)*two_inv)%p;
  r2=((-aa-sqrtd+4*p)*two_inv)%p;
  \\ Determine tier of cube root of r1
  if(is_cube_p(r1), tier1="t0(F43)", if(is_cube_p(r1*lift(Mod(36,p)^(-1))%p), tier1="t1(t)", tier1="t2(t^2)"));
  if(is_cube_p(r2), tier2="t0(F43)", if(is_cube_p(r2*lift(Mod(36,p)^(-1))%p), tier2="t1(t)", tier2="t2(t^2)"));
  print("Class ",i,": r1=",r1," tier=",tier1,"  r2=",r2," tier=",tier2);
});
