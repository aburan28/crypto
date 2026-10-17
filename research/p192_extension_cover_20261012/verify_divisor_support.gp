\\ Independent PARI replay of the 160 normalized x-polynomials.
read("research/p192_extension_cover_20261012/type_i_divisor_support_input.gp");
verify()={
  my(rootcounts=vector(4,i,vector(4,j,0)),squarefree=0,shared=0,pairs=0);
  if(#X != 160, error("wrong class count"));
  if(P != 6277101735386680763835789423207666416083908700390324961279,
    error("wrong P-192 field prime"));
  for(i=1,#X,
    my(f=X[i],factors,roots=0,group);
    if(poldegree(f)!=3 || pollead(f)!=1, error("nonmonic cubic"));
    if(poldegree(gcd(Mod(1,P)*f,Mod(1,P)*deriv(f)))!=0,
      error("nonsquarefree cubic"));
    squarefree++;
    factors=factormod(f,P);
    for(k=1,matsize(factors)[1],
      if(factors[k,2]!=1, error("repeated factor"));
      if(poldegree(factors[k,1])==1, roots++);
    );
    if(roots!=0 && roots!=1 && roots!=3, error("unexpected cubic split"));
    group=if(i<=16,1,if(i<=80,2,if(i<=96,3,4)));
    rootcounts[group][roots+1]++;
  );
  for(i=1,#X-1,
    for(j=i+1,#X,
      pairs++;
      if(poldegree(gcd(Mod(1,P)*X[i],Mod(1,P)*X[j]))>0, shared++);
    );
  );
  print("DIVISOR_SUPPORT_GP_PASS classes=",#X,
        " squarefree=",squarefree,
        " pairs=",pairs,
        " shared_x_pairs=",shared,
        " root_counts=",rootcounts);
};
verify();
quit();
