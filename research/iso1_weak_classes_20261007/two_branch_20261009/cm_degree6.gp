\\ Constant-degree quadratic CM intersections; class-polynomial construction is separate.
norm_remainder(N,F)={
  my(y=Mod(Mod(1,7)*'Y,F),one=Mod(Mod(1,7),F),m=(N-1)/2,M=[y,-one;one,0],v=M^(m-1)*[y,one]~);lift(v[1]+v[2])
};
run_microfactor()={
  my(p=7,q=p^2,Q=p^6,nm=q^2-q+1,maxdeg=0,totalfactors=0);
  for(i=1,5,
    my(t=[10,38,610,674,682][i],dk=coredisc(t^2-4*Q),f=sqrtint((t^2-4*Q)/dk),depth=valuation(f,2),div=divisors(f),combined=Mod(1,p)+0*'Y);
    for(ci=1,#div,
      my(c=div[ci]);if(valuation(c,2)!=depth,next);my(H=polclass(dk*c^2),fac=factor(Mod(1,p)*H));
      for(j=1,matsize(fac)[1],
        my(h=fac[j,1],d=poldegree(h));if(6%d,error("CM factor degree does not divide six"));if(d!=6,next);my(A=256*('Y-1)^3,B='Y-2,F=sum(k=0,d,polcoef(h,k)*A^k*B^(d-k)),C=norm_remainder(nm,F),g=gcd(F,C));
        maxdeg=max(maxdeg,poldegree(F));totalfactors++;combined=combined*g/gcd(combined,g);
      );
    );
    if(type(combined)!="t_POL",error("unexpected union polynomial"));print("MICROFACTOR_CLASS|trace=",t,"|minus_parameters=",2*poldegree(combined));
  );
  if(maxdeg>18,error("constant degree bound failed"));print("ALGEBRA_COMPLETE|CM_microfactors=",totalfactors,"|largest_decision_degree=",maxdeg);
};
run_microfactor();
