\\ Independent CM-order implementation of the quadratic branch criterion.
jleg(v)=256*(v^2-v+1)^3/(v^2*(v-1)^2);
run_cm_controls()={
  my(p=7,q=p^2,Q=p^6,nm=q^2-q+1,z=ffgen(Mod('z^6+5*'z^5+5*'z^4+3*'z^3+5*'z^2+6*'z+5,p),'z),w=ffgen([p,12],'w),emb=ffembed(z,w),NmPoly=('x^nm-1)/('x-1));
  for(i=1,5,
    my(t=[10,38,610,674,682][i],delta=t^2-4*Q,dk=coredisc(delta),f=sqrtint(delta/dk),depth=valuation(f,2),div=divisors(f),total=0,roots_total=0);
    for(ci=1,#div,
      my(c=div[ci]);if(valuation(c,2)!=depth,next);
      my(H=polclass(dk*c^2),num=256*('x^2-'x+1)^3,den='x^2*('x-1)^2,h=poldegree(H),composition=sum(j=0,h,polcoef(H,j)*num^j*den^(h-j)),g=gcd(Mod(1,p)*composition,Mod(1,p)*NmPoly),deg=poldegree(g),js=polrootsmod(Mod(1,p)*H,z));
      if(h!=qfbclassno(dk*c^2),error("CM polynomial degree mismatch"));total+=deg;roots_total+=#js;
      my(groots=polrootsmod(g,w));if(#groots!=deg,error("quadratic CM gcd did not split"));
      for(j=1,#groots,if(groots[j]^nm!=1||groots[j]==1||subst(ffmap(emb,Mod(1,p)*H),'x,jleg(groots[j]))!=0,error("quadratic CM parameter verification failed")));
      print("CM_ORDER|trace=",t,"|D_K=",dk,"|f_pi=",f,"|conductor=",c,"|class_number=",h,"|j_roots=",#js,"|minus_parameters=",deg);
    );
    print("CM_CLASS|trace=",t,"|minus_parameters=",total,"|weak=",total>0,"|selected_j_roots=",roots_total);
  );
  print("ALGEBRA_COMPLETE|quadratic_CM_classes=5");
};
run_cm_controls();
