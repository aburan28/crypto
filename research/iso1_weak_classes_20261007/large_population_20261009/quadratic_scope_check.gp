\\ Verify the quadratic h branch separately from the cubic norm-one census.
cv(v,k)=vector(k,i,Str(lift(polcoef(v.pol,i-1))));
run_check()={
  my(p=7,q=p^2,Q=p^6,N=(Q-1)/(q-1),seed=2026100977);
  setrand(seed);my(z=ffgen([p,6],'z),d=random(z)^N);while(d==0||issquare(d),d=random(z)^N);
  if(d^q!=d||issquare(d),error("quadratic h must be nonsplit over the base field"));
  print("FIELD|p=",p,"|degree=6|Q=",Q,"|seed=",seed,"|modulus=",vector(6,i,Str(lift(polcoef(z.mod,i-1)))),"|d=",cv(d,6));
  my(samples=List(),admitted=0,excluded=0,trace610=0);
  for(i=1,100,
    my(alpha=random(z));while(alpha^q==alpha,alpha=random(z));
    my(eq='y^2-('x^2-d)*('x-alpha)*('x-alpha^q),E=ellinit(ellfromeqn(eq),z),order=ellcard(E),t=Q+1-order,pass=(t-(Q+1))%16==0||(t+(Q+1))%16==0);
    if(t^2>4*Q||order%4,error("Hasse or quartic divisibility check failed"));
    admitted+=pass;excluded+=!pass&&t%p!=0;trace610+=abs(t)==610;
    if((!pass&&t%p!=0&&excluded<=2)||abs(t)==610,listput(samples,[alpha,E,order,t,i]));
    print("QUARTIC|",i,"|",t,"|",order,"|",pass,"|",cv(alpha,6));
  );
  for(i=1,#samples,
    my(row=samples[i],alpha=row[1],E=row[2],total=2);
    for(code=0,Q-1,
      my(xx=sum(j=0,5,(code\p^j)%p*z^j),v=(xx^2-d)*(xx-alpha)*(xx-alpha^q));
      total+=if(v==0,1,if(issquare(v),2,0));
    );
    if(total!=row[3],error("independent quartic all-x count mismatch"));
    my(f=('x^2-d)*('x-alpha)*('x-alpha^q),c3=subst(deriv(f),'x,alpha),c2=subst(deriv(deriv(f)),'x,alpha)/2,c1=subst(deriv(deriv(deriv(f))),'x,alpha)/6,F=ellinit([0,c2,0,c3*c1,c3^2],z));
    if(c3==0||ellcard(F)!=total||F.j!=E.j,error("explicit quartic coordinate model failed"));
    for(j=1,3,
      my(xx=random(z),value=subst(f,'x,xx));while(xx==alpha||!issquare(value),xx=random(z);value=subst(f,'x,xx));my(yy=sqrt(value),X=c3/(xx-alpha),Y=c3*yy/(xx-alpha)^2);
      if(!ellisoncurve(F,[X,Y])||alpha+c3/X!=xx||c3*Y/X^2!=yy,error("quartic point isomorphism failed"));
    );
    my(delta=row[4]^2-4*Q,dk=coredisc(delta),fpi=sqrtint(delta/dk));
    print("ALL_X|sample=",row[5],"|trace=",row[4],"|ellcard=",row[3],"|quartic_count=",total,"|alpha=",cv(alpha,6),"|short_a=",cv(F.a4-F.a2^2/3,6),"|short_b=",cv(F.a6-F.a2*F.a4/3+2*F.a2^3/27,6),"|j=",cv(F.j,6),"|D_K=",dk,"|f_pi=",fpi);
    print("QUARTIC_MAP|sample=",row[5],"|c3=",cv(c3,6),"|c2=",cv(c2,6),"|c1=",cv(c1,6),"|formula=X=c3/(x-alpha);Y=c3*y/(x-alpha)^2");
  );
  print("QUADRATIC_SCOPE_COMPLETE|samples=100|cubic_admitted=",admitted,"|cubic_excluded=",excluded,"|absolute_trace610=",trace610,"|all_x_controls=",#samples);
};
run_check();
