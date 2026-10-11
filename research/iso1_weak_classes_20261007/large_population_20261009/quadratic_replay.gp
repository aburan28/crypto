\\ Replay fixed quadratic inputs across GP versions; regenerate no sample seeds.
element(v,z)=sum(i=1,#v,eval(v[i])*z^(i-1));
scope_replay(modulus,av,trace,order,dv,shortav,shortbv)={
  my(p=7,q=49,Q=117649,poly='z^6+sum(i=1,6,eval(modulus[i])*'z^(i-1)),z=ffgen(Mod(poly,p),'z),alpha=element(av,z),d=element(dv,z));
  if(alpha^q==alpha||d^q!=d||issquare(d),error("invalid published quadratic-family input"));
  my(f=('x^2-d)*('x-alpha)*('x-alpha^q),c3=subst(deriv(f),'x,alpha),c2=subst(deriv(deriv(f)),'x,alpha)/2,c1=subst(deriv(deriv(deriv(f))),'x,alpha)/6,E=ellinit([0,c2,0,c3*c1,c3^2],z),total=2);
  if(E.a4-E.a2^2/3!=element(shortav,z)||E.a6-E.a2*E.a4/3+2*E.a2^3/27!=element(shortbv,z),error("quadratic model identity mismatch"));
  for(code=0,Q-1,my(xx=sum(j=0,5,(code\p^j)%p*z^j),v=subst(f,'x,xx));total+=if(v==0,1,if(issquare(v),2,0)));
  if(total!=order||ellcard(E)!=order||Q+1-order!=trace,error("quadratic exact count replay failed"));
  print("QUADRATIC_REPLAY|trace=",trace,"|exact_order=",order);
};
read(Str(getenv("ISO1_POPULATION_STUDY"),"/quadratic_inputs.gp"));
