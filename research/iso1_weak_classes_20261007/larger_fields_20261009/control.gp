\\ One field cell; each printed record is an exact completed check.
coeff_vector(v,d)=vector(d,i,lift(polcoef(v.pol,i-1)));
run_cell()={
  my(p=eval(getenv("ISO1_P")),a=eval(getenv("ISO1_A")),d=eval(getenv("ISO1_N")),samples=eval(getenv("ISO1_SAMPLES")),seed=eval(getenv("ISO1_SEED")));
  if(!isprime(p)||p<=3||d<3||d%2!=1,error("invalid field cell"));
  my(q=p^a,Q=q^d,N=(Q-1)/(q-1),degree=a*d);
  if(q%4!=1,error("base field must equal 1 mod 4"));
  setrand(seed);
  my(z=ffgen([p,degree],'z),u=lift(1/Mod(4,N)));
  print("FIELD|",p,"|",a,"|",d,"|",Q,"|",vector(degree,i,lift(polcoef(z.mod,i-1))));
  for(i=1,samples,
    setrand(seed+i);
    my(start=getwalltime(),lam=random(z)^(q-1));
    while(lam==0||lam==1,lam=random(z)^(q-1));
    my(mu=lam^u,s=mu^2,E=ellinit([0,-(1+lam),0,lam,0],z),Ep=ellinit([0,2*(1+lam),0,(1-lam)^2,0],z));
    if(lam^N!=1||mu^4!=lam,error("norm or fourth root mismatch"));
    my(roots=[0*z,-(1+s)^2,-(1-s)^2]);
    for(j=1,2,
      my(other=select(k->k!=j,[1,2,3]),aa=sqrt(roots[j]-roots[other[1]]),bb=sqrt(roots[j]-roots[other[2]]));
      my(P=[roots[j]+aa*bb,aa*bb*(aa+bb)]);
      if(!ellisoncurve(Ep,P)||ellmul(Ep,P,2)!=[roots[j],0*z]||ellmul(Ep,P,4)!=[0],error("rational halving failed"));
    );
    my(ii=sqrt(-z^0),P0=[-s,-ii*s*(1+s)],X=P0[1]+E.a2+E.a4/P0[1],Y=P0[2]*(1-E.a4/P0[1]^2));
    if(!ellisoncurve(E,P0)||ellmul(E,P0,2)!=[0*z,0*z]||!ellisoncurve(Ep,[X,Y])||[X,Y]!=[roots[2],0*z],error("explicit isogeny point check failed"));
    print("GEOMETRY|",p,"|",i,"|",getwalltime()-start,"|",coeff_vector(lam,degree));
    if(getenv("ISO1_MODE")=="geometry",next);
    my(c0=getwalltime(),ne=ellcard(E),c1=getwalltime(),np=ellcard(Ep),c2=getwalltime(),t=Q+1-ne);
    if(ne!=np||np%16||t^2>4*Q,error("cardinality or Hasse failure"));
    print("COUNT|",p,"|",a,"|",d,"|",i,"|",ne,"|",t,"|",c1-c0,"|",c2-c1,"|",t%p!=0,"|",coeff_vector(E.a4-E.a2^2/3,degree),"|",coeff_vector(E.a6-E.a2*E.a4/3+2*E.a2^3/27,degree),"|",coeff_vector(E.j,degree),"|",coeff_vector(Ep.a4-Ep.a2^2/3,degree),"|",coeff_vector(Ep.a6-Ep.a2*Ep.a4/3+2*Ep.a2^3/27,degree),"|",coeff_vector(Ep.j,degree));
  );
  print("COMPLETE|",p,"|",a,"|",d,"|",samples);
};
run_cell();
