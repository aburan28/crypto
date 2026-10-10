\\ Exact two-torus parameter enumeration; errors fail closed in the native driver.
cv(v,k)=vector(k,i,Str(lift(polcoef(v.pol,i-1))));
key(v)=Str(v.pol);
jleg(v)=256*(v^2-v+1)^3/(v^2*(v-1)^2);
gen_sub(z,N)={
  my(g=random(z)^((z.p^poldegree(z.mod)-1)/N),fac=factor(N));
  while(g==0||g==1||sum(i=1,matsize(fac)[1],g^(N/fac[i,1])==1),g=random(z)^((z.p^poldegree(z.mod)-1)/N));g
};
run_census()={
  my(p=eval(getenv("ISO1_P")),q=p^2,Q=p^6,np=q^2+q+1,nm=q^2-q+1,seed=202610090701+p);
  if(p<=3||!isprime(p),error("requires prime p>3"));
  setrand(seed);my(z=ffgen([p,6],'z),w=ffgen([p,12],'w),embed=ffembed(z,w),back=ffinvmap(embed));
  my(d=random(z)^np);while(d==0||issquare(d),d=random(z)^np);
  my(beta=sqrt(ffmap(embed,d)));if(beta^q!=-beta,error("nonsplit base-field pair failed"));
  print("FIELD|p=",p,"|q=",q,"|Q=",Q,"|seed=",seed,"|modulus=",vector(6,i,Str(lift(polcoef(z.mod,i-1)))),"|extension_modulus=",vector(12,i,Str(lift(polcoef(w.mod,i-1)))),"|d=",cv(d,6));
  my(gp=gen_sub(z,np),gm=gen_sub(w,nm),seenp=vector(np),seenm=vector(nm),plus=Map(),minus=Map(),jminus=Map(),cp=0,cm=0);
  for(e=1,np-1,
    if(seenp[e+1],next);cp++;my(lambda=gp^e,E=ellinit([0,-(1+lambda),0,lambda,0],z),t=Q+1-ellcard(E),orbit=List(),a=e);
    for(i=1,6,for(s=1,2,my(k=if(s==1,a,(-a)%np));if(!seenp[k+1],seenp[k+1]=1;listput(orbit,k)));a=(a*p)%np);
    if(lambda^np!=1||lambda==1||E.j!=jleg(lambda),error("cubic parameter check failed"));
    mapput(plus,abs(t),1);
    print("PLUS|e=",e,"|orbit_size=",#orbit,"|trace=",t,"|lambda=",cv(lambda,6),"|short_a=",cv(E.a4-E.a2^2/3,6),"|short_b=",cv(2*E.a2^3/27-E.a2*E.a4/3,6),"|j=",cv(E.j,6));
  );
  for(e=1,nm-1,
    if(seenm[e+1],next);cm++;my(gamma=gm^e,lambda=gamma^(q+1),alphaL=beta*(gamma+1)/(gamma-1),alpha=ffmap(back,alphaL),orbit=List(),a=e);
    for(i=1,12,if(!seenm[a+1],seenm[a+1]=1;listput(orbit,a));a=(a*p)%nm);
    if(#orbit!=12||alpha==[]||alpha^q==alpha||alphaL^Q!=alphaL||lambda^nm!=1,error("quadratic torus or field check failed"));
    my(f=('x^2-d)*('x-alpha)*('x-alpha^q),c3=subst(deriv(f),'x,alpha),c2=subst(deriv(deriv(f)),'x,alpha)/2,c1=subst(deriv(deriv(deriv(f))),'x,alpha)/6,E=ellinit([0,c2,0,c3*c1,c3^2],z),t=Q+1-ellcard(E));
    if(ffmap(embed,E.j)!=jleg(lambda)||t^2>4*Q||(Q+1-t)%4,error("quadratic invariant or order check failed"));
    mapput(minus,abs(t),1);for(i=0,5,mapput(jminus,key(E.j^(p^i)),abs(t)));
    print("MINUS|e=",e,"|orbit_size=",#orbit,"|trace=",t,"|alpha=",cv(alpha,6),"|short_a=",cv(E.a4-E.a2^2/3,6),"|short_b=",cv(E.a6-E.a2*E.a4/3+2*E.a2^3/27,6),"|j=",cv(E.j,6));
  );
  if(cp!=(p^4+3*p^2+8)/12||cm!=(p^4-p^2)/12,error("orbit formula failed"));
  my(brute=0);if(getenv("ISO1_BRUTE")=="1",
    for(code=0,Q-1,
      my(alpha=sum(j=0,5,(code\p^j)%p*z^j));if(alpha^q==alpha,next);
      my(ga=(ffmap(embed,alpha)+beta)/(ffmap(embed,alpha)-beta),lambda=ga^(q+1),jj=ffmap(back,jleg(lambda)));
      if(lambda^nm!=1||lambda==1||jj==[]||!mapisdefined(jminus,key(jj)),error("alpha completeness check failed"));brute++;
    );if(brute!=Q-q,error("alpha completeness cardinality failed"));
  );
  my(op=0,om=0,ob=0,zero=0,total=0);
  forstep(t=-2*p^3+4,2*p^3-4,4,
    if(t%p==0,next);total++;my(a=mapisdefined(plus,abs(t)),b=mapisdefined(minus,abs(t)));op+=a;om+=b;ob+=a&&b;zero+=!a&&!b;
    print("CLASS|trace=",t,"|plus=",a,"|minus=",b,"|combined=",a||b);
  );
  print("TWO_BRANCH_COMPLETE|p=",p,"|plus_orbits=",cp,"|minus_orbits=",cm,"|ordinary_classes=",total,"|plus_support=",op,"|minus_support=",om,"|overlap=",ob,"|combined_zero=",zero,"|brute_alpha=",brute);
};
run_census();
