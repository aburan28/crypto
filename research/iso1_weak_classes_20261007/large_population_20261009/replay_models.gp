\\ Reconstruct from recorded moduli and roots, independently of generation seeds.
FIELDS=Map();MODEL_COUNT=0;ROUTE_COUNT=0;H90_COUNT=0;setrand(2026100944);
element(v,z)=sum(i=1,#v,eval(v[i])*z^(i-1));
cv(v,k)=vector(k,i,Str(lift(polcoef(v.pol,i-1))));
field_generator(p,k,c)={
  my(key=Str(p,"/",k,"/",c),z);
  if(mapisdefined(FIELDS,key,&z),return(z));
  my(poly='z^k+sum(i=1,k,eval(c[i])*'z^(i-1)));
  if(!isprime(p)||!polisirreducible(Mod(poly,p)),error("invalid recorded field"));
  z=ffgen(Mod(poly,p),'z);mapput(FIELDS,key,z);return(z);
};
root_curve(r,z)=ellinit([0,-vecsum(r),0,r[1]*r[2]+r[1]*r[3]+r[2]*r[3],-r[1]*r[2]*r[3]],z);
norm_test(r,p,k)={my(N=(p^k-1)/(p^2-1),a=(r[2]-r[1])^N,b=(r[3]-r[1])^N,c=(r[3]-r[2])^N);return(a==b||c==-a||c==b);};
canonical_sqrt(v,p,k)={my(h=sqrt(v));for(i=0,k-1,my(a=lift(polcoef(h.pol,i)));if(a!=0,if(a>p/2,h=-h);break));return(h);};
literal_weak_map(E,r,p,k,z)={
  my(q=p^2,n=k/2,N=(p^k-1)/(q-1),perm=[]);
  for(i=1,3,my(other=select(j->j!=i,[1,2,3]));if((r[other[1]]-r[i])^N==(r[other[2]]-r[i])^N,perm=[i,other[1],other[2]];break));
  if(!#perm,error("no equal-norm root ordering"));
  my(origin=r[perm[1]],u=r[perm[2]]-origin,v=r[perm[3]]-origin,lambda=v/u,alpha=0*z);
  for(j=0,k-1,
    my(c=z^j,coefficient=1+0*z);alpha=c;
    for(i=1,n-1,coefficient=coefficient^q/lambda;alpha+=coefficient*c^(q^i));
    if(alpha!=0,break);
  );
  if(alpha==0,error("Hilbert-90 basis search failed"));
  if(alpha^q!=lambda*alpha,error("constructive Hilbert-90 check failed"));
  if(!issquare(u/alpha),my(d=2);while(kronecker(d,p)!=-1,d++);my(w=canonical_sqrt(d+0*z,p,k),a=0,ns=w);while(issquare(ns),a++;ns=w+a);if(ns^q!=ns,error("base-field scaling failed"));alpha*=ns);
  my(h=canonical_sqrt(u/alpha,p,k),W=ellinit([0,-alpha-alpha^q,0,alpha*alpha^q,0],z),change=[h,origin,0,0]);
  if(h^2*alpha!=u||h^2*alpha^q!=v||ellchangecurve(E,change)[1..5]!=W[1..5],error("literal weak-model coordinate change failed"));
  for(i=1,3,my(P=random(E),P1=ellchangepoint(P,change));if(!ellisoncurve(W,P1)||ellchangepointinv(P1,change)!=P,error("weak-model point transport failed")));
  H90_COUNT++;print("H90_MAP|p=",p,"|degree=",k,"|root_order=",perm,"|alpha=",cv(alpha,k),"|h=",cv(h,k),"|origin=",cv(origin,k));
};
verify_record(p,k,modulus,av,bv,jv,rv,order,expected_weak,expected_four)={
  my(z=field_generator(p,k,modulus),r=vector(3,i,element(rv[i],z)),E=root_curve(r,z),a=element(av,z),b=element(bv,z));
  if(cv(E.a4-E.a2^2/3,k)!=av||cv(E.a6-E.a2*E.a4/3+2*E.a2^3/27,k)!=bv||cv(E.j,k)!=jv,error("recorded short-model mismatch"));
  my(S=ellinit([a,b],z));if(S.j!=E.j,error("short-model j mismatch"));
  for(i=1,3,if(!ellisoncurve(E,[r[i],0*z])||!ellisoncurve(S,[r[i]+E.a2/3,0*z]),error("root-coordinate translation failed")));
  if(expected_weak>=0&&norm_test(r,p,k)!=expected_weak,error("direct norm label failed"));
  if(expected_weak==1,literal_weak_map(E,r,p,k,z));
  if(expected_four==1,
    if(!issquare(r[2]-r[1])||!issquare(r[3]-r[1])||!issquare(r[3]-r[2]),error("full-4 root differences failed"));
    for(i=1,2,my(other=select(j->j!=i,[1,2,3]),u=sqrt(r[i]-r[other[1]]),v=sqrt(r[i]-r[other[2]]),P=[r[i]+u*v,u*v*(u+v)]);if(!ellisoncurve(E,P)||ellmul(E,P,2)!=[r[i],0*z],error("independent halving replay failed")));
  );
  for(i=1,3,my(P=random(E));if(ellmul(E,P,order)!=[0],error("recorded order failed point replay")));
  MODEL_COUNT++;
};
verify_route(p,k,modulus,source_roots,target_roots,degree,kernelv,order)={
  my(z=field_generator(p,k,modulus),rs=vector(3,i,element(source_roots[i],z)),rt=vector(3,i,element(target_roots[i],z)),E=root_curve(rs,z),F=root_curve(rt,z),kernel=element(kernelv,z));
  if(degree==2,if(!setsearch(Set(rs),kernel),error("2-isogeny kernel not a source root")),if(degree!=3||subst(elldivpol(E,3),'x,kernel)!=0,error("3-isogeny kernel failed")));
  my(iso=ellisogeny(E,'x-kernel),actual=ellinit(iso[1],z));if(actual[1..5]!=F[1..5],error("saved isogeny target mismatch"));
  for(i=1,3,my(P=random(E),image=ellisogenyapply(iso[2],P));if(!ellisoncurve(F,image)||ellmul(F,image,order)!=[0],error("route image replay failed")));
  ROUTE_COUNT++;
};
read(Str(getenv("ISO1_POPULATION_STUDY"),"/model_inputs.gp"));
print("REPLAY_TOTAL|models=",MODEL_COUNT,"|routes=",ROUTE_COUNT);
print("H90_TOTAL|literal_weak_isomorphisms=",H90_COUNT);
