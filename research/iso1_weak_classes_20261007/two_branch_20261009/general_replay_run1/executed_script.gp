\\ Independent recorded-coordinate replay; does not run the construction search.
FIELDS=Map();MODELS=Map();setrand(202610091111);
element(v,z)=sum(i=1,#v,eval(v[i])*z^(i-1));
cv(v,z)={v+=0*z;vector(poldegree(z.mod),i,Str(lift(polcoef(v.pol,i-1))));};
field(p,k,m)={my(key=Str(p,m),z);if(mapisdefined(FIELDS,key,&z),return(z));my(f='z^k+sum(i=1,k,eval(m[i])*'z^(i-1)));if(!polisirreducible(Mod(f,p)),error("recorded modulus reducible"));z=ffgen(Mod(f,p),'z);mapput(FIELDS,key,z);z};
check_model(p,k,m,av,bv,ord)={
  my(z=field(p,k,m),E=ellinit([element(av,z),element(bv,z)],z),key=Str(p,m,av,bv,ord));if(mapisdefined(MODELS,key),return(E));
  if(ellcard(E)!=ord,error("independent recorded model point count failed"));for(i=1,3,my(P=random(E));if(ellmul(E,P,ord)!=[0],error("independent recorded order replay failed")));
  mapput(MODELS,key,1);print("MODEL|p=",p,"|degree=",k,"|trace=",p^k+1-ord,"|order=",ord,"|modulus=",m,"|short_a=",av,"|short_b=",bv,"|j=",cv(E.j,z));E
};
check_route(p,k,m,as,bs,at,bt,kv,mode,degree,ord)={
  my(z=field(p,k,m),A=check_model(p,k,m,as,bs,ord),B=check_model(p,k,m,at,bt,ord),kernel=if(mode,sum(i=1,#kv,element(kv[i],z)*'x^(i-1)),'x-element(kv,z)));
  if(poldegree(kernel)!=if(degree==2,1,(degree-1)/2),error("recorded kernel degree mismatch"));my(iso=ellisogeny(A,kernel));if(ellinit(iso[1],z)[1..5]!=B[1..5],error("recorded isogeny coefficients failed"));
  for(i=1,3,my(P=random(A),P1=ellisogenyapply(iso[2],P));if(!ellisoncurve(B,P1),error("independent isogeny image failed"));print("POINT_MAP|degree=",degree,"|source=",if(#P==1,[],vector(2,j,cv(P[j],z))),"|target=",if(#P1==1,[],vector(2,j,cv(P1[j],z)))););
};
check_plus(p,k,m,av,bv,alpha,h,origin,ord)={
  my(z=field(p,k,m),A=check_model(p,k,m,av,bv,ord),a=element(alpha,z),s=element(h,z),r=element(origin,z),W=ellinit([0,-a-a^(p^2),0,a*a^(p^2),0],z),change=[s,r,0,0]);if(a^(p^2)==a||ellchangecurve(A,change)[1..5]!=W[1..5],error("recorded cubic conversion failed"));
  check_model(p,k,m,cv(W.a4-W.a2^2/3,z),cv(2*W.a2^3/27-W.a2*W.a4/3,z),ord);
  for(i=1,3,my(P=random(A),P1=ellchangepoint(P,change));if(!ellisoncurve(W,P1)||ellchangepointinv(P1,change)!=P,error("independent cubic transport failed")));
};
check_minus(p,k,m,av,bv,alpha,d,twist,h,c3,a2,ord)={
  my(z=field(p,k,m),A=check_model(p,k,m,av,bv,ord),a=element(alpha,z),dd=element(d,z),e=element(twist,z),s=element(h,z),c=element(c3,z),b2=element(a2,z),f=e*('x^2-dd)*('x-a)*('x-a^(p^2)),c2=subst(deriv(deriv(f)),'x,a)/2,c1=subst(deriv(deriv(deriv(f))),'x,a)/6,W=ellinit([0,c2,0,c*c1,e*c^2],z));
  if(dd^(p^2)!=dd||issquare(dd)||a^(p^2)==a||e^(p^2)!=e||subst(deriv(f),'x,a)!=c||W.a2!=b2,error("recorded quadratic coefficients failed"));
  my(wa=W.a4-W.a2^2/3,wb=W.a6-W.a2*W.a4/3+2*W.a2^3/27);if(s^4*wa!=A.a4||s^6*wb!=A.a6,error("recorded quadratic scaling failed"));check_model(p,k,m,cv(wa,z),cv(wb,z),ord);
  for(i=1,3,my(P=random(A));if(#P==1,next);my(X=P[1]/s^2-b2/3,Y=P[2]/s^3);if(X==0,next);my(x=a+c/X,y=c*Y/X^2);if(y^2!=subst(f,'x,x),error("independent quartic point transport failed"));print("QUARTIC_POINT|x=",cv(x,z),"|y=",cv(y,z)););
};
read(getenv("ISO1_INPUTS"));
