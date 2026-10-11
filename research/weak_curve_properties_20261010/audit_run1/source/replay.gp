\\ Independent finite-field properties for the frozen catalogue only.
default(parisizemax,2^31);setrand(202610102603);FIELDS=Map();COUNTER=0;
element(v,z)=sum(i=1,#v,v[i]*z^(i-1))+0*z;
cv(v,z,n)={v+=0*z;vector(n,i,lift(polcoef(v.pol,i-1)));};
field(p,n,m)={my(k=Str(p,n,m),z);if(mapisdefined(FIELDS,k,&z),return(z));if(!isprime(p),error("characteristic not prime"));my(f='z^n+sum(i=1,n,m[i]*'z^(i-1)));if(!polisirreducible(Mod(f,p)),error("modulus reducible"));z=ffgen(Mod(f,p),'z);mapput(FIELDS,k,z);z};
depth(t,Q)={if(t%2,return(0));my(k=1,r=0,a);while(k<=128,a=(t/2)%2^(k-1);if((Q-a*t+a^2)%2^(2*k)==0||(Q-(a+2^(k-1))*t+(a+2^(k-1))^2)%2^(2*k)==0,r=k;k++,break));r};
audit(i,p,n,m,av,bv,jv,t,N)={
  my(z=field(p,n,m),A=element(av,z),B=element(bv,z),J=element(jv,z),Q=p^n,E=ellinit([A,B],z));
  if(E.j!=J||N!=Q+1-t||t^2>4*Q,error("input consistency failed"));
  my(jd=0,cd=0,ds=divisors(n));for(k=1,#ds,my(d=ds[k]);if(!jd&&J^(p^d)==J,jd=d);if(!cd&&A^(p^d)==A&&B^(p^d)==B,cd=d));
  my(kind,tw,twist_order);if(J==0,kind="sextic";tw=B;twist_order=6,if(J==1728,kind="quartic";tw=A;twist_order=4,kind="quadratic";tw=3*B/(2*A);twist_order=2));
  my(iso=0);for(k=1,#ds,my(d=ds[k]);if(J^(p^d)==J&&tw^((Q-1)/gcd(twist_order,(Q-1)/(p^d-1)))==1,iso=d;break));if(!iso,error("missing definition degree"));
  my(roots=polrootsmod('x^3+A*'x+B,z),rvec=vecsort(vector(#roots,k,cv(roots[k],z,n))));
  my(f='x^3+A*'x+B,dr=divrem(elldivpol(E,4),f));if(dr[2]!=0||poldegree(dr[1])!=6,error("division polynomial quotient failed"));
  my(r4=polrootsmod(dr[1],z),full4=#r4==6);if(full4,for(k=1,#r4,if(!issquare(subst(f,'x,r4[k])),full4=0;break)));
  my(ordinary=t%p!=0,dk2=-1,f2=-1,split=-2,Delta=t^2-4*Q);if(ordinary,
    f2=depth(t,Q);my(localdisc=Delta/2^(2*f2));dk2=valuation(localdisc,2);split=kronecker(localdisc,2);
    if(abs(Delta)<2^40,my(D=quaddisc(Delta));if(valuation(Delta/D,2)!=2*f2||valuation(D,2)!=dk2||kronecker(D,2)!=split,error("full-discriminant cross-check failed")));
  );
  print("PROPERTY|index=",i,"|j_degree=",jd,"|coefficient_degree=",cd,"|isomorphic_degree=",iso,"|twist_kind=",kind,"|roots=",rvec,"|full4=",full4,"|f2=",f2,"|dk2=",dk2,"|split2=",split);
  COUNTER++;
};
read(getenv("WEAK_PROPERTIES_INPUTS"));
if(COUNTER!=814,error("incomplete model replay"));print("PROPERTY_REPLAY_COMPLETE|models=814");
