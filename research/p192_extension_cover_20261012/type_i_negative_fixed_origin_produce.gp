\\ Native PARI producer for fixed-origin representatives of the 16 negative
\\ P-192 Type I classes. The old A/B pair is swapped, avoiding a double
\\ intersection at the three-point origin.
read("research/p192_extension_cover_20261012/type_i_fixed_origin_native_input.gp");
v=varhigher("v",x);
P(z)=sum(i=1,#z,Mod(z[i],p)*x^(i-1));
polyvec(z,n)=vector(n,i,lift(polcoeff(z,i-1,x)));
normal(z)={for(j=1,#z,if(z[j]!=0,return(vector(#z,i,z[i]/z[j]))));error("zero cubic")};
section(c)=sum(i=1,#c,Mod(c[i],p)*x^exponents[i][1]*v^exponents[i][2]);
curve=sum(i=1,#terms,Mod(terms[i][3],p)*x^terms[i][1]*v^terms[i][2]);
originX=vector(3,i,Mod(anchors[i][1],p)/Mod(anchors[i][3],p));
originY=vector(3,i,Mod(anchors[i][2],p)/Mod(anchors[i][3],p));
hT=prod(i=1,3,x-originX[i]);
yT=polinterpolate(originX,originY,x);
if(subst(curve,v,yT)%hT!=0,error("origin off quartic"));
trials=[[1,1,1,1],[1,2,3,5],[1,3,7,11],[2,5,11,17],[1,0,0,0],[0,1,0,0]];
evalrow(h,yy,row,col)=polcoeff((x^exponents[col][1]*yy^exponents[col][2])%h,row,x);

residual(cubic,known,degree)={
  my(res,qr,h,common,yy);
  res=polresultant(curve,cubic,v);
  if(poldegree(res,x)!=12,return([]));
  qr=divrem(res,known,x);
  if(qr[2]!=0,return([]));
  h=qr[1]/pollead(qr[1]);
  if(poldegree(h,x)!=degree,return([]));
  if(poldegree(gcd(h,deriv(h,x)),x)!=0,return([]));
  if(poldegree(gcd(h,known),x)!=0,return([]));
  common=gcd(subst(curve,x,Mod(x,h)),subst(cubic,x,Mod(x,h)));
  if(poldegree(common,v)!=1,return([]));
  yy=lift(-polcoeff(common,0,v)/polcoeff(common,1,v));
  if(poldegree(yy,x)>=degree,return([]));
  if(subst(curve,v,yy)%h!=0||subst(cubic,v,yy)%h!=0,return([]));
  return([h,yy]);
};

attempt(k,t)={
  my(c,hA,yA,hB,yB,M,ker,gcoeff,G,U,hU,yU,M2,ker2,fcoeff,F,D);
  c=cases[k];
  hA=P(c[3]);yA=P(c[4]);hB=P(c[1]);yB=P(c[2]);
  if(poldegree(gcd(hB,hT),x)!=0,error(Str("base/origin overlap ",k)));
  M=matrix(6,10,i,j,if(i<=3,evalrow(hB,yB,i-1,j),evalrow(hT,yT,i-4,j)));
  if(matrank(M)!=6,error(Str("denominator interpolation rank ",k)));
  ker=matker(M);
  if(matsize(ker)[2]!=4,error("denominator kernel dimension"));
  gcoeff=normal(vector(10,i,sum(j=1,4,Mod(trials[t][j],p)*ker[i,j])));
  G=section(gcoeff);
  U=residual(G,hB*hT,6);
  if(#U==0,return([]));
  hU=U[1];yU=U[2];
  if(poldegree(gcd(hA,hU),x)!=0,return([]));
  M2=matrix(9,10,i,j,if(i<=3,evalrow(hA,yA,i-1,j),evalrow(hU,yU,i-4,j)));
  if(matrank(M2)!=9,return([]));
  ker2=matker(M2);
  fcoeff=normal(vector(10,i,ker2[i,1]));
  F=section(fcoeff);
  D=residual(F,hA*hU,3);
  if(#D==0,return([]));
  return([k,vector(10,i,lift(gcoeff[i])),polyvec(hU,7),polyvec(yU,6),vector(10,i,lift(fcoeff[i])),polyvec(D[1],4),polyvec(D[2],3),trials[t]]);
};

produce(k)={
  my(ans);
  for(t=1,#trials,
    ans=attempt(k,t);
    if(#ans,return(ans));
  );
  error(Str("negative fixed-origin reduction failed at base index ",k));
};

{
results=List();
for(k=1,16,
  listput(results,produce(k));
  print("negative fixed-origin candidate ",k,"/16");
);
output=fileopen("research/p192_extension_cover_20261012/type_i_negative_fixed_origin_candidates.gp","w");
filewrite(output,Str("negative_cases=",Vec(results),";"));
fileclose(output);
print("NEGATIVE_FIXED_ORIGIN_PRODUCE_PASS classes=16 cuts=32");
}
