\\ Independent native replay of the 64 guided negative-class cubic sections.
\\ The candidate file supplies coefficients; this script recomputes every
\\ resultant, residual coordinate, interpolation rank, and class identity.
read("research/p192_extension_cover_20261012/type_i_fixed_origin_native_input.gp");
read("research/p192_extension_cover_20261012/type_i_guided_negative_fixed_origin_candidates.gp");
v=varhigher("v",x);
P(z)=sum(i=1,#z,Mod(z[i],p)*x^(i-1));
monic(z)=z/pollead(z);
section(c)=sum(i=1,#c,Mod(c[i],p)*x^exponents[i][1]*v^exponents[i][2]);
curve=sum(i=1,#terms,Mod(terms[i][3],p)*x^terms[i][1]*v^terms[i][2]);
originX=vector(3,i,Mod(anchors[i][1],p)/Mod(anchors[i][3],p));
originY=vector(3,i,Mod(anchors[i][2],p)/Mod(anchors[i][3],p));
hT=prod(i=1,3,x-originX[i]);
yT=polinterpolate(originX,originY,x);
if(poldegree(hT,x)!=3||subst(curve,v,yT)%hT!=0,error("origin"));

checkPoint(h,yy,label)={
  if(poldegree(h,x)<1||poldegree(yy,x)>=poldegree(h,x),error(label));
  if(subst(curve,v,yy)%h!=0,error(Str(label," curve")));
  if(poldegree(gcd(h,deriv(h,x)),x)!=0,error(Str(label," squarefree")));
};

checkCut(cubic,known,expectedH,expectedY,label)={
  my(res,qr,h,common,root);
  res=polresultant(curve,cubic,v);
  if(poldegree(res,x)!=12,error(Str(label," resultant degree")));
  qr=divrem(res,known,x);
  if(qr[2]!=0,error(Str(label," known factor")));
  h=monic(qr[1]);
  if(h!=expectedH,error(Str(label," residual x")));
  if(poldegree(gcd(h,known),x)!=0,error(Str(label," overlap")));
  checkPoint(h,expectedY,label);
  if(subst(cubic,v,expectedY)%h!=0,error(Str(label," cubic point")));
  common=gcd(subst(curve,x,Mod(x,h)),subst(cubic,x,Mod(x,h)));
  if(poldegree(common,v)!=1,error(Str(label," common y degree")));
  root=-polcoeff(common,0,v)/polcoeff(common,1,v);
  if(root!=Mod(expectedY,h),error(Str(label," residual y")));
};

{
if(#negative_cases!=64,error("guided negative candidate count"));
for(k=17,80,
  my(c,s,hA,yA,hB,yB,G,hU,yU,F,hD,yD,M,ker,Fgen);
  c=negative_cases[k-16];s=cases[k];
  if(c[1]!=k,error("candidate index"));
  hA=P(s[3]);yA=P(s[4]);hB=P(s[1]);yB=P(s[2]);
  G=section(c[2]);hU=P(c[3]);yU=P(c[4]);
  F=section(c[5]);hD=P(c[6]);yD=P(c[7]);
  checkPoint(hA,yA,"negative A");checkPoint(hB,yB,"negative B");
  if(poldegree(hA,x)!=3||poldegree(hB,x)!=3,error("input degree"));
  if(poldegree(gcd(hB,hT),x)!=0,error("origin overlap"));
  if(subst(G,v,yB)%hB!=0||subst(G,v,yT)%hT!=0,error("denominator base"));
  checkCut(G,hB*hT,hU,yU,"denominator" );
  if(poldegree(hU,x)!=6||poldegree(gcd(hA,hU),x)!=0,error("intermediate"));
  M=matrix(9,10,i,j,if(i<=3,
      polcoeff((x^exponents[j][1]*yA^exponents[j][2])%hA,i-1,x),
      polcoeff((x^exponents[j][1]*yU^exponents[j][2])%hU,i-4,x)));
  if(matrank(M)!=9,error("numerator interpolation rank"));
  ker=matker(M);Fgen=section(vector(10,i,ker[i,1]));
  if(Fgen/F!=pollead(Fgen)/pollead(F),error("numerator interpolation"));
  if(subst(F,v,yA)%hA!=0||subst(F,v,yU)%hU!=0,error("numerator base"));
  checkCut(F,hA*hU,hD,yD,"numerator");
  if(poldegree(hD,x)!=3,error("representative degree"));
  print("guided negative fixed-origin replay ",k-16,"/64");
);
print("GUIDED_NEGATIVE_FIXED_ORIGIN_VERIFY_PASS classes=64 cuts=128 prime_bits=192");
}
