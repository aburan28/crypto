\\ Independent PARI/GP replay of 64 two-pivot P-192 two-cubic reduction.
read("research/p192_extension_cover_20261012/type_i_two_pivot_reduced64_input.gp");
v=varhigher("v",x);
P(z)=sum(i=1,#z,Mod(z[i],p)*x^(i-1));
section(c)=sum(i=1,#c,Mod(c[i],p)*x^exponents[i][1]*v^exponents[i][2]);
curve=sum(i=1,#terms,Mod(terms[i][3],p)*x^terms[i][1]*v^terms[i][2]);
originX=vector(3,i,Mod(anchors[i][1],p)/Mod(anchors[i][3],p));
originY=vector(3,i,Mod(anchors[i][2],p)/Mod(anchors[i][3],p));
hT=prod(i=1,3,x-originX[i]);
yT=polinterpolate(originX,originY,x);
if(#cases!=64 || poldegree(hT,x)!=3 || subst(curve,v,yT)%hT!=0,error("origin or count"));

checkPoint(h,yy,label)={
  if(poldegree(h,x)<1 || poldegree(yy,x)>=poldegree(h,x),error(label));
  if(subst(curve,v,yy)%h!=0,error(Str(label," curve")));
  if(poldegree(gcd(h,deriv(h,x)),x)!=0,error(Str(label," squarefree")));
};

checkCut(cubic,known,expectedH,expectedY,label)={
  my(res,qr,h,common,root);
  res=polresultant(curve,cubic,v);
  if(poldegree(res,x)!=12,error(Str(label," resultant degree")));
  qr=divrem(res,known,x);
  if(qr[2]!=0,error(Str(label," known factor")));
  h=qr[1]/pollead(qr[1]);
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
for(k=1,#cases,
  my(c,hA,yA,hB,yB,G,hU,yU,F,hD,yD,M,ker,Fgen);
  c=cases[k];
  hA=P(c[1]); yA=P(c[2]); hB=P(c[3]); yB=P(c[4]);
  G=section(c[5]); hU=P(c[6]); yU=P(c[7]);
  F=section(c[8]); hD=P(c[9]); yD=P(c[10]);
  checkPoint(hA,yA,"A"); checkPoint(hB,yB,"B");
  if(poldegree(hA,x)!=3 || poldegree(hB,x)!=3,error("input degree"));
  if(subst(G,v,yB)%hB!=0 || subst(G,v,yT)%hT!=0,error("denominator base"));
  checkCut(G,hB*hT,hU,yU,"denominator");
  if(poldegree(hU,x)!=6 || poldegree(gcd(hA,hU),x)!=0,error("intermediate"));
  M=matrix(9,10,i,j,if(i<=3,
      polcoeff((x^exponents[j][1]*yA^exponents[j][2])%hA,i-1,x),
      polcoeff((x^exponents[j][1]*yU^exponents[j][2])%hU,i-4,x)));
  if(matrank(M)!=9,error("numerator interpolation rank"));
  ker=matker(M);
  Fgen=section(vector(10,i,ker[i,1]));
  if(Fgen/F!=pollead(Fgen)/pollead(F),error("numerator interpolation"));
  if(subst(F,v,yA)%hA!=0 || subst(F,v,yU)%hU!=0,error("numerator base"));
  checkCut(F,hA*hU,hD,yD,"numerator");
  if(poldegree(hD,x)!=3,error("representative degree"));
  if(polcoeff(yD,2,x)==0,error("collinear representative"));
  if(k%512==0,print("TWO_PIVOT_FIXED_ORIGIN_NATIVE_PROGRESS ",k,"/",#cases));
);
print("TWO_PIVOT_FIXED_ORIGIN_NATIVE_PASS cases=",#cases," prime_bits=192");
}
