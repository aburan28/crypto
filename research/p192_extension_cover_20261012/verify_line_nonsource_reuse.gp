\\ Native replay of the first 128 frozen non-source quartic line cuts.
read("research/p192_extension_cover_20261012/type_i_line_nonsource_reuse_old141_gp_input.gp");
z='z;
if(#cases!=128,error("case count"));
for(k=1,#cases,{
  my(c,P,Q,R,S,m,d,a,b,ff,known,qr);
  c=cases[k];P=c[3];Q=c[4];R=c[5];S=c[6];
  if(P==Q || R==S || P==R || P==S || Q==R || Q==S,error(Str("distinct points ",k)));
  d=Mod(Q[1]-P[1],p);
  if(d==0,error(Str("vertical pair ",k)));
  m=Mod(Q[2]-P[2],p)/d;
  for(j=4,6,if(Mod(c[j][2]-P[2],p)!=m*Mod(c[j][1]-P[1],p),error(Str("line incidence ",k))));
  ff=sum(j=1,#terms,Mod(terms[j][3],p)*(Mod(P[1],p)+z)^terms[j][1]*(Mod(P[2],p)+m*z)^terms[j][2]);
  if(poldegree(ff,z)!=4,error(Str("quartic degree ",k)));
  known=z*(z-d)*(z-Mod(R[1]-P[1],p))*(z-Mod(S[1]-P[1],p));
  qr=divrem(ff,known,z);
  if(qr[2]!=0 || poldegree(qr[1],z)!=0 || qr[1]==0,error(Str("complete cut ",k)));
});
print("NONSOURCE_REUSE_GP_REPLAY_PASS cases=128 prime_bits=192");
