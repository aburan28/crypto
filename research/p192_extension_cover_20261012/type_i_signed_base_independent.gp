\\ Independent PARI elliptic arithmetic replay of the signed-base census.
\\ Input is decimal-transcribed from the three frozen JSON receipts by jq.
read("research/p192_extension_cover_20261012/type_i_signed_base_gp_input.gp");
E=ellinit([-3,2455155546008943817740293915197451784769108058161191238065],p);
G=[Mod(602046282375688656758213480587526111916698976636884684818,p),Mod(174050332293622031404857552280219410364023488927386650641,p)];
B=vector(16,i,[Mod(B[i][1],p),Mod(B[i][2],p)]);
V=concat(B,vector(16,i,ellmul(E,B[i],-1)));
for(i=1,32,if(!ellisoncurve(E,V[i]),error("signed point off curve")));
if(#Set(V)!=32,error("signed base collision"));
char(P)={
  my(x,v);
  if(#P==1,return(0));
  x=P[1];v=Mod(0,p);
  for(k=1,7,v+=S[k]*x^(k-1));
  if(v==0,return(0));
  return(kronecker(lift(v),p));
};
T=List();cancellations=0;
addtrip(i,j,k)={
  if((i+16==j||j+16==i||i+16==k||k+16==i||j+16==k||k+16==j),cancellations++);
  listput(T,elladd(E,elladd(E,V[i],V[j]),V[k]));
};
for(i=1,32,for(j=i,32,for(k=j,32,addtrip(i,j,k))));
if(#T!=5984||cancellations!=512,error("signed triple count"));
U=Set(Vec(T));
if(#U!=5504,error("signed support count"));
cells=vector(4,i,0);accepted=0;
checktarget(Q)={
  my(P,c1,c2,pos);
  P=ellmul(E,Q,-((r+1)/2));
  if(elladd(E,P,P)!=ellmul(E,Q,-1),error("known-scalar half equation"));
  c1=char(P);c2=char(elladd(E,P,P));
  if(c1==0||c2==0,error("unexpected branch target"));
  pos=1+2*(c1==-1)+(c2==-1);cells[pos]++;
  if(c1==1&&c2==1,accepted++);
};
for(i=1,#U,checktarget(U[i]));
if(cells!=[1344,1416,1392,1352]||accepted!=1344,error("character partition"));
hits=0;
checkquery(s)={
  my(P,Q);
  P=ellmul(E,G,s);
  if(char(P)!=1||char(elladd(E,P,P))!=1,error("accepted query character"));
  Q=ellmul(E,elladd(E,P,P),-1);
  if(setsearch(U,Q),hits++);
};
for(i=1,#A,checkquery(A[i]));
if(#A!=64||hits!=0,error("frozen query hits"));
print("SIGNED_BASE_INDEPENDENT_PASS triples=5984 distinct=5504 accepted=1344 cancel=512 queries=64 hits=0");
