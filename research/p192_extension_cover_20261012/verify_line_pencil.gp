\\ Independent PARI elliptic-group replay of the 32-line quotient pencil.
read("research/p192_extension_cover_20261012/type_i_line_pencil_gp_input.gp");
E=ellinit([0,0,0,-3,Mod(0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1,p)],p);
point(z)=[Mod(z[1],p),Mod(z[2],p)];
if(#oldPoints!=141 || #newPoints!=64 || #sourceLines!=32,error("input sizes"));
for(k=1,32,{
  if(sourceLines[k][1]!=1,error("common point"));
  if(!ellisoncurve(E,point(oldPoints[sourceLines[k][2]])),error("old point"));
  if(!ellisoncurve(E,point(newPoints[2*k-1])),error("first new point"));
  if(!ellisoncurve(E,point(newPoints[2*k])),error("second new point"));
});
targets=vector(32,k,elladd(E,elladd(E,point(oldPoints[sourceLines[k][2]]),point(newPoints[2*k-1])),point(newPoints[2*k])));
for(k=2,32,if(targets[k]!=targets[1],error(Str("three-sum ",k))));
M=matrix(31,205,i,j,0);
for(k=2,32,{
  M[k-1,sourceLines[k][2]]+=1;
  M[k-1,141+2*k-1]+=1;
  M[k-1,141+2*k]+=1;
  M[k-1,sourceLines[1][2]]-=1;
  M[k-1,142]-=1;
  M[k-1,143]-=1;
});
if(matrank(matrix(31,205,i,j,Mod(M[i,j],r)))!=31,error("collision rank"));
print("PENCIL_GP_REPLAY_PASS lines=32 triples=32 rank=31 prime_bits=192");
