\\ Independent PARI elliptic-group replay of all source-line quotients.
read("research/p192_extension_cover_20261012/type_i_line_full_pencil_gp_input.gp");
E=ellinit([0,0,0,-3,Mod(0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1,p)],p);
point(z)=[Mod(z[1],p),Mod(z[2],p)];
if(#oldPoints!=141 || #newPoints!=9784 || #sourceLines!=4892,error("input sizes"));
section=elladd(E,elladd(E,point(oldPoints[sourceLines[1][1]]),point(oldPoints[sourceLines[1][2]])),elladd(E,point(newPoints[1]),point(newPoints[2])));
stars=vector(141,i,0);
for(k=1,4892,{
  my(i,j,x,y,s,ti,tj);
  i=sourceLines[k][1];j=sourceLines[k][2];
  if(i<1 || j>141 || i>=j,error(Str("source indices ",k)));
  x=point(newPoints[2*k-1]);y=point(newPoints[2*k]);
  if(!ellisoncurve(E,x) || !ellisoncurve(E,y),error(Str("new points ",k)));
  s=elladd(E,elladd(E,point(oldPoints[i]),point(oldPoints[j])),elladd(E,x,y));
  if(s!=section,error(Str("line section ",k)));
  ti=elladd(E,point(oldPoints[j]),elladd(E,x,y));
  tj=elladd(E,point(oldPoints[i]),elladd(E,x,y));
  if(stars[i]==0,stars[i]=ti,if(stars[i]!=ti,error(Str("star i ",k))));
  if(stars[j]==0,stars[j]=tj,if(stars[j]!=tj,error(Str("star j ",k))));
  if(k%512==0,print("FULL_PENCIL_GP_PROGRESS ",k,"/4892"));
});
for(i=1,141,if(stars[i]==0 || !ellisoncurve(E,point(oldPoints[i])),error(Str("old point or star ",i))));
print("FULL_PENCIL_GP_REPLAY_PASS lines=4892 stars=141 prime_bits=192");
