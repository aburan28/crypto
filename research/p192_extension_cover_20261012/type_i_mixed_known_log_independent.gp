\\ Independent PARI replay of known-only support and both known-scalar laws.
\\ The GP query input is decimal-transcribed from the frozen Rust receipt.
read("research/p192_extension_cover_20261012/type_i_signed_base_gp_input.gp");
read("research/p192_extension_cover_20261012/type_i_mixed_known_log_gp_queries.gp");
E=ellinit([-3,2455155546008943817740293915197451784769108058161191238065],p);
G=[Mod(602046282375688656758213480587526111916698976636884684818,p),Mod(174050332293622031404857552280219410364023488927386650641,p)];
if(!ellisoncurve(E,G)||ellmul(E,G,r)!=[0],error("P-192 generator/order"));
if(#A!=64||#steered!=64||#fresh!=64,error("query count"));

char(P)={
  my(x,v);
  if(#P==1,return(0));
  x=P[1];v=Mod(0,p);
  for(k=1,7,v+=S[k]*x^(k-1));
  if(v==0,return(0));
  return(kronecker(lift(v),p));
};

{
K=concat(A,vector(64,i,r-A[i]));
if(#Set(K)!=128,error("known signed scalars not distinct"));
knownsums=List();
for(i=1,128,for(j=i,128,for(k=j,128,
  listput(knownsums,(K[i]+K[j]+K[k])%r)
)));
if(#knownsums!=357760,error("known-only triple count"));
knownset=Set(Vec(knownsums));
if(#knownset!=349696,error(Str("known-only target count ",#knownset)));

for(i=1,64,
  my(s,P,Q);
  s=steered[i][1];
  if(s<=0||s>=r||!setsearch(knownset,s),error(Str("steered scalar/support ",i)));
  P=ellmul(E,G,s);Q=ellmul(E,P,-2);
  if(#Q!=2||lift(Q[1])!=steered[i][2]||lift(Q[2])!=steered[i][3],error(Str("steered point ",i)));
  if(char(P)!=1||char(elladd(E,P,P))!=1,error(Str("steered character ",i)));
);

fresh_known_only_hits=0;
for(i=1,64,
  my(s,P);
  s=fresh[i];
  if(s<=0||s>=r||setsearch(Set(A),s),error(Str("fresh scalar ",i)));
  P=ellmul(E,G,s);
  if(char(P)!=1||char(elladd(E,P,P))!=1,error(Str("fresh character ",i)));
  if(setsearch(knownset,s),fresh_known_only_hits++);
);

print("MIXED_KNOWN_LOG_INDEPENDENT_PASS signed_known=128 known_triples=357760 known_targets=349696 steered=64 fresh=64 fresh_known_only_hits=",fresh_known_only_hits);
}
