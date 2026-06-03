\\ igusa_7classes_full.gp
\\ Full Igusa-Clebsch invariants for the 7 Z/3Z classes over F_43.
\\ Pre-sets parisize=256MB before read to avoid stack-reset abort.
\\ Run: gp -q igusa_7classes_full.gp

default(parisize,256000000); default(timer,0);
\\ Load transvectant engine (the parisize inside it matches, so no reset)
read("igusa_clebsch.gp");

\\ Clebsch A,B,C,D and Igusa (I2,I4,I6,I10) from igusa_clebsch_complete recipe
clebsch_ABCD_fn(h) = {
  my(hv,iv,Delta,y1,y2,y3,A,B,C,D);
  hv=poly_to_binary_form(h,6);
  iv=transvectant(hv,6,hv,6,4);
  A=binary_form_constant(transvectant(hv,6,hv,6,6));
  B=binary_form_constant(transvectant(iv,4,iv,4,4));
  Delta=transvectant(iv,4,iv,4,2);
  C=binary_form_constant(transvectant(iv,4,Delta,4,4));
  y1=transvectant(hv,6,iv,4,4);
  y2=transvectant(iv,4,y1,2,2);
  y3=transvectant(iv,4,y2,2,2);
  D=binary_form_constant(transvectant(y3,2,y1,2,2));
  [A,B,C,D]};

igusa_full(h,p) = {
  my(abcd,A,B,C,D,I2,I4,I6,I10);
  abcd=clebsch_ABCD_fn(h);
  A=abcd[1]; B=abcd[2]; C=abcd[3]; D=abcd[4];
  I2=lift(Mod(-120*A,p));
  I4=lift(Mod(-720*A^2+6750*B,p));
  I6=lift(Mod(8640*A^3-108000*A*B+202500*C,p));
  I10=lift(Mod(-62208*A^5+972000*A^3*B+1620000*A^2*C-3037500*A*B^2-6075000*B*C-4556250*D,p));
  [I2,I4,I6,I10]};

p=43;
ca=[1,2,3,4,6,8,16]; cb=[2,8,42,32,39,42,39];

print("");
print("Igusa invariants (I2,I4,I6,I10) mod 43 for 7 Z/3Z classes");
print("char poly = T^4-83T^2+1849 = (T^2-13T+43)(T^2+13T+43)");
print("I2 = 3a^2-120b mod p  (direct formula check)");
print("");

for(i=1,7,{
  my(aa,bb,f,ig,I2chk);
  aa=ca[i]; bb=cb[i];
  f=x^6+aa*x^3+bb;
  ig=igusa_full(f,p);
  I2chk=lift(Mod(3*aa^2-120*bb,p));
  print("Class ",i," (a=",aa," b=",bb,"): I2=",ig[1],"(",I2chk,") I4=",ig[2]," I6=",ig[3]," I10=",ig[4])
});

print("");
print("Naive cover (a=20,b=5): NOT in target isogeny class");
f0=x^6+20*x^3+5;
ig0=igusa_full(f0,p);
I2c=lift(Mod(3*20^2-120*5,p));
print("  I2=",ig0[1],"(",I2c,") I4=",ig0[2]," I6=",ig0[3]," I10=",ig0[4]);
