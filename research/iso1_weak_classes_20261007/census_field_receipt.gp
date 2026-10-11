\\ Reproduce the retained constructor's explicit field moduli in PARI.
field_receipt(p)={
  my(w=2);
  while(kronecker(w,p)!=-1,w++);
  my(u=ffgen(Mod(1,p)*(x^2-w),'u),k=1,s,zeta);
  while(1,
    s=(k%p)+(k\p)*u;
    if(s^((p^2-1)/3)!=1,break);
    k++;
  );
  zeta=s^((p^2-1)/3);
  if(zeta^3!=1||zeta==1,error("noncube or Frobenius multiplier mismatch"));
  print("p=",p,",u2_minus_w=",w,",theta3_minus_s_coefficients=[",k%p,",",k\p,"],zeta=",zeta);
};
for(i=1,6,field_receipt([7,13,53,59,61,199][i]));
