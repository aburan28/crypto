\\ Independent PARI arithmetic, no solver or factor-base implementation reused.
allocatemem(100000000);
rows=readstr("research/conductor_gap_cm_orbits_20261015/prepared/target-replay.gpvec");
verify_target(v,id)={
  my(p=v[1],r=v[4],E,G,Q);
  if(!isprime(p)||!isprime(r)||r^2<=16*p||(p+1-r)^2>4*p,error("prime/Hasse prerequisites failed"));
  E=ellinit([v[2],v[3]],p);
  G=[Mod(v[5][1],p),Mod(v[5][2],p)];Q=[Mod(v[6][1],p),Mod(v[6][2],p)];
  if(E.disc==0||!ellisoncurve(E,G)||!ellisoncurve(E,Q)||ellmul(E,G,v[7])!=Q||ellmul(E,G,r)!=[0],error("target replay failed"));
  if(ellmul(E,G,(v[7]+1)%r)==Q,error("negative scalar mutation accepted"));
  print("{\"row\":",id,",\"target_scalar_verified\":true,\"hasse_unique_prime_order\":true,\"mutated_scalar_rejected\":true}");
};
for(i=1,#rows,verify_target(eval(rows[i]),i));
quit;
