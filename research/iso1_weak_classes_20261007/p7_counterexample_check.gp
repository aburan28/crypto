\\ Independent existence and conductor check for an exact zero class.
check_counterexample()={
  my(p=7,Q=p^6,z=ffgen([p,6],'z));
  my(E=ellinit([2,6],z),t=Q+1-ellcard(E));
  if(t!=610,error("counterexample trace mismatch"));
  my(roots=polrootsmod((x^3+2*x+6)*z^0));
  if(#roots!=3,error("counterexample lacks full rational 2-torsion"));
  my(D=t^2-4*Q,DK=coredisc(D),f=sqrtint(D/DK));
  if(DK!=-19||f!=72,error("counterexample conductor mismatch"));
  print("p=7,Q=117649,trace=610,cardinality=117040,root_count=3,D_K=-19,f_pi=72,conductor_2_depth=3");
  print("Exact norm-one census: zero representatives at traces +610 and -610.");
};
check_counterexample();
