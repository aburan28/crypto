\\ thread11_weil_sweep.gp
\\ Sweep primes p<=PLIMIT with p%6==1, compute Weil poly for pair (g^1,g^2).
\\ Output CSV: p,a1,a2
\\ Usage: echo "PLIMIT=1000" | cat - thread11_weil_sweep.gp | gp -q

if(!defined(PLIMIT), PLIMIT=600);

forprime(p=7,PLIMIT,
  if(p%6!=1,next());
  g=znprimroot(p);
  b1=lift(g^1);
  b2=lift(g^2);
  t=ffgen(p,'t);
  P=t^0*x^6+(b1+b2)*t^0*x^3+(b1*b2%p)*t^0;
  f=hyperellcharpoly(P);
  a1_neg=polcoeff(f,3);
  a2=polcoeff(f,2);
  a1=-a1_neg;
  print(p,",",a1,",",a2)
);
quit();
