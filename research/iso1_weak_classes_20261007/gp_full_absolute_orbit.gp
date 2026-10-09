\\ Exact PARI/GP point counts on the norm-one torus, one per full
\\ absolute-Frobenius/inversion orbit. ISO1_P must be an odd prime >3.
\\ Each output row is positive_trace,orbit_size. Quadratic twisting
\\ supplies both trace signs, so the positive-trace row gets this weight.
p=eval(getenv("ISO1_P"));if(p<=3||!isprime(p),error("ISO1_P must be a prime >3"));q=p^2;n=p^6;N=q^2+q+1;z=ffgen([p,6],'z);g=ffprimroot(z)^((n-1)/N);lam=g;for(k=1,N-1,orb=vector(12);for(v=0,5,e=lift(Mod(k*p^v,N));orb[2*v+1]=e;orb[2*v+2]=lift(Mod(-e,N)));if(k==vecmin(orb),E=ellinit([0,-(1+lam),0,lam,0],z);t=abs(n+1-ellcard(E));print(t,",",#Set(orb)));lam*=g)
