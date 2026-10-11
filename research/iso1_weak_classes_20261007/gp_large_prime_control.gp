\\ Independent norm-one point-count samples for a completed larger prime.
\\ ISO1_P is recorded by the caller; output has no header.
p=eval(getenv("ISO1_P"));if(p<=3||!isprime(p),error("ISO1_P must be a prime >3"));
setrand(1501);q=p^2;n=p^6;z=ffgen([p,6],'z);
for(i=1,5000,lam=random(z)^(q-1);while(lam==0||lam==1,lam=random(z)^(q-1));E=ellinit([0,-(1+lam),0,lam,0],z);print(p,",",n+1-ellcard(E)));
