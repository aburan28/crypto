\\ Independent positive-label control. Native PARI/GP ellcard uses a
\\ different finite-field representation and point-count algorithm.
for(j=1,2,p=[41,43][j];setrand(1501);q=p^2;n=p^6;z=ffgen([p,6],'z);for(i=1,5000,lam=random(z)^(q-1);while(lam==0||lam==1,lam=random(z)^(q-1));E=ellinit([0,-(1+lam),0,lam,0],z);print(p,",",n+1-ellcard(E))))
