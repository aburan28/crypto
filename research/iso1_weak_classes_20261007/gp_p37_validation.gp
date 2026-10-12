\\ Independent PARI/GP trace check for random norm-one Legendre curves.
\\ The chosen seed and field generator are fixed; ellcard is independent
\\ of the Rust two-point baby-step counter used in the main census.
setrand(1501);
p=37;
q=p^2;
n=p^6;
z=ffgen([p,6],'z);
for(i=1,100, lam=random(z)^(q-1); while(lam==0 || lam==1, lam=random(z)^(q-1)); E=ellinit([0,-(1+lam),0,lam,0],z); print(n+1-ellcard(E)));
