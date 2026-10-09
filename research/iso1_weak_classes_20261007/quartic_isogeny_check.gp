\\ Independent finite-field check of the algebraic fourth-power proof.
setrand(1502);p=11;q=p^2;n=p^6;z=ffgen([p,6],'z);bad_fourth=0;
for(i=1,100,lam=random(z)^4;while(lam==0||lam==1,lam=random(z)^4);E=ellinit([0,-(1+lam),0,lam,0],z);if(ellcard(E)%16!=0,bad_fourth++));
print("p=11 random fourth-power parameters: n=100 nondivisible_by_16=",bad_fourth);

setrand(1516);bad_norm=0;
for(i=1,100,a=random(z);while(a==0,a=random(z));lam=a^(q-1);while(lam==1,a=random(z);lam=a^(q-1));mu=lam^((q^2+q+2)/4);E=ellinit([0,-(1+lam),0,lam,0],z);E2=ellinit([0,2*(1+lam),0,(1-lam)^2,0],z);c1=ellcard(E);c2=ellcard(E2);if(mu^4!=lam||c1!=c2||c2%16!=0,bad_norm++));
print("p=11 norm-one 2-isogeny checks: n=100 failures=",bad_norm)
