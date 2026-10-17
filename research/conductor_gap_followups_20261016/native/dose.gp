\\ Independently construct four chosen-trace conductor strata per class.
default(parisizemax, 2000000000);
{
for(ci=1,4,
 f=[35,55,65,77][ci]; t=134217729+2*ci;
 while(1,p=(t^2+19*f^2)/4;N=p+1-t;if(isprime(p)&&isprime(N),break);t+=2);
 for(si=1,4,
 c=[1,factor(f)[1,1],factor(f)[2,1],f][si];
 H=polclass(-19*c^2);js=polrootsmod(H,p);if(#js==0,error("missing stratum"));j=lift(js[1]);
 E=ellinit(ellfromj(Mod(j,p)));tw=0;
 for(tr=0,1,
 if(tr,E=ellinit(elltwist(E)));G=0;
 for(x=0,1000,ys=ellordinate(E,Mod(x,p));if(#ys,P=[Mod(x,p),ys[1]];if(ellmul(E,P,N)==[0],G=P;break)));
 if(G!=0,tw=tr;break));
 if(G==0,error("no chosen trace model"));
 if(ellmul(E,G,N)!=[0]||N<=4*sqrtint(p),error("order witness"));
 print("{\"class\":",ci,",\"stratum\":",si,",\"p\":",p,",\"N\":",N,",\"t\":",t,",\"fpi\":",f,",\"conductor\":",c,",\"D\":-19,\"j\":",j,",\"a\":",lift(E.a4),",\"b\":",lift(E.a6),",\"Gx\":",lift(G[1]),",\"Gy\":",lift(G[2]),",\"classpoly_degree\":",poldegree(H),",\"classpoly_root\":true,\"prime_proofs\":true,\"hasse_unique_order\":true}");
 ));
}
quit
