\\ Candidate search and rigorous prime/kernel certificates. No ECDLP scalar.
\\ Family q=(|D|+1)/4, p=q*l^2-l+1, N=q*l^2-2*l+4.
\\ Twist cardinality q*l^2 and pi == 1 mod l give rational l-torsion.
allocatemem(1000000000);
make(ell0,D,id)={
 my(q=(D+1)/4,l=nextprime(ell0),p,N,t,j,k,E,T,P,K,dd,aa,bb,done=0);
 while(!done,
  if(l%4==3,p=q*l^2-l+1;N=q*l^2-2*l+4;if(isprime(p)&&isprime(N),done=1));
  if(!done,l=nextprime(l+1))
 );
 t=l-2;j=if(D==19,-884736,-884736000);k=Mod(j,p)/(1728-j);aa=lift(3*k);bb=lift(2*k);E=ellinit([aa,bb],p);
 P=[];for(x=0,100,if(kronecker(lift(Mod(x,p)^3+aa*x+bb),p)==1,P=[Mod(x,p),sqrt(Mod(x,p)^3+aa*x+bb)];break));
 if(ellmul(E,P,N)!=[0],dd=2;while(kronecker(dd,p)!=-1,dd++);aa=lift(Mod(aa*dd^2,p));bb=lift(Mod(bb*dd^3,p));E=ellinit([aa,bb],p));
 for(x=0,100,if(kronecker(lift(Mod(x,p)^3+aa*x+bb),p)==1,P=[Mod(x,p),sqrt(Mod(x,p)^3+aa*x+bb)];if(ellmul(E,P,N)==[0],break)));
 if(P==[0]||ellmul(E,P,N)!=[0],error("crater certificate failed"));
 dd=2;while(kronecker(dd,p)!=-1,dd++);my(at=lift(Mod(aa/dd^2,p)),bt=lift(Mod(bb/dd^3,p)));T=ellinit([at,bt],p);
 for(x=1,100,if(kronecker(lift(Mod(x,p)^3+at*x+bt),p)==1,K=ellmul(T,[Mod(x,p),sqrt(Mod(x,p)^3+at*x+bt)],q);if(K!=[0]&&ellmul(T,K,l)==[0],break)));
 if(K==[0]||ellmul(T,K,l)!=[0],error("kernel certificate failed"));
 print(id," ",l," ",p," ",at," ",bt," ",lift(K[1])," ",lift(K[2])," ",dd," ",aa," ",bb," -1 -1");
 write("research/conductor_gap_cm_orbits_20261015/prepared/d43-control-crater.jsonl",Str("{\"id\":",id,",\"D\":-",D,",\"ell\":\"",l,"\",\"p\":\"",p,"\",\"N\":\"",N,"\",\"t\":\"",t,"\",\"a\":\"",aa,"\",\"b\":\"",bb,"\",\"Gx\":\"",lift(P[1]),"\",\"Gy\":\"",lift(P[2]),"\",\"prime_proofs\":true,\"crater_generator_order\":true,\"twist_kernel_order\":true}"));
};
make(263,43,9000);
quit;
