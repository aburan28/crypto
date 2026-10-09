\\ Independent finite-field controls for the generalized norm-one theorem.
\\ Each case is [characteristic, base-field degree, odd extension degree].
setrand(20261008);
cases=[[5,1,3],[5,1,5],[5,1,7],[3,2,3],[3,2,5],[13,1,3],[17,1,3],[5,2,3],[7,2,3],[11,2,3],[13,2,3]];
for(ci=1,#cases,p=cases[ci][1];a=cases[ci][2];d=cases[ci][3];q=p^a;Q=q^d;N=(Q-1)/(q-1);u=lift(1/Mod(4,N));z=ffgen([p,a*d],'z);ordinary=0;for(i=1,100,lam=random(z)^(q-1);while(lam==0||lam==1,lam=random(z)^(q-1));mu=lam^u;if(mu^4!=lam,error("fourth-root check failed"));s=mu^2;E=ellinit([0,-(1+lam),0,lam,0],z);Ep=ellinit([0,2*(1+lam),0,(1-lam)^2,0],z);ne=ellcard(E);np=ellcard(Ep);if(ne!=np,error("isogeny count mismatch"));if(np%16,error("target cardinality not divisible by 16"));gr=ellgroup(Ep);if(#gr!=2||gr[2]%4,error("target lacks full rational 4-torsion"));t=Q+1-ne;if(t%p,ordinary++;D=t^2-4*Q;DK=coredisc(D);f=sqrtint(D/DK);if(f^2*DK!=D||f%4,error("ordinary conductor check failed"))));print("q=",q,",n=",d,",samples=100,ordinary=",ordinary,",fourth_root_failures=0,isogeny_failures=0,full4_failures=0,conductor_failures=0"));
\\ Check the trace/conductor equivalence on every ordinary Hasse trace
\\ for several nonsquare and square cardinalities with Q=1 (mod 4).
Qs=[[5,5],[13,13],[17,17],[5,25],[3,81],[5,125],[13,169],[17,289],[3,729],[5,3125]];
for(ci=1,#Qs,p=Qs[ci][1];Q=Qs[ci][2];m=sqrtint(4*Q);checked=0;for(t=-m,m,if(t%p&&t^2<4*Q,D=t^2-4*Q;DK=coredisc(D);f=sqrtint(D/DK);left=(f%4==0);right=(lift(Mod(t-(Q+1),16))==0||lift(Mod(t+(Q+1),16))==0);if(left!=right,error("trace/conductor equivalence failed"));checked++));print("Q=",Q,",ordinary_hasse_traces=",checked,",equivalence_failures=0"));
