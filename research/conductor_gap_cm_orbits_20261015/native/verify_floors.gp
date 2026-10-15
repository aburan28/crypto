\\ Input is a GP vector per line, produced from exact decimal strings.
\\ [id,p,a,b,N,ell,D,crater_a,crater_b]
allocatemem(200000000);
lines=readstr("research/conductor_gap_cm_orbits_20261015/prepared/floors-for-verification.gpvec");
for(i=1,#lines, v=eval(lines[i]); id=v[1];p=v[2];a=v[3];b=v[4];N=v[5];l=v[6];D=v[7];t=p+1-N; if(!isprime(p)||!isprime(N)||!isprime(l),error("prime certificate failed")); if(t^2-4*p!=D*l^2||gcd(t,p)!=1,error("discriminant failed"));E=ellinit([a,b],p);if(E.disc==0,error("singular floor"));j=lift(E.j);jr=lift(ellinit([v[8],v[9]],p).j);if(j==jr||jr!=lift(Mod(if(D==-19,-884736,-884736000),p)),error("floor classification failed"));P=[0];for(x=0,1000,if(kronecker(lift(Mod(x,p)^3+a*x+b),p)==1,P=[Mod(x,p),sqrt(Mod(x,p)^3+a*x+b)];if(ellmul(E,P,N)==[0],break)));if(P==[0]||ellmul(E,P,N)!=[0]||N<=4*sqrtint(p),error("Hasse-unique order witness failed"));print("{\"id\":",id,",\"p\":\"",p,"\",\"N\":\"",N,"\",\"ell\":\"",l,"\",\"D\":",D,",\"a\":\"",a,"\",\"b\":\"",b,"\",\"Gx\":\"",lift(P[1]),"\",\"Gy\":\"",lift(P[2]),"\",\"j\":\"",j,"\",\"prime_proofs\":true,\"hasse_unique_point_order\":true,\"floor_conductor\":\"",l,"\"}"));
quit;
