default(parisizemax,2000000000);
{
 rows=readstr("research/conductor_gap_followups_20261016/prepared/dose-cert.gpvec");
 for(k=1,#rows,
 v=eval(rows[k]);cl=v[1];s=v[2];p=v[3];N=v[4];t=v[5];f=v[6];c=v[7];D=v[8];
 if(!isprime(p)||!isprime(N)||gcd(t,p)!=1||t^2-4*p!=D*f^2,error("prime or discriminant failed"));
 E=ellinit([v[9],v[10]],p);P=[Mod(v[11],p),Mod(v[12],p)];
 if(E.disc==0||!ellisoncurve(E,P)||ellmul(E,P,N)!=[0]||N^2<=16*p||t^2>4*p,error("group witness failed"));
 ds=divisors(f);matches=List();
 for(i=1,#ds,H=polclass(D*ds[i]^2);if(subst(H,variable(H),E.j)==0,listput(matches,ds[i])));
 if(Vec(matches)!=[c],error("ambiguous conductor"));
 print("{\"class\":",cl,",\"stratum\":",s,",\"p\":",p,",\"N\":",N,",\"conductor\":",c,",\"all_divisor_classpolys_checked\":true,\"prime_proofs\":true,\"hasse_unique_order\":true,\"ok\":true}");
 );
}
quit
