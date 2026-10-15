ell=129149576694539;t=2105;p=79228162511798403263235114731;N=p+1-t;
if(!isprime(ell)||!isprime(p)||!isprime(N)||t^2-4*p!=-19*ell^2,error("candidate certificates failed"));
s=Mod(t,ell)/2;d=znorder(s);dt=znorder(-s);dx=d/gcd(d,2);dtx=dt/gcd(dt,2);
if(s^d!=1||(-s)^dt!=1,error("order witness failed"));
for(i=1,matsize(factor(d))[1],if(s^(d/factor(d)[i,1])==1,error("order not minimal")));
for(i=1,matsize(factor(dt))[1],if((-s)^(dt/factor(dt)[i,1])==1,error("twist order not minimal")));
print("{\"ell\":\"",ell,"\",\"p\":\"",p,"\",\"N\":\"",N,"\",\"t\":",t,",\"frobenius_scalar_mod_ell\":\"",lift(s),"\",\"point_extension_degree\":\"",d,"\",\"twist_point_extension_degree\":\"",dt,"\",\"x_extension_degree\":\"",dx,"\",\"twist_x_extension_degree\":\"",dtx,"\",\"primes_and_minimal_orders_checked\":true,\"dense_Fp_basis_bytes_per_coordinate_at_96_bits\":\"",12*min(d,dt),"\",\"floor_constructed\":false}");
quit;
