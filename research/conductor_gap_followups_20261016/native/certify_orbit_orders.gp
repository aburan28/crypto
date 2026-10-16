{
 rows=readstr("research/conductor_gap_followups_20261016/prepared/orbit-order.gpvec");
 for(k=1,#rows,v=eval(rows[k]);id=v[1];r=v[2];lam=v[3];if(!isprime(r),error("nonprime subgroup"));ord=znorder(Mod(lam,r));fa=factor(ord);if(lift(Mod(lam,r)^ord)!=1,error("order witness"));for(i=1,matsize(fa)[1],if(!isprime(fa[i,1]),error("factor prime proof"));if(lift(Mod(lam,r)^(ord/fa[i,1]))==1,error("minimality failed")));signed=ord;if(ord%2==0&&lift(Mod(lam,r)^(ord/2))==r-1,signed=ord/2);print("{\"pair\":",id,",\"r\":",r,",\"lambda\":",lam,",\"multiplicative_order\":",ord,",\"signed_period\":",signed,",\"order_prime_factor_witnesses\":true,\"ok\":true}"));
}
quit
