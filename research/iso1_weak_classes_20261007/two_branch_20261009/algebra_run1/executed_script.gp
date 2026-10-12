\\ Symbolic verification is separate from the finite-field census.
jleg(v)=256*(v^2-v+1)^3/(v^2*(v-1)^2);
run_algebra()={
  my(u='u,v='v,b='b,alpha=b*(u+1)/(u-1),alphaq=-b*(v+1)/(v-1),f=('x^2-b^2)*('x-alpha)*('x-alphaq));
  my(c3=subst(deriv(f),'x,alpha),c2=subst(deriv(deriv(f)),'x,alpha)/2,c1=subst(deriv(deriv(deriv(f))),'x,alpha)/6,b2=4*c2,b4=2*c3*c1,b6=4*c3^2,c4=b2^2-24*b4,c6=-b2^3+36*b2*b4-216*b6,j=1728*c4^3/(c4^3-c6^2));
  if(j!=jleg(u*v),error("symbolic quadratic cross-ratio invariant failed"));
  my(l='l,Y=l+1/l);if(jleg(l)!=256*(Y-1)^3/(Y-2),error("symbolic reciprocal invariant failed"));
  if(('q+1)*('q^2-'q+1)!='q^3+1,error("alternating norm cardinality failed"));
  print("SYMBOLIC|quartic_j_equals_legendre_cross_ratio=1|reciprocal_j_formula=1|alternating_norm_factor=1");
  forprime(p=5,199,
    my(N=p^4-p^2+1);if(gcd(p^2+1,N)!=1||gcd(N,p^6-1)!=1,error("quadratic torus section or field failed"));
    fordiv(12,d,if(d<12&&gcd(N,p^d-1)!=1,error("quadratic orbit length failed")));
  );
  print("ALGEBRA_COMPLETE|symbolic_identities=3|orbit_primes_through=199");
};
run_algebra();
