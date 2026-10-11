\\ Loaded after construct_general.gp; controls are separate from the population.
general_control(p,degree,modulus,seed)={
  setup(p,degree,modulus);setrand(seed);my(q=p^2,Q=p^degree,nm=(Q+1)/(q+1),w=ffgen(BETA),lambda=random(w)^((Q^2-1)/nm));while(lambda==0||lambda==1,lambda=random(w)^((Q^2-1)/nm));
  my(alpha=general_reconstruct(lambda,q),f=('x^2-D)*('x-alpha)*('x-alpha^q),c3=subst(deriv(f),'x,alpha),c2=subst(deriv(deriv(f)),'x,alpha)/2,c1=subst(deriv(deriv(deriv(f))),'x,alpha)/6,W=ellinit([0,c2,0,c3*c1,c3^2],Z),E=ellinit([W.a4-W.a2^2/3,W.a6-W.a2*W.a4/3+2*W.a2^3/27],Z),ord=ellcard(E));
  if(ord%4||E.j!=ffmap(BACK,256*(lambda^2-lambda+1)^3/(lambda^2*(lambda-1)^2)),error("generalized quadratic invariant or count failed"));my(r=two_roots(E),hit=minus_hit(E,r,q));if(#r!=1||!#hit,error("generalized quadratic endpoint test failed"));
  print("GENERAL_CONTROL|p=",p,"|degree=",degree,"|cover_genus=",1+2^(degree/2-2)*(degree/2-2),"|seed=",seed,"|trace=",Q+1-ord,"|order=",ord,"|modulus=",modulus,"|short_a=",cv(E.a4),"|short_b=",cv(E.a6),"|j=",cv(E.j),"|norm_section_gcd=",gcd(q+1,nm));
  verify_order(E,ord);literal_minus(E,hit,q,ord);
};
read(getenv("ISO1_GENERAL_CONTROL_INPUTS"));
