\\ Loaded after construct.gp defines the direct test and literal conversion.
run_controls()={
  my(p=7,q=p^2,Q=p^6,modulus=["5","6","5","3","5","5"]);setup(p,modulus);setrand(202610090999);
  my(quad=0,cubic=0);
  for(i=1,64,
    my(alpha=random(Z));while(alpha^q==alpha,alpha=random(Z));my(f=('x^2-D)*('x-alpha)*('x-alpha^q),c3=subst(deriv(f),'x,alpha),c2=subst(deriv(deriv(f)),'x,alpha)/2,c1=subst(deriv(deriv(deriv(f))),'x,alpha)/6,W=ellinit([0,c2,0,c3*c1,c3^2],Z),E=ellinit([W.a4-W.a2^2/3,W.a6-W.a2*W.a4/3+2*W.a2^3/27],Z));
    my(r=two_roots(E),hit=minus_hit(E,r,q));if(#r!=1||!#hit||plus_hit(r,q),error("quadratic direct criterion control failed"));
    literal_minus(E,hit,q,ellcard(E));quad++;
    my(lambda=random(Z)^(q-1));while(lambda==0||lambda==1,lambda=random(Z)^(q-1));W=ellinit([0,-1-lambda,0,lambda,0],Z);E=ellinit([W.a4-W.a2^2/3,2*W.a2^3/27-W.a2*W.a4/3],Z);r=two_roots(E);
    if(#r!=3||!plus_hit(r,q)||#minus_hit(E,r,q),error("cubic direct criterion control failed"));literal_plus(E,r,p,q,ellcard(E));cubic++;
  );
  print("DIRECT_COMPLETE|quadratic_controls=",quad,"|cubic_controls=",cubic);
};
run_controls();
