\\ Independent symbolic replay and exhaustive direct-locus controls.
print("duplication_difference=",(3*x^2+2*a*x+b)^2-4*(x^3+a*x^2+b*x)*(a+2*x)-(x^2-b)^2);
forstep(n=3,7,2,print("norm_fiber_degree_",n,"=",(q-1)*sum(i=0,n-1,q^i)-(q^n-1)));
norm_locus(p,n)={
  my(k=2*n,Q=p^k,q=p^2,N=(Q-1)/(q-1),z=ffgen([p,k],'z),one=0,second=0,third=0,union=0,pair=0);
  for(i=0,Q-1,
    my(d=digits(i,p),v=0*z);for(j=1,#d,v=v*z+d[j]);
    if(v==0||v==1,next);
    my(A=v^N==1,B=(v-1)^N==-1,C=((v-1)/v)^N==1);
    one+=A;second+=B;third+=C;union+=(A||B||C);pair+=(A&&B);
    if(A&&B&&C,error("unexpected triple intersection"));
  );
  if(one!=N-1||second!=N-1||third!=N-1||union>3*(N-1),error("norm-fiber count failed"));
  print("direct_locus|p=",p,"|n=",n,"|Q=",Q,"|fiber_events=",[one,second,third],"|union=",union,"|union_bound=",3*(N-1),"|first_pair=",pair);
};
norm_locus(5,3);norm_locus(7,3);norm_locus(3,5);
