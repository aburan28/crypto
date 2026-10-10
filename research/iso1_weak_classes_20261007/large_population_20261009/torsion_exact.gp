\\ Independent exhaustive controls for the one-edge full-4 theorem.
square_differences(r)=#r==3&&issquare(r[2]-r[1])&&issquare(r[3]-r[1])&&issquare(r[3]-r[2]);
audit_field(p,k)={
  my(z=ffgen([p,k],'z),Q=p^k,elements=vector(Q,i,sum(j=0,k-1,((i-1)\p^j)%p*z^j)),pairs=0,admitted=0,already=0,edge=0);
  for(i=2,Q,for(j=2,Q,
    if(i==j,next);my(u=elements[i],v=elements[j],E=ellinit([0,-u-v,0,u*v,0],z),r=[0*z,u,v],order=ellcard(E),found=square_differences(r));pairs++;
    if(found,already++);
    for(h=1,3,
      my(F=ellinit(ellisogeny(E,'x-r[h],1),z),roots=polrootsmod('x^3+F.a2*'x^2+F.a4*'x+F.a6));
      if(ellcard(F)!=order,error("independent quotient count mismatch"));
      if(square_differences(roots),found=1);
    );
    if(found!=(order%16==0),error("one-edge torsion equivalence failed"));
    admitted+=(order%16==0);edge+=(found&&!square_differences(r));
  ));
  print("ONE_EDGE_EXACT|p=",p,"|degree=",k,"|Q=",Q,"|pairs=",pairs,"|order_div16=",admitted,"|already_full4=",already,"|requires_edge=",edge);
};
audit_field(5,1);audit_field(3,2);audit_field(13,1);audit_field(17,1);
audit_field(5,2);audit_field(29,1);audit_field(37,1);audit_field(7,2);
