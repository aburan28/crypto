\\ Bounded two-branch search from frozen independent source curves.
FIELDS=Map();BUDGET=256;SEARCH_MS=90000;setrand(202610090988);
element(v,z)=sum(i=1,#v,eval(v[i])*z^(i-1));
cv(v)={v+=0*Z;vector(6,i,Str(lift(polcoef(v.pol,i-1))));};
key(v)=Str(v.pol);
setup(p,modulus)={
  my(k=Str(p,"/",modulus),saved);if(mapisdefined(FIELDS,k,&saved),[Z,BETA,EMBED,BACK,D]=saved;return);
  my(poly='z^6+sum(i=1,6,eval(modulus[i])*'z^(i-1)));if(!isprime(p)||!polisirreducible(Mod(poly,p)),error("invalid frozen field"));
  Z=ffgen(Mod(poly,p),'z);my(q=p^2,np=q^2+q+1);D=random(Z)^np;while(D==0||issquare(D),D=random(Z)^np);
  my(ext=ffextend(Z,'x^2-D,'w));BETA=ext[1];EMBED=ext[2];BACK=ffinvmap(EMBED);
  if(D^q!=D||BETA^q!=-BETA,error("nonsplit base-field setup failed"));
  mapput(FIELDS,k,[Z,BETA,EMBED,BACK,D]);
};
two_roots(E)=polrootsmod('x^3+E.a2*'x^2+E.a4*'x+E.a6);
plus_hit(r,q)={
  if(#r!=3,return(0));my(N=q^2+q+1,a=(r[2]-r[1])^N,b=(r[3]-r[1])^N,c=(r[3]-r[2])^N);return(a==b||c==-a||c==b);
};
minus_hit(E,r,q)={
  if(#r!=1,return([]));my(nm=q^2-q+1,ys=polrootsmod(256*('x-1)^3-E.j*('x-2)));
  for(i=1,#ys,
    my(dis=ys[i]^2-4);if(dis==0||issquare(dis),next);
    my(lambda=(ffmap(EMBED,ys[i])+sqrt(ffmap(EMBED,dis)))/2);
    if(lambda^nm!=1,next);
    my(gamma=lambda^lift(Mod(q+1,nm)^(-1)),alpha=ffmap(BACK,BETA*(gamma+1)/(gamma-1)));
    if(alpha==[]||alpha^q==alpha,error("quadratic reconstruction failed"));
    return([alpha,lambda]);
  );[]
};
verify_map(E,F,f)={for(i=1,3,my(P=random(E),image=ellisogenyapply(f,P));if(!ellisoncurve(F,image),error("isogeny image check failed")));};
verify_order(E,ord)={for(i=1,3,my(P=random(E));if(ellmul(E,P,ord)!=[0],error("point-order replay failed")));};
odd_kernels(E,ell,Q)={
  my(fac=factor(elldivpol(E,ell)),half=(ell-1)/2,result=List(),seen=Map());
  for(i=1,matsize(fac)[1],
    my(f=fac[i,1],deg=poldegree(f));if(half%deg,next);
    my(root,base=Z,embed=[],back=[],F=E);
    if(deg==1,root=-polcoef(f,0)/polcoef(f,1),
      my(ext=ffextend(Z,f,'u));root=ext[1];embed=ext[2];base=ffgen(root);F=ellinit(ffmap(embed,E[1..5]),base);
    );
    my(value=root^3+F.a4*root+F.a6,yy);
    if(issquare(value),yy=sqrt(value),
      my(ext=ffextend(base,'x^2-value,'v));yy=ext[1];my(m=ext[2]);root=ffmap(m,root);F=ellinit(ffmap(m,F[1..5]),ffgen(yy));embed=if(embed==[],m,ffcompomap(m,embed));
    );
    my(P=[root,yy]);if(!ellisoncurve(F,P)||ellmul(F,P,ell)!=[0],error("odd-kernel torsion point failed"));
    my(kernel=1+0*root);for(j=1,half,my(T=ellmul(F,P,j));if(#T==1,error("odd-kernel subgroup size failed"));kernel*=('x-T[1]));
    if(embed!=[],kernel=ffmap(ffinvmap(embed),kernel));if(kernel==[],next);
    if(poldegree(kernel)!=half,error("odd-kernel degree failed"));
    my(k=Str(kernel));if(!mapisdefined(seen,k),mapput(seen,k,1);listput(result,kernel));
  );Vec(result)
};
literal_plus(E,r,p,q,ord)={
  my(N=q^2+q+1,perm=[]);for(i=1,3,my(o=select(j->j!=i,[1,2,3]));if((r[o[1]]-r[i])^N==(r[o[2]]-r[i])^N,perm=[i,o[1],o[2]];break));
  if(!#perm,error("equal-norm root ordering failed"));my(origin=r[perm[1]],u=r[perm[2]]-origin,v=r[perm[3]]-origin,lambda=v/u,alpha=0*Z);
  for(j=0,5,my(c=Z^j,coef=1+0*Z);alpha=c;for(i=1,2,coef=coef^q/lambda;alpha+=coef*c^(q^i));if(alpha!=0,break));
  if(alpha==0||alpha^q!=lambda*alpha,error("cubic Hilbert-90 failed"));
  if(!issquare(u/alpha),alpha*=D);my(h=sqrt(u/alpha),W=ellinit([0,-alpha-alpha^q,0,alpha*alpha^q,0],Z),change=[h,origin,0,0]);
  if(ellchangecurve(E,change)[1..5]!=W[1..5]||ellcard(W)!=ord,error("cubic conversion or count failed"));
  for(i=1,3,my(P=random(E),P1=ellchangepoint(P,change));if(!ellisoncurve(W,P1)||ellchangepointinv(P1,change)!=P,error("cubic transport failed")));
  print("LITERAL_PLUS|alpha=",cv(alpha),"|h=",cv(h),"|origin=",cv(origin),"|short_a=",cv(W.a4-W.a2^2/3),"|short_b=",cv(2*W.a2^3/27-W.a2*W.a4/3),"|j=",cv(W.j));
};
literal_minus(E,hit,q,ord)={
  my(alpha=hit[1],f=('x^2-D)*('x-alpha)*('x-alpha^q),c3=subst(deriv(f),'x,alpha),c2=subst(deriv(deriv(f)),'x,alpha)/2,c1=subst(deriv(deriv(deriv(f))),'x,alpha)/6,W=ellinit([0,c2,0,c3*c1,c3^2],Z));
  my(a=W.a4-W.a2^2/3,b=W.a6-W.a2*W.a4/3+2*W.a2^3/27,twist=1+0*Z,h2=E.a6*a/(E.a4*b));
  if(!issquare(h2),twist=D;c3*=D;c2*=D;c1*=D;W=ellinit([0,c2,0,c3*c1,D*c3^2],Z);a=W.a4-W.a2^2/3;b=W.a6-W.a2*W.a4/3+2*W.a2^3/27;h2=E.a6*a/(E.a4*b));
  my(h=sqrt(h2),S=ellinit([a,b],Z));
  if(E.j!=W.j||h^4*a!=E.a4||h^6*b!=E.a6||ellcard(W)!=ord,error("quadratic conversion or twist count failed"));
  for(i=1,3,
    my(P=random(E));if(#P==1,next);my(X=P[1]/h^2-W.a2/3,Y=P[2]/h^3);if(X==0,next);
    my(xx=alpha+c3/X,yy=c3*Y/X^2);
    if(yy^2!=twist*(xx^2-D)*(xx-alpha)*(xx-alpha^q)||c3/(xx-alpha)+W.a2/3!=P[1]/h^2,error("quadratic point transport failed"));
  );
  print("LITERAL_MINUS|alpha=",cv(alpha),"|d=",cv(D),"|twist=",cv(twist),"|h=",cv(h),"|c3=",cv(c3),"|a2=",cv(W.a2),"|short_a=",cv(a),"|short_b=",cv(b),"|j=",cv(W.j));
};
construct_source(p,modulus,rv,ord,seed,icv1)={
  my(start=getwalltime());setup(p,modulus);my(setup_ms=getwalltime()-start,q=p^2,Q=p^6,r=vector(3,i,element(rv[i],Z)),a2=-vecsum(r),a4=r[1]*r[2]+r[1]*r[3]+r[2]*r[3],a6=-r[1]*r[2]*r[3],E=ellinit([a4-a2^2/3,a6-a2*a4/3+2*a2^3/27],Z));
  if(Q+1-ord==0||(Q+1-ord)%p==0,error("ordinary source required"));
  my(count_start=getwalltime());if(ellcard(E)!=ord,error("frozen source point count failed"));verify_order(E,ord);my(source_check_ms=getwalltime()-count_start);
  print("SOURCE|seed=",seed,"|p=",p,"|trace=",Q+1-ord,"|order=",ord,"|icv1=",icv1,"|setup_ms=",setup_ms,"|source_check_ms=",source_check_ms,"|modulus=",modulus,"|short_a=",cv(E.a4),"|short_b=",cv(E.a6));
  my(nodes=List(),seen=Map(),head=1,expanded=0,edges2=0,edges3=0,edges5=0,edges7=0,witness=0,branch="",hit=[],limit=0,search_start=getwalltime());
  listput(nodes,[E,0,0,0*Z]);mapput(seen,key(E.j),1);
  while(head<=#nodes,
    if(getwalltime()-search_start>=SEARCH_MS,limit=2;break);my(F=nodes[head][1],roots=two_roots(F));
    if(plus_hit(roots,q),witness=head;branch="PLUS";break);
    hit=minus_hit(F,roots,q);if(#hit,witness=head;branch="MINUS";break);
    expanded++;if(#nodes>=BUDGET,limit=1;break);
    for(di=1,4,
      my(deg=[2,3,5,7][di],kernels=if(deg==2,vector(#roots,i,'x-roots[i]),if(deg==3,my(xs=polrootsmod(elldivpol(F,3)));vector(#xs,i,'x-xs[i]),odd_kernels(F,deg,Q))));
      for(i=1,#kernels,
        my(iso=ellisogeny(F,kernels[i]),G=ellinit(iso[1],Z));if(deg==2,edges2++,if(deg==3,edges3++,if(deg==5,edges5++,edges7++)));
        verify_map(F,G,iso[2]);if(mapisdefined(seen,key(G.j)),next);mapput(seen,key(G.j),#nodes+1);listput(nodes,[G,head,deg,kernels[i]]);
        if(#nodes>=BUDGET,limit=1;break);
      );if(limit,break);
    );if(limit,break);head++;
  );
  my(search_ms=getwalltime()-search_start,verify_start=getwalltime(),route_len=0);
  if(witness,
    my(W=nodes[witness][1],route=List(),index=witness);if(ellcard(W)!=ord,error("endpoint point count failed"));verify_order(W,ord);
    print("ENDPOINT|branch=",branch,"|short_a=",cv(W.a4),"|short_b=",cv(W.a6),"|j=",cv(W.j),"|order=",ord);
    while(index>1,listput(route,index);index=nodes[index][2]);route_len=#route;
    forstep(i=#route,1,-1,
      my(child=route[i],node=nodes[child],parent=node[2],A=nodes[parent][1],B=node[1],iso=ellisogeny(A,node[4]));
      if(ellinit(iso[1],Z)[1..5]!=B[1..5],error("route reconstruction failed"));verify_map(A,B,iso[2]);if(p==7&&ellcard(B)!=ord,error("small route point count failed"));
      print("ROUTE|parent=",parent,"|child=",child,"|degree=",node[3],"|kernel_polynomial=",vector(poldegree(node[4])+1,i,cv(polcoef(node[4],i-1))),"|source_a=",cv(A.a4),"|source_b=",cv(A.a6),"|target_a=",cv(B.a4),"|target_b=",cv(B.a6));
    );
    if(branch=="PLUS",literal_plus(W,two_roots(W),p,q,ord),literal_minus(W,hit,q,ord));
  );
  print("CONSTRUCTION|seed=",seed,"|p=",p,"|status=",if(witness,"VERIFIED_WITNESS",if(limit==2,"TIME_CAP",if(limit==1,"VERTEX_CAP","RESTRICTED_COMPONENT_EXHAUSTED"))),"|branch=",branch,"|visited=",#nodes,"|tested=",if(witness||limit,head,head-1),"|expanded=",expanded,"|edges2=",edges2,"|edges3=",edges3,"|edges5=",edges5,"|edges7=",edges7,"|route_length=",route_len,"|search_ms=",search_ms,"|verification_ms=",getwalltime()-verify_start);
};
read(getenv("ISO1_INPUTS"));
