\\ Exact finite-field checks. A COMPLETE marker and empty stderr are required.
cv(v)=vector(DEGREE,i,Str(lift(polcoef(v.pol,i-1))));
root_curve(r)=ellinit([0,-vecsum(r),0,r[1]*r[2]+r[1]*r[3]+r[2]*r[3],-r[1]*r[2]*r[3]],Z);
emit_curve(role,E,r)={
  print("CURVE|",role,"|",cv(E.a4-E.a2^2/3),"|",cv(E.a6-E.a2*E.a4/3+2*E.a2^3/27),"|",cv(E.j),"|",vector(3,i,cv(r[i])));
};
weak(r)={
  my(a=(r[2]-r[1])^NORMEXP,b=(r[3]-r[1])^NORMEXP,c=(r[3]-r[2])^NORMEXP);
  return(a==b||c==-a||c==b);
};
full4(r)={return(issquare(r[2]-r[1])&&issquare(r[3]-r[1])&&issquare(r[3]-r[2]));};
verify_roots(E,r)={
  if(#r!=3||#Set(r)!=3,error("2-torsion roots distinctness failed"));
  for(i=1,3,if(!ellisoncurve(E,[r[i],0*Z])||ellmul(E,[r[i],0*Z],2)!=[0],error("2-torsion root check failed")));
};
verify_points(E,ord)={for(i=1,3,my(P=random(E));if(!ellisoncurve(E,P)||ellmul(E,P,ord)!=[0],error("order replay failed")));};
verify_halves(E,r)={
  for(i=1,2,
    my(other=select(j->j!=i,[1,2,3]),a=sqrt(r[i]-r[other[1]]),b=sqrt(r[i]-r[other[2]]),P=[r[i]+a*b,a*b*(a+b)]);
    if(!ellisoncurve(E,P)||ellmul(E,P,2)!=[r[i],0*Z]||ellmul(E,P,4)!=[0],error("full-4 halving failed"));
  );
};
verify_map(E,F,f)={for(i=1,3,my(P=random(E),P1=ellisogenyapply(f,P));if(!ellisoncurve(F,P1),error("isogeny image failed")));};
neighbor2(E,r,i)={
  my(other=select(j->j!=i,[1,2,3]),a=-(r[other[1]]+r[other[2]]-2*r[i]),b=(r[other[1]]-r[i])*(r[other[2]]-r[i]));
  if(!issquare(b),return([]));
  my(iso=ellisogeny(E,'x-r[i]),F=ellinit(iso[1],Z),s=sqrt(b),rf=[r[i]-a,r[i]+2*s,r[i]-2*s]);
  verify_roots(F,rf);verify_map(E,F,iso[2]);
  return([F,rf,iso[2],r[i]]);
};
neighbor3(E,r,kernel)={
  if(subst(elldivpol(E,3),'x,kernel)!=0,error("degree-3 kernel not on division polynomial"));
  my(iso=ellisogeny(E,'x-kernel),F=ellinit(iso[1],Z),rf=vector(3,i,ellisogenyapply(iso[2],[r[i],0*Z])[1]));
  verify_roots(F,rf);verify_map(E,F,iso[2]);
  return([F,rf,iso[2],kernel]);
};
conductor_depth(t)={
  my(d=(t/2)^2-Q,v=valuation(d,2),u=d/2^v);
  if(v%2,return((v-1)/2));
  return(v/2+if(u%4==1,1,0));
};
full4_neighbor(E,r,ord)={
  if(full4(r),verify_halves(E,r);return([1,0,0*Z,E,r]));
  for(i=1,3,my(nb=neighbor2(E,r,i));if(#nb&&full4(nb[2]),verify_halves(nb[1],nb[2]);return([1,2,r[i],nb[1],nb[2]])));
  return([0,0,0*Z,E,r]);
};
search_class(E,r,ord)={
  my(start=getwalltime(),nodes=List(),seen=Map(),head=1,edges2=0,edges3=0,discarded2=0,expanded=0,witness=0,limit=0);
  listput(nodes,[E,r,0,0,0*Z]);mapput(seen,cv(E.j),1);
  while(head<=#nodes,
    if(getwalltime()-start>=SEARCH_MS,limit=2;break);
    my(node=nodes[head],F=node[1],rf=node[2]);
    if(weak(rf),witness=head;break);
    if(#nodes>=BUDGET,limit=1;break);
    expanded++;
    for(i=1,3,
      my(nb=neighbor2(F,rf,i));
      if(!#nb,discarded2++;next);
      edges2++;
      if(mapisdefined(seen,cv(nb[1].j)),next);
      mapput(seen,cv(nb[1].j),#nodes+1);listput(nodes,[nb[1],nb[2],head,2,nb[4]]);
      if(weak(nb[2]),witness=#nodes;break);
      if(#nodes>=BUDGET,limit=1;break);
    );
    if(witness||limit,break);
    my(kernels=polrootsmod(elldivpol(F,3)));
    for(i=1,#kernels,
      my(nb=neighbor3(F,rf,kernels[i]));edges3++;
      if(mapisdefined(seen,cv(nb[1].j)),next);
      mapput(seen,cv(nb[1].j),#nodes+1);listput(nodes,[nb[1],nb[2],head,3,nb[4]]);
      if(weak(nb[2]),witness=#nodes;break);
      if(#nodes>=BUDGET,limit=1;break);
    );
    if(witness||limit,break);
    head++;
  );
  my(status=if(witness,"WITNESS",if(limit==2,"TIME_CAP",if(limit==1,"VERTEX_CAP","COMPONENT_EXHAUSTED"))));
  if(witness,
    my(W=nodes[witness][1],rw=nodes[witness][2],count=ellcard(W));
    if(count!=ord||!weak(rw),error("witness count or norm failed"));
    verify_points(W,ord);emit_curve("witness",W,rw);
    my(route=List(),idx=witness);
    while(idx>1,listput(route,idx);idx=nodes[idx][3]);
    forstep(i=#route,1,-1,
      my(child=route[i],parent=nodes[child][3],node=nodes[child]);
      print("ROUTE|",parent,"|",child,"|",node[4],"|",cv(node[5]));
      emit_curve(Str("route",child),node[1],node[2]);
      if(SMALL_AUDIT&&ellcard(node[1])!=ord,error("small route count failed"));
    );
  );
  print("SEARCH|",status,"|",#nodes,"|",expanded,"|",edges2,"|",edges3,"|",discarded2,"|",getwalltime()-start);
  return(witness!=0);
};
run_sample()={
  my(p=eval(getenv("ISO1_P")),n=eval(getenv("ISO1_N")),seed=eval(getenv("ISO1_SEED")),mode=getenv("ISO1_MODE"));
  DEGREE=2*n;Q=p^DEGREE;NORMEXP=(Q-1)/(p^2-1);BUDGET=eval(getenv("ISO1_BUDGET"));SEARCH_MS=90000;SMALL_AUDIT=p==7;
  if(!isprime(p)||n<3||n%2==0,error("invalid field"));
  setrand(2026100900);Z=ffgen([p,DEGREE],'z);setrand(seed);
  print("FIELD|",p,"|",n,"|",Q,"|",vector(DEGREE,i,Str(lift(polcoef(Z.mod,i-1)))));
  my(r,E,control_attempts=0);
  if(mode=="control",
    while(1,
      control_attempts++;my(lam=random(Z)^(p^2-1));if(lam==0||lam==1,next);
      my(wr=[0*Z,1+0*Z,lam],WE=root_curve(wr),nb=neighbor2(WE,wr,1));
      if(!#nb,error("norm-one neighbor missing"));
      if(weak(nb[2]),next);
      if(control_attempts>100,error("control off-locus generation cap"));
      E=nb[1];r=nb[2];emit_curve("control_known_weak",WE,wr);break;
    );
  ,
    my(u=random(Z),v=random(Z));while(u==0||v==0||u==v,u=random(Z);v=random(Z));
    r=[0*Z,u,v];E=root_curve(r);
  );
  verify_roots(E,r);emit_curve("source",E,r);
  my(start=getwalltime(),ord=ellcard(E),t=Q+1-ord,count_ms=getwalltime()-start,ordinary=t%p!=0,admitted=(t-(Q+1))%16==0||(t+(Q+1))%16==0,direct=weak(r),depth=if(ordinary,conductor_depth(t),-1));
  if(t^2>4*Q||ord%4,error("Hasse or full-2 count failed"));
  verify_points(E,ord);
  if(ordinary&&(admitted!=(depth>=2)||direct&&!admitted),error("conductor condition failed"));
  my(epsilon=if(ord%16==0,1,-1),R=r,F=E,selected_ord=ord);
  if(admitted&&epsilon==-1,
    my(ns=random(Z));while(ns==0||issquare(ns),ns=random(Z));
    R=ns*r;F=root_curve(R);selected_ord=2*Q+2-ord;verify_roots(F,R);verify_points(F,selected_ord);emit_curve("selected_twist",F,R);
  );
  print("SAMPLE|",seed,"|",mode,"|",t,"|",ord,"|",ordinary,"|",admitted,"|",direct,"|",depth,"|",full4(r),"|",epsilon,"|",count_ms,"|",control_attempts);
  if(ordinary&&admitted,
    my(four=full4_neighbor(F,R,selected_ord));
    print("TORSION|",four[1],"|",four[2],"|",cv(four[3]));
    if(!four[1],error("admitted class has no immediate full-4 neighbor"));
    emit_curve("full4",four[4],four[5]);
    search_class(F,R,selected_ord);
  ,print("SEARCH|NOT_ADMITTED|0|0|0|0|0|0"));
  print("COMPLETE|",p,"|",n,"|",seed,"|",mode);
};
run_sample();
