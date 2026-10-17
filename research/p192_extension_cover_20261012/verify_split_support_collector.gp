\\ Independent PARI point-law replay of the 20-class split-divisor collector.
read("research/p192_extension_cover_20261012/type_i_split_support_input.gp");
E=ellinit([0,0,0,-3,B],P);
G=[Mod(GX,P),Mod(GY,P)];
chiF(q)={
  my(value);
  if(#q==1,return("I"));
  value=sum(i=1,#C,Mod(C[i],P)*q[1]^(i-1));
  if(value==0,return("0"));
  if(kronecker(lift(value),P)==1,return("+"),return("-"));
};
pointkey(q)={
  if(#q==1,return("O"));
  return(Str(lift(q[1]),":",lift(q[2])));
};
cellidx(label)={
  if(label=="++",return(1));
  if(label=="+-",return(2));
  if(label=="-+",return(3));
  if(label=="--",return(4));
  error("unexpected character cell");
};
verify()={
  my(pts=vector(#POINTS,i,[Mod(POINTS[i][1],P),Mod(POINTS[i][2],P)]),
     bearing=Map(),anchor=Map(),bc=vector(4),ac=vector(4),
     bearing_count=0,anchor_count=0,overlap=0,repeated=0,
     fresh_bearing=0,fresh_anchor=0,triples=0);
  if(#POINTS!=20 || #OLDINDEX!=20 || #KNOWN!=20 || #FRESH!=64,
    error("wrong input sizes"));
  if(ellmul(E,G,R)!=[0],error("wrong subgroup order"));
  for(i=1,#pts,if(!ellisoncurve(E,pts[i]),error("point off curve")));
  for(i=1,#pts,
    for(j=i,#pts,
      for(k=j,#pts,
        my(q=elladd(E,elladd(E,pts[i],pts[j]),pts[k]),row=vector(16),
           key,label,knownsum);
        triples++;
        if(q==[0],error("identity target"));
        for(t=1,3,
          my(idx=OLDINDEX[if(t==1,i,if(t==2,j,k))]);
          if(idx<16,row[idx+1]++);
          if(idx>=16 && idx<32,row[idx-15]--);
        );
        key=pointkey(q);
        label=Str(chiF(ellmul(E,q,(R-1)/2)),chiF(q));
        if(sum(t=1,16,abs(row[t]))>0,
          if(mapisdefined(bearing,key),repeated++);
          if(mapisdefined(anchor,key),overlap++);
          mapput(bearing,key,1);
          bearing_count++;
          bc[cellidx(label)]++;
        ,
          if(mapisdefined(anchor,key),repeated++);
          if(mapisdefined(bearing,key),overlap++);
          knownsum=(KNOWN[i]+KNOWN[j]+KNOWN[k])%R;
          if(ellmul(E,G,knownsum)!=q,error("known-scalar anchor mismatch"));
          mapput(anchor,key,1);
          anchor_count++;
          ac[cellidx(label)]++;
        );
      );
    );
  );
  for(i=1,#FRESH,
    my(pre=ellmul(E,G,FRESH[i]),q);
    if(Str(chiF(pre),chiF(ellmul(E,pre,2)))!="++",
      error("fresh query character mismatch"));
    q=ellmul(E,G,-2*FRESH[i]);
    fresh_bearing+=mapisdefined(bearing,pointkey(q));
    fresh_anchor+=mapisdefined(anchor,pointkey(q));
  );
  print("SPLIT_SUPPORT_GP_PASS triples=",triples,
        " bearing=",bearing_count," anchored=",anchor_count,
        " overlap=",overlap," repeated=",repeated,
        " bearing_cells=",bc," anchored_cells=",ac,
        " fresh_bearing=",fresh_bearing," fresh_anchor=",fresh_anchor);
};
verify();
quit();
