\\ Exact second-generation line-incidence screen on 141 old and first 64 new points.
read("research/p192_extension_cover_20261012/type_i_fixed_origin_native_input.gp");
read("research/p192_extension_cover_20261012/type_i_divisor_support_input.gp");
read("research/p192_extension_cover_20261012/type_i_nonspecial_support_input.gp");
read("research/p192_extension_cover_20261012/type_i_line_second_gen_input.gp");
v=varhigher("v",x);
curve=sum(i=1,#terms,Mod(terms[i][3],p)*x^terms[i][1]*v^terms[i][2]);

check()={
  my(selected=concat(line_old,vector(64,k,line_new[k])));
  my(existing=vecsort(concat(line_old,line_new)));
  my(quads=vecsort(line_original_quads));
  my(total=0,split=0,nonsplit=0,ramified=0,degother=0,oneexisting=0,closed=0,original=0,novel=0);
  my(witnesses=List());
  if(#line_old!=141 || #line_new!=9784 || #selected!=205 || #quads!=32,error("wrong input dimensions"));
  if(#existing!=#Set(existing),error("source points not distinct"));
  for(i=142,205,
    for(j=1,i-1,
      my(x1=selected[i][1],y1=selected[i][2],x2=selected[j][1],y2=selected[j][2]);
      if(x1==x2,error("vertical selected pair"));
      total++;
      my(m=Mod(y2-y1,p)/Mod(x2-x1,p),b=Mod(y1,p)-m*Mod(x1,p));
      my(f=subst(curve,v,m*x+b),h=(x-Mod(x1,p))*(x-Mod(x2,p)),qr=divrem(f,h));
      if(qr[2]!=0,error("selected pair off line section"));
      my(res=qr[1]);
      if(poldegree(res,x)!=2,degother++;next());
      my(a=polcoeff(res,2),bb=polcoeff(res,1),c=polcoeff(res,0));
      my(chi=kronecker(lift(bb^2-4*a*c),p));
      if(chi==0,ramified++;next());
      if(chi!=1,nonsplit++;next());
      split++;
      my(roots=polrootsmod(res,p));
      if(#roots!=2,error("wrong split residual roots"));
      my(r1=[lift(roots[1]),lift(m*roots[1]+b)]);
      my(r2=[lift(roots[2]),lift(m*roots[2]+b)]);
      my(n=(setsearch(existing,r1)!=0)+(setsearch(existing,r2)!=0));
      if(n==1,oneexisting++);
      if(n==2,
        closed++;
        my(q=vecsort([selected[i],selected[j],r1,r2]));
        if(setsearch(quads,q),original++,novel++;listput(witnesses,[i,j,q]));
      );
    );
  );
  if(total!=11040,error("wrong pair count"));
  print("LINE_SECOND_GEN_PASS pairs=",total," split=",split," nonsplit=",nonsplit," ramified=",ramified,
    " degree_other=",degother," one_existing=",oneexisting," closed=",closed,
    " original=",original," novel=",novel);
  print("LINE_SECOND_GEN_NOVEL_WITNESSES=",Vec(witnesses));
};
check();
quit();
