\\ Exact line-section residual screen for the frozen Type I support divisors.
read("research/p192_extension_cover_20261012/type_i_fixed_origin_native_input.gp");
read("research/p192_extension_cover_20261012/type_i_divisor_support_input.gp");
read("research/p192_extension_cover_20261012/type_i_nonspecial_support_input.gp");
v=varhigher("v",x);
curve=sum(i=1,#terms,Mod(terms[i][3],p)*x^terms[i][1]*v^terms[i][2]);
check()={
  my(points=List(),newpts=List(),pairids=List(),counts=vector(4,i,0),ramified=0,splitcounts=vector(2,i,0));
  if(P!=p,error("field mismatch"));
  for(i=1,#X,
    my(roots=polrootsmod(X[i],p));
    for(j=1,#roots,
      my(xx=lift(roots[j]),yy=lift(Mod(subst(Y[i],z,lift(roots[j])),p)));
      if(subst(subst(curve,x,Mod(xx,p)),v,Mod(yy,p))!=0,error("point off quartic"));
      listput(points,[i,xx,yy]);
    );
  );
  points=Vec(points);
  if(#points!=141,error("wrong point count"));
  for(i=1,#points,write("research/p192_extension_cover_20261012/type_i_line_base_points.tsv",
      points[i][1]," ",points[i][2]," ",points[i][3]));
  for(i=1,#points-1,
    for(j=i+1,#points,
      my(x1=points[i][2],y1=points[i][3],x2=points[j][2],y2=points[j][3]);
      my(m=Mod(y2-y1,p)/Mod(x2-x1,p),b=Mod(y1,p)-m*Mod(x1,p));
      my(h=(x-Mod(x1,p))*(x-Mod(x2,p)),f=subst(curve,v,m*x+b),qr=divrem(f,h));
      if(qr[2]!=0,error("pair not on line section"));
      my(res=qr[1],degree=poldegree(res,x));
      if(degree<0 || degree>2,error("bad residual degree"));
      counts[degree+1]++;
      if(degree==2,
        my(a=polcoeff(res,2),bb=polcoeff(res,1),c=polcoeff(res,0));
        my(disc=bb^2-4*a*c,chi=kronecker(lift(disc),p));
        if(chi==0,ramified++);
        if(chi==1,
          my(roots=polrootsmod(res,p));
          if(#roots!=2,error("wrong residual rational roots"));
          for(k=1,#roots,
            my(rx=lift(roots[k]),ry=lift(m*roots[k]+b));
            listput(newpts,[rx,ry]);
            listput(pairids,[i,j]);
          );
          splitcounts[2]++;
        ,splitcounts[1]++);
      );
    );
  );
  newpts=Vec(newpts);
  pairids=Vec(pairids);
  if(#newpts!=2*splitcounts[2],error("wrong new point count"));
  if(#pairids!=#newpts,error("pair identity count"));
  my(all=vecsort(concat(vector(#points,i,[points[i][2],points[i][3]]),newpts)));
  for(i=2,#all,if(all[i]==all[i-1],error(Str("reused rational point ",i))));
  print("LINE_RESIDUAL_GP_PASS points=",#points," pairs=9870 degree_counts=",counts,
    " repeated_quadratic_roots=",ramified," nonsplit_split=",splitcounts,
    " new_rational_points=",#newpts," all_points_distinct=",#all);
  print("FIRST_NEW_RATIONAL_POINT=",newpts[1]);
  for(i=1,#newpts,
    write("research/p192_extension_cover_20261012/type_i_line_new_points.tsv",newpts[i][1]," ",newpts[i][2]);
    write("research/p192_extension_cover_20261012/type_i_line_new_point_pairs.tsv",pairids[i][1]," ",pairids[i][2]);
  );
};
check();
quit();
