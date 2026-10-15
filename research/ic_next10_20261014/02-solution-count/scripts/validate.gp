\\ Validation arm: for each small base and target, (1) count grid points of V^m on S_{m+1} = 0 by direct
\\ evaluation of the summation polynomial (no group law), and (2) write msolve input files for the
\\ combined system {S_{m+1}(x_1..x_m, x_R), L(x_1), .., L(x_m)}, plus a random non-split L of the same
\\ degree and the partial system without L(x_m).
\\ Called with globals: P, A, B (curve), DIR (validation dir), TMAX2, TMAX3, NT (targets), MS (msolve tmax3).
S3(u,v,w) = (u-v)^2*w^2 - 2*((u+v)*(u*v+A)+2*B)*w + (u*v-A)^2 - 4*B*(u+v);
\\ expanded string of polynomial f (integer coeffs mod P) in variables vars, for msolve
monorec(g, vars, i, pre, ~L) = {
  if(i > #vars, my(c = lift(Mod(g, P))); if(c, listput(L, concat(Str(c), pre))); return);
  my(d = poldegree(g, vars[i]));
  for(e = 0, d, my(cf = polcoef(g, e, vars[i]));
    if(cf != 0, monorec(cf, vars, i+1, if(e == 0, pre, concat([pre, "*", Str(vars[i]), "^", Str(e)])), ~L)));
}
mono(f, vars) = {
  my(L = List());
  monorec(f, vars, 1, "", ~L);
  if(#L == 0, "0", my(s = L[1]); for(k = 2, #L, s = concat([s, "+", L[k]])); s);
}
readvec_txt(fn) ={ my(v = externstr(concat("cat ", fn))); vector(#v, i, eval(v[i])); }
{
  my(arms = ["kummer","subgroup","dickson","dickson-torus","smallx","random","onlycurve","dloginterval"]);
  my(tg = externstr(concat(["tail -n +2 ", DIR, "/targets.csv"])));
  my(tx = vector(min(NT,#tg), i, eval(strsplit(tg[i], ",")[2])));
  my(out = fileopen(concat(DIR, "/grid_counts.csv"), "w"));
  filewrite(out, "arm,m,t,target,xR,grid_count");
  mkdir_cmd = concat(["mkdir -p ", DIR, "/msolve_in"]); system(mkdir_cmd);
  for(t = 2, TMAX2,
    foreach(arms, arm,
      my(fn = concat([DIR, "/base_", arm, "_t", Str(t), ".txt"]));
      if(system(concat(["test -f ", fn])) != 0, next);
      my(V = readvec_txt(fn), Bsz = #V, Lpol = Mod(prod(i=1, Bsz, 'z - V[i]), P));
      for(m = 2, 3,
        if(m == 3 && t > TMAX3, next);
        for(k = 1, #tx,
          my(xR = tx[k], cnt = 0);
          if(setsearch(Set(V), xR), next);
          if(m == 2,
            for(i=1,Bsz, for(j=1,Bsz, if(Mod(S3(V[i],V[j],xR),P) == 0, cnt++))),
            \\ m = 3: S4 = Res_U(S3(x1,x2,U), S3(x3,xR,U)) evaluated pointwise
            my(Q3 = vector(Bsz, l, Mod(S3(V[l], xR, 'U), P)));
            for(i=1,Bsz, for(j=1,Bsz, my(q = Mod(S3(V[i],V[j],'U), P));
              for(l=1,Bsz, if(polresultant(q, Q3[l], 'U) == 0, cnt++)))));
          filewrite(out, Str(arm, ",", m, ",", t, ",", k-1, ",", xR, ",", cnt));
          if((m == 2 && t <= MS2) || (m == 3 && t <= MS),
            my(vars = if(m==2, ['x1,'x2], ['x1,'x2,'x3]));
            my(S = if(m==2, S3('x1,'x2,xR), polresultant(S3('x1,'x2,'U), S3('x3,xR,'U), 'U)));
            my(Ls = vector(m, i, subst(lift(Lpol), 'z, vars[i])));
            my(hdr = concat([if(m==2,"x1,x2","x1,x2,x3"), "\n", Str(P), "\n"]));
            my(base = concat([DIR, "/msolve_in/", arm, "_m", Str(m), "_t", Str(t), "_k", Str(k-1)]));
            my(full = concat([hdr, mono(S, vars)]));
            for(i=1,m, full = concat([full, ",\n", mono(Ls[i], vars)]));
            my(f1 = fileopen(concat(base, "_full.ms"), "w")); filewrite(f1, full); fileclose(f1);
            my(part = concat([hdr, mono(S, vars)]));
            for(i=1,m-1, part = concat([part, ",\n", mono(Ls[i], vars)]));
            my(f2 = fileopen(concat(base, "_partial.ms"), "w")); filewrite(f2, part); fileclose(f2);
            if(arm == "random",
              \\ non-split control: random monic L of degree B, seeded by (t, k)
              setrand(1000*t + k); my(Ln = 'z^Bsz + sum(e=0, Bsz-1, random(P)*'z^e));
              my(nsp = concat([hdr, mono(S, vars)]));
              for(i=1,m, nsp = concat([nsp, ",\n", mono(subst(Ln,'z,vars[i]), vars)]));
              my(f3 = fileopen(concat([DIR, "/msolve_in/random-nonsplit_m", Str(m), "_t", Str(t), "_k", Str(k-1), "_full.ms"]), "w"));
              filewrite(f3, nsp); fileclose(f3);
              my(fr = fileopen(concat([DIR, "/msolve_in/random-nonsplit_m", Str(m), "_t", Str(t), "_k", Str(k-1), "_roots.txt"]), "w"));
              filewrite(fr, Str(#polrootsmod(Ln, P))); fileclose(fr)))))));
  fileclose(out);
}
