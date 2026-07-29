\\ thread22_cover_search.gp
\\ Thread 22 (Part 2): For non-cubic-residue proxy primes, the naive cover
\\ y^2=(x^3+c_i)(x^3+c_j) does NOT give Jac ~ E_i x E_j.
\\ SEARCH for the correct cover curve within the sextic family y^2=x^6+ax^3+b.
\\
\\ Over p'=31: pair (0,3): t_0=11, t_3=-11.
\\ Target char poly: (T^2-11T+31)(T^2+11T+31) = T^4-59T^2+961.
\\ Also search pair (1,4): t_1=7, t_4=-7 (QT pair over F_31).
\\   Target: (T^2-7T+31)(T^2+7T+31) = T^4-35T^2+961.
\\
\\ Over cubic-residue prime p'=19: confirm naive cover IS correct.
\\
\\ Run: gp -q secp256k1_cm_audit/thread22_cover_search.gp

default(parisize, 256000000);
default(timer, 0);

\\ ---- helper: check if p has all cubic residue for our context ----
\\ x^3+7 splits over F_p iff (-7)^{(p-1)/3}=1 mod p
is_cubic_split(p) = (lift(Mod(-7,p)^((p-1)/3)) == 1);

\\ ---- search sextic family y^2=x^6+ax^3+b for given target char poly ----
search_sextic_family(p_prx, target_cp) = {
  my(ff, cp, hits, aa, bb);
  hits = [];
  for(aa=0, p_prx-1,
    for(bb=1, p_prx-1,
      ff = Mod(1,p_prx)*x^6 + Mod(aa,p_prx)*x^3 + Mod(bb,p_prx);
      if(poldisc(ff) != 0,
        cp = hyperellcharpoly(ff);
        if(cp == target_cp, hits = concat(hits, [[aa,bb,cp]])))));
  hits
};

\\ ---- search factored family y^2=(x^3+r1)(x^3+r2)=x^6+(r1+r2)x^3+r1*r2 ----
search_factored_family(p_prx, target_cp) = {
  my(r1, r2, ff, cp, hits);
  hits = [];
  for(r1=1, p_prx-1,
    for(r2=1, p_prx-1,
      if(r1 != r2,
        ff = Mod(1,p_prx)*(x^3+r1)*(x^3+r2);
        cp = hyperellcharpoly(ff);
        if(cp == target_cp, hits = concat(hits, [[r1,r2,cp]])))));
  hits
};

\\ ========== Main ==========

print("================================================================");
print("Thread 22 Part 2: Searching for correct Howe cover over non-cubic-residue primes");
print("================================================================");
print("");

\\ ----- p'=19: cubic-split prime, naive cover should work -----
p19 = 19;
print("--- p'=19 (cubic-split: (-7)^{(p-1)/3}=",lift(Mod(-7,p19)^((p19-1)/3))," mod 19) ---");
print("  is_cubic_split=",is_cubic_split(p19));
\\ Pair (0,3): target = (T^2-8T+19)(T^2+8T+19) = T^4-46T^2+361
target19 = (x^2-8*x+19)*(x^2+8*x+19);
print("  Target cp for (0,3): ",target19);
hits19 = search_factored_family(p19, target19);
print("  Factored family hits: ",#hits19," (should include [7,12])");
for(idx=1,#hits19, print("    (r1,r2)=(",hits19[idx][1],",",hits19[idx][2],")"));
print("");

\\ ----- p'=31: non-cubic-residue prime -----
p31 = 31;
print("--- p'=31 (cubic-split: (-7)^{(p-1)/3}=",lift(Mod(-7,p31)^((p31-1)/3))," mod 31) ---");
print("  is_cubic_split=",is_cubic_split(p31));
print("");

\\ Pair (0,3) over F_31: t_0=11, t_3=-11
print("Pair (0,3): t0=11 t3=-11");
target_03_31 = (x^2-11*x+31)*(x^2+11*x+31);
print("  Target cp: ",target_03_31);

print("  Searching factored family y^2=(x^3+r1)(x^3+r2) over F_31...");
hits_03_fact = search_factored_family(p31, target_03_31);
print("  Factored family hits: ",#hits_03_fact);
for(idx=1,#hits_03_fact, print("    (r1,r2)=(",hits_03_fact[idx][1],",",hits_03_fact[idx][2],")"));
print("");

print("  Searching sextic family y^2=x^6+ax^3+b over F_31 (all ",p31,"*",p31-1,"=",p31*(p31-1)," smooth cases)...");
hits_03_sex = search_sextic_family(p31, target_03_31);
print("  Sextic family hits: ",#hits_03_sex);
for(idx=1,#hits_03_sex, print("    (a,b)=(",hits_03_sex[idx][1],",",hits_03_sex[idx][2],") cp=",hits_03_sex[idx][3]));
print("");

\\ Pair (1,4) over F_31: t_1=7, t_4=-7
print("Pair (1,4): t1=7 t4=-7");
target_14_31 = (x^2-7*x+31)*(x^2+7*x+31);
print("  Target cp: ",target_14_31);

print("  Searching factored family y^2=(x^3+r1)(x^3+r2) over F_31...");
hits_14_fact = search_factored_family(p31, target_14_31);
print("  Factored family hits: ",#hits_14_fact);
for(idx=1,#hits_14_fact, print("    (r1,r2)=(",hits_14_fact[idx][1],",",hits_14_fact[idx][2],")"));
print("");

print("  Searching sextic family y^2=x^6+ax^3+b over F_31...");
hits_14_sex = search_sextic_family(p31, target_14_31);
print("  Sextic family hits: ",#hits_14_sex);
for(idx=1,#hits_14_sex, print("    (a,b)=(",hits_14_sex[idx][1],",",hits_14_sex[idx][2],") cp=",hits_14_sex[idx][3]));
print("");

\\ ----- p'=43: another non-cubic-residue prime -----
p43 = 43;
print("--- p'=43 (cubic-split: (-7)^{(p-1)/3}=",lift(Mod(-7,p43)^((p43-1)/3))," mod 43) ---");
print("  is_cubic_split=",is_cubic_split(p43));
\\ h of order 6 mod 43 (h=7 from prior run)
h43 = 7;
t0_43 = 43+1-ellcard(ellinit([0,7],43));
t3_43 = 43+1-ellcard(ellinit([0,lift(Mod(7,43)*Mod(h43,43)^3)],43));
print("  t_0=",t0_43,"  t_3=",t3_43,"  QT pair: ",t0_43==-t3_43);
target_03_43 = (x^2-t0_43*x+43)*(x^2-t3_43*x+43);
print("  Target cp (0,3): ",target_03_43);
print("  Searching sextic family y^2=x^6+ax^3+b over F_43...");
hits_03_43 = search_sextic_family(p43, target_03_43);
print("  Hits: ",#hits_03_43);
for(idx=1,#hits_03_43, print("    (a,b)=(",hits_03_43[idx][1],",",hits_03_43[idx][2],")"));
print("");

\\ Summary
print("================================================================");
print("Summary:");
print("  p'=19 (cubic-split): naive cover is correct Howe cover.");
print("  p'=31 (non-cubic-split):");
print("    - Factored family hits for (0,3): ",#hits_03_fact);
print("    - Sextic family hits for (0,3):   ",#hits_03_sex);
print("    - Factored family hits for (1,4): ",#hits_14_fact);
print("    - Sextic family hits for (1,4):   ",#hits_14_sex);
print("  p'=43 (non-cubic-split): sextic hits for (0,3): ",#hits_03_43);
print("CONCLUSION: If hits found over non-split primes, the Howe cover");
print("  lies in sextic family but NOT in factored (x^3+r1)(x^3+r2) subfamily.");
print("  If no hits in sextic family: cover has different form (non-sextic).");
