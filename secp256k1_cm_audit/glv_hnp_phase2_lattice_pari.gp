\\======================================================================
\\ GLV-aware HNP Phase 2 toy attack -- PARI/GP native LLL
\\======================================================================
\\
\\ Threat model (Phase 2 from RESEARCH_GLV_HNP_PHASE2.md):
\\   k = k1 + lambda * k2 mod n    (GLV decomposition)
\\   k1 in [0, K1_BOUND)            strongly biased (side-channel leak)
\\   k2 in [0, K2_BOUND)            GLV domain bound: K2_BOUND ~= ceil(sqrt(n))
\\
\\ Standard HNP on k_full fails: k_full = k1 + lambda*k2 is nearly uniform
\\ because lambda*k2 dominates.  Phase 2 lattice exploits BOTH bounds:
\\   (K1_BOUND * K2_BOUND / n)^m << 1  for m >= info-theoretic threshold.
\\
\\ Toy curve: y^2 = x^3 + 2 over F_211,  n=199,  lambda=106.
\\ Same curve as glv_hnp_phase2_attack.py.
\\
\\ Run: gp -q glv_hnp_phase2_lattice_pari.gp

default(parisize, 256000000);
default(timer, 0);

P_FIELD  = 211;
B_COEFF  = 2;
N        = 199;
LAM      = 106;
K1_BOUND = 2;
K2_BOUND = 15;

\\ Column-diagonal scaling: all planted-vector entries are O(N) ~200
S_K1     = N \ K1_BOUND;   \\ 99
S_D      = 1;
S_K2     = N \ K2_BOUND;   \\ 13
S_KANNAN = N;               \\ 199

print("=================================================================");
print("GLV-aware HNP Phase 2 -- PARI/GP native LLL");
print("=================================================================");
print("Curve: y^2 = x^3 + 2 over F_211,  n=199,  lambda=106");
print("K1_BOUND=2  K2_BOUND=15  S_K1=99  S_D=1  S_K2=13  S_KANNAN=199");
print("");

\\----------------------------------------------------------------------
\\ EC arithmetic (a=0, y^2 = x^3 + b over F_p)
\\----------------------------------------------------------------------
{ec_add(P1,P2) =
  if(P1==[0,0],return(P2));
  if(P2==[0,0],return(P1));
  local(x1,y1,x2,y2,s,x3,y3);
  x1=P1[1]; y1=P1[2]; x2=P2[1]; y2=P2[2];
  if(x1==x2,
    if((y1+y2)%P_FIELD==0,return([0,0]));
    s=Mod(3*x1^2,P_FIELD)/Mod(2*y1,P_FIELD),
    s=Mod(y2-y1,P_FIELD)/Mod(x2-x1,P_FIELD));
  x3=lift(s^2-Mod(x1+x2,P_FIELD));
  y3=lift(s*(Mod(x1,P_FIELD)-Mod(x3,P_FIELD))-Mod(y1,P_FIELD));
  [x3%P_FIELD, y3%P_FIELD]}

{ec_mul(P,k) =
  local(R,Q); R=[0,0]; Q=P;
  while(k>0, if(k%2==1,R=ec_add(R,Q)); Q=ec_add(Q,Q); k=k\2);
  R}

\\ Generator (first point of order N)
{G=[0,0]; found_G=0;
 for(x=0,P_FIELD-1,
   if(found_G,break);
   rhs=Mod(x^3+B_COEFF,P_FIELD);
   for(y=0,P_FIELD-1,
     if(Mod(y^2,P_FIELD)==rhs,
       Pt=[x,y];
       if(ec_mul(Pt,N)==[0,0],G=Pt;found_G=1;break))))}
if(!found_G,error("generator not found"));
print("G = ",G);

\\----------------------------------------------------------------------
\\ Generate m ECDSA signatures with GLV-domain bias.
\\ Returns vector of [A, B, k1, k2, kfull] with A+B*d = kfull mod n.
\\----------------------------------------------------------------------
{gen_sigs(d_secret,m_sigs,seed_val) =
  local(sigs,k1,k2,kfull,R,r_sig,h_msg,s_sig,A,B);
  sigs=vector(m_sigs); setrand(seed_val);
  i=1; att=0;
  while(i<=m_sigs && att<100000,
    att++;
    k1=random(K1_BOUND); k2=random(K2_BOUND);
    kfull=lift(Mod(k1+LAM*k2,N));
    if(kfull==0,next);
    R=ec_mul(G,kfull);
    if(R==[0,0],next);
    r_sig=R[1]%N; if(r_sig==0,next);
    h_msg=random(N);
    s_sig=lift(Mod(h_msg+d_secret*r_sig,N)/Mod(kfull,N));
    if(s_sig==0,next);
    A=lift(Mod(h_msg,N)/Mod(s_sig,N));
    B=lift(Mod(r_sig,N)/Mod(s_sig,N));
    if(lift(Mod(A+B*d_secret,N))!=kfull,error("sig sanity"));
    sigs[i]=[A,B,k1,k2,kfull]; i++);
  if(i<=m_sigs,error("too few sigs"));
  sigs}

\\----------------------------------------------------------------------
\\ Phase 2 lattice (rows = integer basis vectors, dim = 2m+2).
\\
\\ Rows 1..m       : mod-n:   N*S_K1 at column i
\\ Row  m+1        : d-row:   B_i*S_K1 at col i,  S_D at col m+1
\\ Rows m+2..2m+1  : k2-rows: -LAM*S_K1 at col i,  S_K2 at col m+1+i
\\ Row  2m+2       : Kannan:  A_i*S_K1 at col i,  S_KANNAN at last col
\\----------------------------------------------------------------------
{build_lattice(sigs) =
  local(m,dim,M);
  m=#sigs; dim=2*m+2; M=matrix(dim,dim);
  for(i=1,m, M[i,i]=N*S_K1);
  for(i=1,m, M[m+1,i]=sigs[i][2]*S_K1);
  M[m+1,m+1]=S_D;
  for(i=1,m, M[m+1+i,i]=-LAM*S_K1; M[m+1+i,m+1+i]=S_K2);
  for(i=1,m, M[2*m+2,i]=sigs[i][1]*S_K1);
  M[2*m+2,dim]=S_KANNAN;
  M}

\\----------------------------------------------------------------------
\\ Verify planted combination is in the lattice and gives the right vector.
\\ Combination:  -q_i*(row i)  +  d*(row m+1)  +  k2_i*(row m+1+i)  +  row 2m+2
\\ Expected entries: k1_i*S_K1, d*S_D, k2_i*S_K2, S_KANNAN.
\\----------------------------------------------------------------------
{verify_planted(M,sigs,d_secret) =
  local(m,dim,v,ured,qi,ok);
  m=#sigs; dim=2*m+2; v=vector(dim);
  for(i=1,m,
    ured=sigs[i][1]+sigs[i][2]*d_secret-LAM*sigs[i][4];
    qi=(ured-sigs[i][3])\N;
    if(ured-sigs[i][3]!=qi*N,print("qi not integer"); return(0));
    for(j=1,dim,v[j]-=qi*M[i,j]));
  for(j=1,dim,v[j]+=d_secret*M[m+1,j]);
  for(i=1,m, for(j=1,dim,v[j]+=sigs[i][4]*M[m+1+i,j]));
  for(j=1,dim,v[j]+=M[2*m+2,j]);
  ok=1;
  for(i=1,m,
    if(v[i]!=sigs[i][3]*S_K1,print("MISMATCH k1_",i);ok=0));
  if(v[m+1]!=d_secret*S_D,print("MISMATCH d-slot");ok=0);
  for(i=1,m,
    if(v[m+1+i]!=sigs[i][4]*S_K2,print("MISMATCH k2_",i);ok=0));
  if(v[dim]!=S_KANNAN,print("MISMATCH Kannan");ok=0);
  ok}

\\----------------------------------------------------------------------
\\ LLL-reduce (rows as basis).
\\ qflll(A,1): A has COLUMNS as basis; returns T s.t. A*T is LLL-reduced.
\\ M has rows => feed M~, then (M~*T)~ has rows = reduced vectors.
\\----------------------------------------------------------------------
{lll_reduce_rows(M) = (M~*qflll(M~,1))~}

\\----------------------------------------------------------------------
\\ Scan LLL-reduced rows for d.
\\ Kannan marker: last entry = ±S_KANNAN.
\\ d-slot at column m+1: d_cand = ±row[m+1] mod N.
\\----------------------------------------------------------------------
{recover_d(Mred,m,d_secret) =
  local(dim,v,last,sgn,dc);
  dim=2*m+2;
  for(row=1,matsize(Mred)[1],
    v=vector(dim,j,Mred[row,j]);
    last=v[dim];
    if(abs(last)!=S_KANNAN,next);
    sgn=if(last>0,1,-1);
    dc=lift(Mod(sgn*v[m+1],N));
    if(dc==0,next);
    if(dc==d_secret,return(dc)));
  -1}

\\----------------------------------------------------------------------
\\ One experiment: generate sigs, build lattice, LLL, recover.
\\----------------------------------------------------------------------
{run_experiment(d_secret,m_val,seed_val,verbose) =
  local(sigs,M,Mred,d_rec,v,pvn);
  sigs=gen_sigs(d_secret,m_val,seed_val);
  M=build_lattice(sigs);
  if(!verify_planted(M,sigs,d_secret),
    print("ERROR: planted not in lattice"); return(0));
  if(verbose,
    \\ compute planted vector and its norm
    v=vector(2*m_val+2);
    for(i=1,m_val,
      ured=sigs[i][1]+sigs[i][2]*d_secret-LAM*sigs[i][4];
      qi=(ured-sigs[i][3])\N;
      for(j=1,2*m_val+2,v[j]-=qi*M[i,j]));
    for(j=1,2*m_val+2,v[j]+=d_secret*M[m_val+1,j]);
    for(i=1,m_val,for(j=1,2*m_val+2,v[j]+=sigs[i][4]*M[m_val+1+i,j]));
    for(j=1,2*m_val+2,v[j]+=M[2*m_val+2,j]);
    pvn=floor(sqrt(sum(j=1,2*m_val+2,v[j]^2))+0.5);
    print("  Signatures (m=",m_val,", seed=",seed_val,", d=",d_secret,"):");
    for(i=1,m_val,
      print("    sig",i,": k1=",sigs[i][3]," k2=",sigs[i][4]," k=",sigs[i][5]));
    print("  Planted vector: [",concat(vector(2*m_val+2,j,Str(v[j])),","),"]");
    print("  Planted norm: ",pvn));
  Mred=lll_reduce_rows(M);
  d_rec=recover_d(Mred,m_val,d_secret);
  if(verbose,
    if(d_rec==d_secret,
      print("  RECOVERED d = ",d_rec,"  ✓"),
      print("  NOT recovered (d = ",d_secret,")")));
  d_rec==d_secret}

\\======================================================================
\\ Main experiment
\\======================================================================

\\ --- Single verbose run at m=6 (reliable threshold) ---
print("---- Single verbose run: m=6, seed=42, d=104 ----");
print("");
ok1 = run_experiment(104, 6, 42, 1);
print("Result: ", if(ok1,"RECOVERED ✓","NOT RECOVERED ✗"));
print("");

\\ --- Confirmation run at m=6, second seed ---
print("---- Confirmation: m=6, seed=1234, d=77 ----");
print("");
ok2 = run_experiment(77, 6, 1234, 1);
print("Result: ", if(ok2,"RECOVERED ✓","NOT RECOVERED ✗"));
print("");

\\ --- Sweep over m=3..7, 5 seeds ---
print("---- Sweep: m=3..7, five seeds each ----");
print("");
seeds3 = [42, 1234, 9999, 7, 314159];
dvals3 = [104, 77, 150, 33, 191];
{for(mm=3,7,
  wins=0;
  for(ss=1,5,
    wins+=run_experiment(dvals3[ss],mm,seeds3[ss],0));
  print("m=",mm,": ",wins,"/5 seeds recovered d"))}
print("");

print("=================================================================");
print("Summary:");
print("  n=199, K1_BOUND=2, K2_BOUND=15, lattice dim=2m+2");
ith=ceil(log(N+0.0)/log(N*1.0/(K1_BOUND*K2_BOUND)));
print("  Info-theoretic threshold: m >= ", ith);
print("  LLL: PARI qflll(M~,1)  (integer LLL, flag=1)");
print("=================================================================");
