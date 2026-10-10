/* nist_rho.cuh -- Pollard-rho r-adding walk for the NIST prime curves,
 * shared by the CUDA kernel, the host replay/solver and the CPU test.
 * Generic over the curve (hence the field and its limb count F::N), so the
 * same code runs on P-256, P-384 and the toy curve of the end-to-end test.
 *
 * Walk (deterministic function of the point; the single step and the
 * batched step agree bit-for-bit):
 *
 *   partition(P)  = x.v[0] & (R-1)                      R = 2^r_bits
 *   step          P <- canonical(P + M[partition(P)])   M[j] = c_j G + d_j Q
 *   distinguished ((x.v[0] >> 8) & dp_mask) == 0
 *
 * canonical() folds by the negation map (fold 2): of {P,-P} keep the one
 * whose internal y is lexicographically smaller.  This makes the walk a
 * function on E/{+-1}, cutting the expected steps to a collision from
 * sqrt(pi n / 2) to sqrt(pi n / 4).  NIST prime curves have j != 0, so +-1
 * is all of Aut(E) and sqrt(2) is the ceiling for this kind of folding.
 * fold 1 (canonical = identity) is kept for the control measurement.
 *
 * x, y are the field's internal representation (Montgomery form on the
 * device fields).  Hashing/comparing it saves a conversion per step and is
 * still a deterministic function of the point.
 *
 * canonical() returns an AUT CODE (here just `negated`); the host replay
 * multiplies a walk's (a,b) coefficients by its scalar action (-1)^negated
 * so a reconstructed point's logarithm stays correct.
 *
 * FRUITLESS CYCLES.  The negation-map fold closes short even-length cycles
 * that carry no information; the commonest is length 2 (same table index
 * twice with a [-1] in between), at rate 1/(fold R) per step.  Each walk
 * keeps the hashes of its last RHO_CYCLE_DEPTH points; a repeat means a
 * cycle of length 2..DEPTH+1 whose members were all just seen, and the walk
 * escapes by DOUBLING THE MEMBER WITH THE SMALLEST HASH -- a function of the
 * cycle alone, so two walks trapped in it leave identically and still
 * collide.  Longer cycles are left to max_steps.
 *
 * A walk is identified by (walk index, restart); its start is a*G+b*Q with
 * (a,b) from a fixed PRNG of that pair, so a distinguished point is reported
 * as {x, walk, restart, steps} and the host recovers (a,b) by replaying the
 * two colliding walks only.  A walk restarts after each DP and after
 * max_steps steps without one.
 */
#ifndef GPU_NIST_RHO_CUH
#define GPU_NIST_RHO_CUH
#include "nist_curve.cuh"

#define RHO_MAX_RBITS 8
#define RHO_DP_SHIFT  8
#define RHO_FOLD_NONE 1
#define RHO_FOLD_NEG  2
#ifndef RHO_CYCLE_DEPTH
#define RHO_CYCLE_DEPTH 5
#endif

struct rho_params {
 uint32_t r_bits;      /* log2 table size, <= RHO_MAX_RBITS */
 uint32_t dp_mask;     /* distinguished iff ((x.v[0]>>8) & dp_mask) == 0 */
 uint32_t fold;        /* RHO_FOLD_NONE or RHO_FOLD_NEG */
 uint32_t max_steps;   /* abandon a walk after this many steps */
 uint32_t table_seed;  /* seed for the c_j,d_j table coefficients */
 uint32_t order_bits;  /* bit-width the seed scalars are masked to */
};

template<class F> struct rho_dp {
 uint32_t x[F::N];
 uint32_t walk, restart, steps, pad;
};

template<class C> struct rho_state {
 nap<typename C::Field> P;
 uint32_t h[RHO_CYCLE_DEPTH];
 uint32_t escape;
};

/* ---- deterministic seeds (splitmix64) --------------------------------- */
NF_HD uint64_t rho_splitmix64(uint64_t&s){
 uint64_t z=(s+=0x9E3779B97F4A7C15ull);
 z=(z^(z>>30))*0xBF58476D1CE4E5B9ull;
 z=(z^(z>>27))*0x94D049BB133111EBull;
 return z^(z>>31);
}
/* An order_bits-wide scalar into a length-N limb array (a*P depends only on
 * a mod n, so masking to the order width spreads toy scalars well and never
 * changes a NIST walk). */
template<int N> NF_HD void rho_scalar_from_seed(uint64_t&s,uint32_t out[N],uint32_t order_bits){
 for(int i=0;i<N;i++) out[i]=0;
 for(int i=0;i<(N+1)/2;i++){uint64_t z=rho_splitmix64(s); if(2*i<N)out[2*i]=(uint32_t)z; if(2*i+1<N)out[2*i+1]=(uint32_t)(z>>32);}
 for(int l=0;l<N;l++){int lo=32*l; if((uint32_t)lo>=order_bits)out[l]=0; else if(lo+32>(int)order_bits)out[l]&=(1u<<(order_bits-lo))-1u;}
}
template<int N> NF_HD void rho_walk_seed(uint32_t walk,uint32_t restart,uint32_t a[N],uint32_t b[N],uint32_t order_bits){
 uint64_t s=((uint64_t)walk<<32|restart)*0xD1B54A32D192ED03ull+0x8CB92BA72F3D8DD7ull;
 rho_scalar_from_seed<N>(s,a,order_bits); rho_scalar_from_seed<N>(s,b,order_bits);
}
template<int N> NF_HD void rho_table_seed(uint32_t table_seed,uint32_t j,uint32_t c[N],uint32_t d[N],uint32_t order_bits){
 uint64_t s=((uint64_t)table_seed<<32|j)*0xA0761D6478BD642Full+0xE7037ED1A0B428DBull;
 rho_scalar_from_seed<N>(s,c,order_bits); rho_scalar_from_seed<N>(s,d,order_bits);
}

template<class C> NF_HD uint32_t rho_partition(const nap<typename C::Field>&P,const rho_params&prm){return P.x.v[0]&((1u<<prm.r_bits)-1u);}
template<class C> NF_HD int rho_is_dp(const nap<typename C::Field>&P,const rho_params&prm){return ((P.x.v[0]>>RHO_DP_SHIFT)&prm.dp_mask)==0;}
template<class C> NF_HD uint32_t rho_hash(const nap<typename C::Field>&P){return P.x.v[0];}

/* Fold P onto its {+-1} representative; return 1 if it negated y. */
template<class C> NF_HD int rho_canonical(nap<typename C::Field>&P,const rho_params&prm){
 typedef typename C::Field F;
 if(prm.fold==RHO_FOLD_NONE||P.inf) return 0;
 typename F::elt ny=F::neg(P.y);
 uint32_t f=(uint32_t)F::lt(ny,P.y);
 F::cmov(P.y,ny,f);
 return (int)f;
}
template<class C> NF_HD nap<typename C::Field> rho_apply_aut(const nap<typename C::Field>&P,int code){
 nap<typename C::Field> r=P; if(!P.inf && (code&1)) r.y=C::Field::neg(P.y); return r;
}

template<class C> NF_HD void rho_reset_cycle(rho_state<C>&st,const nap<typename C::Field>&P){
 uint32_t h=rho_hash<C>(P);
 for(int i=0;i<RHO_CYCLE_DEPTH;i++) st.h[i]=h;
 st.escape=0;
}

/* phase A: choose the addend and the denominator to invert. */
#define RHO_MODE_ADD 0
#define RHO_MODE_DOUBLE 1
#define RHO_MODE_INF 2
#define RHO_MODE_ESCAPE 3
#define RHO_STEP_ADVANCED 1
#define RHO_STEP_SEEK 0

template<class C> NF_HD int rho_phase_a(const rho_state<C>&st,const nap<typename C::Field>*table,
                                        const rho_params&prm,typename C::Field::elt&den,uint32_t&j_out){
 typedef typename C::Field F;
 if(st.escape && rho_hash<C>(st.P)==st.h[0]){
  den=F::dbl(st.P.y); j_out=0;
  if(F::is_zero(den)) return RHO_MODE_INF;
  return RHO_MODE_ESCAPE;
 }
 uint32_t j=rho_partition<C>(st.P,prm); j_out=j;
 typename F::elt d=F::sub(table[j].x,st.P.x);
 if(!F::is_zero(d)){den=d; return RHO_MODE_ADD;}
 if(F::eq(st.P.y,table[j].y)){
  den=F::dbl(st.P.y);
  if(F::is_zero(den)) return RHO_MODE_INF;
  return RHO_MODE_DOUBLE;
 }
 den=F::one(); return RHO_MODE_INF;
}

/* phase B: finish the step given inv = 1/den.  Returns RHO_STEP_ADVANCED on
 * a new counted point, RHO_STEP_SEEK while walking a detected cycle, or
 * -(length) when a cycle of that length was just detected. */
template<class C> NF_HD int rho_phase_b(rho_state<C>&st,const nap<typename C::Field>*table,const rho_params&prm,
                                        int mode,uint32_t j,const typename C::Field::elt&inv,int*aut_out){
 typedef typename C::Field F;
 if(mode==RHO_MODE_ESCAPE){
  nap<F> nxt=C::affine_add_with_inv(st.P,st.P,inv,1);
  *aut_out=rho_canonical<C>(nxt,prm);
  rho_reset_cycle<C>(st,st.P); st.P=nxt; return RHO_STEP_ADVANCED;
 }
 nap<F> nxt=C::affine_add_with_inv(st.P,table[j],inv,mode==RHO_MODE_DOUBLE);
 *aut_out=rho_canonical<C>(nxt,prm);
 uint32_t hn=rho_hash<C>(nxt);
 if(st.escape){ st.P=nxt; if(--st.escape==0) rho_reset_cycle<C>(st,nxt); return RHO_STEP_SEEK; }
 if(prm.fold!=RHO_FOLD_NONE){
  int len=0;
  for(int i=RHO_CYCLE_DEPTH-1;i>=0;i--) if(hn==st.h[i]) len=i+2;
  if(len){
   uint32_t hp=rho_hash<C>(st.P), target=hn<hp?hn:hp;
   for(int i=0;i<RHO_CYCLE_DEPTH;i++) if(i+2<len && st.h[i]<target) target=st.h[i];
   st.h[0]=target; st.P=nxt; st.escape=RHO_CYCLE_DEPTH+1; return -len;
  }
 }
 for(int i=RHO_CYCLE_DEPTH-1;i>0;i--) st.h[i]=st.h[i-1];
 st.h[0]=rho_hash<C>(st.P); st.P=nxt; return RHO_STEP_ADVANCED;
}

/* One unbatched step, no cycle handling: the reference primitive. */
template<class C> NF_HD nap<typename C::Field> rho_step_single(const nap<typename C::Field>&P,
                                 const nap<typename C::Field>*table,const rho_params&prm,int*j_out,int*aut_out){
 typedef typename C::Field F;
 uint32_t j=rho_partition<C>(P,prm); const nap<F>&M=table[j]; nap<F> r;
 if(F::eq(P.x,M.x)){
  if(F::eq(P.y,M.y)){ r=C::affine_add_with_inv(P,M,F::inv(F::dbl(P.y)),1); }
  else { r=C::infinity(); }
 } else { r=C::affine_add_with_inv(P,M,F::inv(F::sub(M.x,P.x)),0); }
 *j_out=(int)j; *aut_out=rho_canonical<C>(r,prm); return r;
}

/* ---- batched walk state (SoA): limb l of walk i at X[l*nwalks+i] -------- */
template<class C> struct rho_ctx {
 uint32_t *X,*Y;           /* [F::N][nwalks] */
 uint32_t *H;              /* [RHO_CYCLE_DEPTH][nwalks] */
 uint32_t *esc,*steps,*restarts;   /* [nwalks] each */
 uint32_t nthreads, walks_per_thread;
 const nap<typename C::Field> *table;
 nap<typename C::Field> G,Q;
 rho_params prm;
 rho_dp<typename C::Field> *dp_out;
 uint32_t *dp_count, dp_cap;
 unsigned long long *cycles, *aborts;   /* optional telemetry */
};
template<class C> NF_HD uint32_t rho_nwalks(const rho_ctx<C>&c){return c.nthreads*c.walks_per_thread;}

template<class C> NF_HD void rho_load(const rho_ctx<C>&c,uint32_t idx,rho_state<C>&st){
 typedef typename C::Field F; uint32_t n=rho_nwalks(c);
 for(int l=0;l<F::N;l++){st.P.x.v[l]=c.X[l*n+idx]; st.P.y.v[l]=c.Y[l*n+idx];}
 st.P.inf=0;
 for(int i=0;i<RHO_CYCLE_DEPTH;i++) st.h[i]=c.H[i*n+idx];
 st.escape=c.esc[idx];
}
template<class C> NF_HD void rho_store(const rho_ctx<C>&c,uint32_t idx,const rho_state<C>&st){
 typedef typename C::Field F; uint32_t n=rho_nwalks(c);
 for(int l=0;l<F::N;l++){c.X[l*n+idx]=st.P.x.v[l]; c.Y[l*n+idx]=st.P.y.v[l];}
 for(int i=0;i<RHO_CYCLE_DEPTH;i++) c.H[i*n+idx]=st.h[i];
 c.esc[idx]=st.escape;
}
NF_HD void rho_count(unsigned long long*ctr){
#ifdef __CUDA_ARCH__
 atomicAdd(ctr,1ull);
#else
 if(ctr)(*ctr)++;
#endif
}
template<class C> NF_HD void rho_emit_dp(const rho_ctx<C>&c,const nap<typename C::Field>&P,uint32_t idx){
 typedef typename C::Field F;
#ifdef __CUDA_ARCH__
 uint32_t slot=atomicAdd(c.dp_count,1u);
#else
 uint32_t slot=(*c.dp_count)++;
#endif
 if(slot<c.dp_cap){ rho_dp<F>&d=c.dp_out[slot]; for(int l=0;l<F::N;l++)d.x[l]=P.x.v[l]; d.walk=idx; d.restart=c.restarts[idx]; d.steps=c.steps[idx]; d.pad=0; }
}

/* (Re)seed walk idx from (idx,restart): P <- a*G + b*Q, folded.  Loops with
 * a fresh restart in the vanishing case a*G+b*Q = O. */
template<class C> NF_HD void rho_seed_walk(const rho_ctx<C>&c,uint32_t idx,rho_state<C>&st){
 typedef typename C::Field F;
 for(;;){
  uint32_t a[F::N],b[F::N];
  rho_walk_seed<F::N>(idx,c.restarts[idx],a,b,c.prm.order_bits);
  nap<F> S=C::double_scalar_mul(c.G,a,c.Q,b,c.prm.order_bits);
  if(!S.inf){ int code; (void)code; rho_canonical<C>(S,c.prm); st.P=S; rho_reset_cycle<C>(st,S); return; }
  c.restarts[idx]++;
 }
}

/* Advance every walk one logical step, batching the inversion across the
 * warp-worth of walks this thread owns.  Matches rho_step_single on a plain
 * step; also does cycle escape and DP emission. */
template<class C> NF_HD void rho_batch_step(const rho_ctx<C>&c,uint32_t tid,
                                            typename C::Field::elt*den,typename C::Field::elt*scratch){
 typedef typename C::Field F;
 uint32_t W=c.walks_per_thread;
 rho_state<C> st; int mode[64]; uint32_t jj[64];   /* W <= 64 */
 /* phase A for each owned walk; collect denominators */
 for(uint32_t w=0;w<W;w++){
  uint32_t idx=tid+w*c.nthreads;
  rho_load(c,idx,st);
  typename F::elt d; uint32_t j;
  int m=rho_phase_a<C>(st,c.table,c.prm,d,j);
  mode[w]=m; jj[w]=j;
  if(m==RHO_MODE_INF){ den[w]=F::one(); }
  else den[w]=d;
 }
 F::batch_inv(den,W,scratch);
 /* phase B for each */
 for(uint32_t w=0;w<W;w++){
  uint32_t idx=tid+w*c.nthreads;
  rho_load(c,idx,st);
  int m=mode[w];
  if(m==RHO_MODE_INF){ c.restarts[idx]++; rho_seed_walk(c,idx,st); c.steps[idx]=0; rho_store(c,idx,st); continue; }
  int aut; int r=rho_phase_b<C>(st,c.table,c.prm,m,jj[w],den[w],&aut);
  if(r<0){ rho_count(c.cycles? &c.cycles[(-r)-2]:0); rho_store(c,idx,st); continue; }
  if(r==RHO_STEP_SEEK){ rho_store(c,idx,st); continue; }
  /* advanced: count, test DP, maybe reseed */
  uint32_t s=++c.steps[idx];
  if(rho_is_dp<C>(st.P,c.prm)){
   rho_emit_dp(c,st.P,idx);
   c.restarts[idx]++; rho_seed_walk(c,idx,st); c.steps[idx]=0;
  } else if(s>=c.prm.max_steps){
   rho_count(c.aborts); c.restarts[idx]++; rho_seed_walk(c,idx,st); c.steps[idx]=0;
  }
  rho_store(c,idx,st);
 }
}

#endif
