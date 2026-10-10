/* Host driver for the NIST rho walk: build the jump table, run the walk on
 * the CPU while carrying each walk's (a,b) with P = a*G + b*Q, and solve the
 * logarithm on the first collision between two walks with different
 * coefficients.  Used by the end-to-end toy test and as the reference the
 * device kernel's distinguished points are checked against.
 *
 * The coefficient update mirrors rho_phase_b exactly: every non-INF step
 * either adds the table entry's (c_j,d_j) or, on an escape, doubles (a,b);
 * a fold then negates both.  So the point and its coefficients never drift,
 * and a distinguished point carries an exact logarithm relation.
 *
 * Modular scalar arithmetic is mod the group order n, kept in 64-bit here
 * (the end-to-end curve is ~2^20); n must be < 2^63.  The device solver for
 * the big curves would carry 256/384-bit coefficients instead -- out of
 * scope for this CPU verification, which never runs a full-size solve.
 */
#ifndef GPU_NIST_RHO_HOST_HPP
#define GPU_NIST_RHO_HOST_HPP
#include "nist_rho.cuh"
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

struct Zn {
 unsigned long long n;
 Zn(unsigned long long n_=1):n(n_){}
 unsigned long long add(unsigned long long a,unsigned long long b)const{a%=n;b%=n;unsigned long long s=a+b;if(s>=n||s<a)s-=n;return s%n;}
 unsigned long long sub(unsigned long long a,unsigned long long b)const{a%=n;b%=n;return a>=b?a-b:a+n-b;}
 unsigned long long mul(unsigned long long a,unsigned long long b)const{return (unsigned long long)(( (unsigned __int128)(a%n)*(b%n))%n);}
 unsigned long long neg(unsigned long long a)const{a%=n;return a?n-a:0;}
 unsigned long long pw(unsigned long long a,unsigned long long e)const{unsigned long long r=1%n;a%=n;while(e){if(e&1)r=mul(r,a);a=mul(a,a);e>>=1;}return r;}
 unsigned long long inv(unsigned long long a)const{return pw(a%n,n-2);}  /* n prime */
 template<int N> unsigned long long from_limbs(const uint32_t*l)const{
  unsigned long long v=0, base=1;
  for(int i=0;i<N;i++){v=add(v,mul(base,(unsigned long long)l[i]%n)); base=mul(base,0x100000000ull%n);}
  return v;
 }
};

template<class C> struct RhoHost {
 typedef typename C::Field F; typedef nap<F> pt;
 struct TableEntry { pt P; unsigned long long c,d; };

 /* M[j] = c_j*G + d_j*Q, with (c_j,d_j) also reduced mod n for replay. */
 static std::vector<TableEntry> build_table(const pt&G,const pt&Q,const rho_params&prm,const Zn&zn){
  uint32_t R=1u<<prm.r_bits; std::vector<TableEntry> t(R);
  for(uint32_t j=0;j<R;j++){
   uint32_t c[F::N],d[F::N];
   rho_table_seed<F::N>(prm.table_seed,j,c,d,prm.order_bits);
   t[j].P=C::double_scalar_mul(G,c,Q,d,prm.order_bits);
   t[j].c=zn.template from_limbs<F::N>(c); t[j].d=zn.template from_limbs<F::N>(d);
  }
  return t;
 }

 struct Walk { rho_state<C> st; unsigned long long a,b; };

 static void seed(Walk&w,uint32_t idx,uint32_t restart,const pt&G,const pt&Q,const rho_params&prm,const Zn&zn){
  for(;;){
   uint32_t ab[F::N],bb[F::N];
   rho_walk_seed<F::N>(idx,restart,ab,bb,prm.order_bits);
   pt S=C::double_scalar_mul(G,ab,Q,bb,prm.order_bits);
   if(!S.inf){
    w.a=zn.template from_limbs<F::N>(ab); w.b=zn.template from_limbs<F::N>(bb);
    int code=rho_canonical<C>(S,prm);
    if(code&1){w.a=zn.neg(w.a);w.b=zn.neg(w.b);}
    w.st.P=S; rho_reset_cycle<C>(w.st,S); return;
   }
   restart++;
  }
 }

 /* Run until a collision yields k with k*G == Q, or `budget` steps pass.
  * Returns true and sets k on success. */
 static bool solve(const pt&G,const pt&Q,const rho_params&prm,unsigned long long n,
                   uint32_t nwalks,unsigned long long budget,unsigned long long&k_out,
                   unsigned long long*steps_out=0){
  Zn zn(n);
  std::vector<TableEntry> tbl=build_table(G,Q,prm,zn);
  std::vector<pt> tptr(tbl.size()); for(size_t i=0;i<tbl.size();i++)tptr[i]=tbl[i].P;
  std::vector<Walk> wk(nwalks);
  std::vector<uint32_t> restart(nwalks,0), wsteps(nwalks,0);
  for(uint32_t i=0;i<nwalks;i++) seed(wk[i],i,restart[i],G,Q,prm,zn);
  struct DP { unsigned long long a,b; };
  std::unordered_map<std::string,DP> seen;
  unsigned long long total=0;
  std::vector<typename F::elt> den(nwalks),scratch(nwalks);
  std::vector<int> mode(nwalks); std::vector<uint32_t> jj(nwalks);
  while(total<budget){
   for(uint32_t i=0;i<nwalks;i++){
    typename F::elt d; uint32_t j; int m=rho_phase_a<C>(wk[i].st,tptr.data(),prm,d,j);
    mode[i]=m; jj[i]=j; den[i]=(m==RHO_MODE_INF)?F::one():d;
   }
   F::batch_inv(den.data(),nwalks,scratch.data());
   for(uint32_t i=0;i<nwalks;i++){
    if(mode[i]==RHO_MODE_INF){ restart[i]++; wsteps[i]=0; seed(wk[i],i,restart[i],G,Q,prm,zn); continue; }
    unsigned long long ca=tbl[jj[i]].c, cb=tbl[jj[i]].d; bool escape=(mode[i]==RHO_MODE_ESCAPE);
    int aut; int r=rho_phase_b<C>(wk[i].st,tptr.data(),prm,mode[i],jj[i],den[i],&aut);
    if(escape){ wk[i].a=zn.add(wk[i].a,wk[i].a); wk[i].b=zn.add(wk[i].b,wk[i].b); }
    else { wk[i].a=zn.add(wk[i].a,ca); wk[i].b=zn.add(wk[i].b,cb); }
    if(aut&1){ wk[i].a=zn.neg(wk[i].a); wk[i].b=zn.neg(wk[i].b); }
    total++;
    if(r<0||r==RHO_STEP_SEEK) continue;        /* inside a fruitless cycle: not counted as a walk step */
    if(rho_is_dp<C>(wk[i].st.P,prm)){
     std::string key((const char*)wk[i].st.P.x.v,sizeof(uint32_t)*F::N);
     auto it=seen.find(key);
     if(it!=seen.end()){
      unsigned long long db=zn.sub(it->second.b,wk[i].b);
      if(db!=0){
       unsigned long long da=zn.sub(wk[i].a,it->second.a);
       unsigned long long k=zn.mul(da,zn.inv(db));
       uint32_t kl[F::N]={0}; kl[0]=(uint32_t)k; if(F::N>1) kl[1]=(uint32_t)(k>>32);
       pt test=C::scalar_mul(G,kl,prm.order_bits);
       if(!test.inf && F::eq(test.x,Q.x) && F::eq(test.y,Q.y)){ k_out=k; if(steps_out)*steps_out=total; return true; }
      }
     } else seen[key]=DP{wk[i].a,wk[i].b};
     restart[i]++; wsteps[i]=0; seed(wk[i],i,restart[i],G,Q,prm,zn);
    } else if(++wsteps[i]>=prm.max_steps){
     restart[i]++; wsteps[i]=0; seed(wk[i],i,restart[i],G,Q,prm,zn);
    }
   }
  }
  if(steps_out)*steps_out=total;
  return false;
 }
};

#endif
