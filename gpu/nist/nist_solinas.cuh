/* Direct generalized-Mersenne candidates for comparison with sparse CIOS.
 * Product uses the same PTX row primitive; reduction folds signed radix-2^32
 * coefficients. Two normalize/fold rounds suffice for a full 2N-limb input.
 */
#ifndef GPU_NIST_SOLINAS_CUH
#define GPU_NIST_SOLINAS_CUH
#include "nist_ptx.cuh"
struct P256SolMod:P256Mod{
 NF_HD static int roff(int i){const int8_t x[4]={7,6,3,0};return x[i];}
 NF_HD static int rcoef(int i){const int8_t x[4]={1,-1,-1,1};return x[i];}
};
struct P384SolMod:P384Mod{
 NF_HD static int roff(int i){const int8_t x[4]={4,3,1,0};return x[i];}
 NF_HD static int rcoef(int i){const int8_t x[4]={1,1,-1,1};return x[i];}
};
template<class M> struct SolinasField{
 static const int N=M::N;typedef nfe<N> elt;
 NF_HD static elt zero(){elt r{};return r;}NF_HD static elt one(){elt r{};r.v[0]=1;return r;}
 NF_HD static elt modulus(){elt r;for(int i=0;i<N;i++)r.v[i]=M::p(i);return r;}
 NF_HD static int eq(const elt&a,const elt&b){return n_eq<N>(a.v,b.v);}
 NF_HD static elt add(const elt&a,const elt&b){elt s,t,p=modulus();uint32_t c=n_add<N>(s.v,a.v,b.v),w=n_sub<N>(t.v,s.v,p.v);n_cmov<N>(s.v,t.v,c|(w^1u));return s;}
 NF_HD static elt sub(const elt&a,const elt&b){elt d,t,p=modulus();uint32_t w=n_sub<N>(d.v,a.v,b.v);n_add<N>(t.v,d.v,p.v);n_cmov<N>(d.v,t.v,w);return d;}
 NF_HD static elt dbl(const elt&a){return add(a,a);}NF_HD static elt neg(const elt&a){return sub(zero(),a);}
 NF_HD static elt mul(const elt&a,const elt&b){
  uint32_t w[2*N+2];
#pragma unroll
  for(int i=0;i<2*N+2;i++)w[i]=0;
#pragma unroll
  for(int i=0;i<N;i++)n_mac_row<N>(w+i,a.v,b.v[i]);
  int64_t x[2*N+1];
#pragma unroll
  for(int i=0;i<2*N+1;i++)x[i]=(int64_t)w[i];
#pragma unroll
  for(int k=2*N-1;k>=N;k--){int64_t h=x[k];x[k]=0;
#pragma unroll
   for(int q=0;q<4;q++)x[k-N+M::roff(q)]+=(int64_t)M::rcoef(q)*h;
  }
  int64_t y[N+1];for(int i=0;i<N;i++)y[i]=x[i];y[N]=0;
#pragma unroll
  for(int round=0;round<2;round++){
#pragma unroll
   for(int i=0;i<N;i++){int64_t q=y[i]>>32;y[i]-=q*(1ll<<32);y[i+1]+=q;}
   int64_t h=y[N];y[N]=0;
#pragma unroll
   for(int q=0;q<4;q++)y[M::roff(q)]+=(int64_t)M::rcoef(q)*h;
  }
#pragma unroll
  for(int i=0;i<N;i++){int64_t q=y[i]>>32;y[i]-=q*(1ll<<32);if(i+1<N)y[i+1]+=q;}
  elt r,s,p=modulus();for(int i=0;i<N;i++)r.v[i]=(uint32_t)y[i];
  uint32_t bw=n_sub<N>(s.v,r.v,p.v);n_cmov<N>(r.v,s.v,bw^1u);return r;
 }
 NF_HD static elt sqr(const elt&a){return mul(a,a);}
 NF_HD static elt from_canonical(const elt&a){return a;}NF_HD static elt to_canonical(const elt&a){return a;}
};
typedef SolinasField<P256SolMod> Fp256Sol;
typedef SolinasField<P384SolMod> Fp384Sol;
#endif
