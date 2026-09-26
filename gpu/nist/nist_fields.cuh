/* P-256/P-384 Montgomery arithmetic for NVIDIA GPUs.
 * 32-bit limbs, little endian. Both primes have p0=0xffffffff so n'=1.
 * Modulus limbs are 0, 1, -1, or -2 in base 2^32, making the CIOS
 * reduction row multiply-free. Host and device compile the same source.
 */
#ifndef GPU_NIST_FIELDS_CUH
#define GPU_NIST_FIELDS_CUH
#include <stdint.h>
#ifdef __CUDACC__
#define NF_HD __host__ __device__ __forceinline__
#else
#define NF_HD inline
#endif
template<int N> struct nfe { uint32_t v[N]; };

template<int N> NF_HD uint32_t n_add(uint32_t*r,const uint32_t*a,const uint32_t*b){
 uint64_t c=0; for(int i=0;i<N;i++){c+=(uint64_t)a[i]+b[i];r[i]=(uint32_t)c;c>>=32;}return (uint32_t)c;
}
template<int N> NF_HD uint32_t n_sub(uint32_t*r,const uint32_t*a,const uint32_t*b){
 uint64_t w=0;for(int i=0;i<N;i++){uint64_t d=(uint64_t)a[i]-b[i]-w;r[i]=(uint32_t)d;w=d>>63;}return (uint32_t)w;
}
template<int N> NF_HD void n_cmov(uint32_t*r,const uint32_t*a,uint32_t f){
 uint32_t m=0u-(f&1u);for(int i=0;i<N;i++)r[i]=(r[i]&~m)|(a[i]&m);
}
template<int N> NF_HD int n_eq(const uint32_t*a,const uint32_t*b){
 uint32_t x=0;for(int i=0;i<N;i++)x|=a[i]^b[i];return x==0;
}

/* Product row. CUDA specializations live in nist_ptx.cuh. */
template<int N> NF_HD void n_mac_row(uint32_t*t,const uint32_t*a,uint32_t b){
 uint64_t c=0;
#pragma unroll
 for(int j=0;j<N;j++){c+=(uint64_t)t[j]+(uint64_t)a[j]*b;t[j]=(uint32_t)c;c>>=32;}
 c+=t[N];t[N]=(uint32_t)c;t[N+1]+=(uint32_t)(c>>32);
}

template<class M> NF_HD uint32_t n_reduce_row_portable(uint32_t*t,uint32_t m){
 uint64_t c=0;
#pragma unroll
 for(int j=0;j<M::N;j++){
  int d=M::digit(j);uint64_t u;
  if(d==-1){u=(uint64_t)t[j]+c+(1ull<<32)-m;t[j]=(uint32_t)u;c=(uint64_t)m+(u>>32)-1;}
  else if(d==-2){u=(uint64_t)t[j]+c+(2ull<<32)-2ull*m;t[j]=(uint32_t)u;c=(uint64_t)m+(u>>32)-2;}
  else if(d==0){u=(uint64_t)t[j]+c;t[j]=(uint32_t)u;c=u>>32;}
  else{u=(uint64_t)t[j]+c+m;t[j]=(uint32_t)u;c=u>>32;}
 }
 return (uint32_t)c;
}

template<class M> NF_HD uint32_t n_reduce_row(uint32_t*t,uint32_t m){return n_reduce_row_portable<M>(t,m);}

struct P256Mod {
 static const int N=8;
 NF_HD static uint32_t p(int i){const uint32_t x[8]={0xffffffffu,0xffffffffu,0xffffffffu,0,0,0,1,0xffffffffu};return x[i];}
 NF_HD static uint32_t r1(int i){const uint32_t x[8]={1,0,0,0xffffffffu,0xffffffffu,0xffffffffu,0xfffffffeu,0};return x[i];}
 NF_HD static uint32_t r2(int i){const uint32_t x[8]={3,0,0xffffffffu,0xfffffffbu,0xfffffffeu,0xffffffffu,0xfffffffdu,4};return x[i];}
 NF_HD static int digit(int i){const int8_t x[8]={-1,-1,-1,0,0,0,1,-1};return x[i];}
};
struct P384Mod {
 static const int N=12;
 NF_HD static uint32_t p(int i){const uint32_t x[12]={0xffffffffu,0,0,0xffffffffu,0xfffffffeu,0xffffffffu,0xffffffffu,0xffffffffu,0xffffffffu,0xffffffffu,0xffffffffu,0xffffffffu};return x[i];}
 NF_HD static uint32_t r1(int i){const uint32_t x[12]={1,0xffffffffu,0xffffffffu,0,1,0,0,0,0,0,0,0};return x[i];}
 NF_HD static uint32_t r2(int i){const uint32_t x[12]={1,0xfffffffeu,0,2,0,0xfffffffeu,0,2,1,0,0,0};return x[i];}
 NF_HD static int digit(int i){const int8_t x[12]={-1,0,0,-1,-2,-1,-1,-1,-1,-1,-1,-1};return x[i];}
};

template<class M> struct NistField {
 static const int N=M::N; typedef nfe<N> elt;
 NF_HD static elt zero(){elt r{};return r;}
 NF_HD static elt modulus(){elt r;for(int i=0;i<N;i++)r.v[i]=M::p(i);return r;}
 NF_HD static elt one(){elt r;for(int i=0;i<N;i++)r.v[i]=M::r1(i);return r;}
 NF_HD static int eq(const elt&a,const elt&b){return n_eq<N>(a.v,b.v);}
 NF_HD static elt add(const elt&a,const elt&b){
  elt s,t,p=modulus();uint32_t c=n_add<N>(s.v,a.v,b.v),w=n_sub<N>(t.v,s.v,p.v);
  n_cmov<N>(s.v,t.v,c|(w^1u));return s;
 }
 NF_HD static elt sub(const elt&a,const elt&b){
  elt d,t,p=modulus();uint32_t w=n_sub<N>(d.v,a.v,b.v);n_add<N>(t.v,d.v,p.v);n_cmov<N>(d.v,t.v,w);return d;
 }
 NF_HD static elt neg(const elt&a){return sub(zero(),a);}
 NF_HD static elt dbl(const elt&a){return add(a,a);}
 NF_HD static elt mul(const elt&a,const elt&b){
  uint32_t t[N+2];
#pragma unroll
  for(int j=0;j<N+2;j++)t[j]=0;
#pragma unroll
  for(int i=0;i<N;i++){
   n_mac_row<N>(t,a.v,b.v[i]);
   uint32_t m=t[0];uint32_t rc=n_reduce_row<M>(t,m);
   uint64_t z=(uint64_t)t[N]+rc;t[N]=(uint32_t)z;t[N+1]+=(uint32_t)(z>>32);
#pragma unroll
   for(int j=0;j<N+1;j++)t[j]=t[j+1];
   t[N+1]=0;
  }
  elt r,s,p=modulus();for(int j=0;j<N;j++)r.v[j]=t[j];
  uint32_t w=n_sub<N>(s.v,r.v,p.v);n_cmov<N>(r.v,s.v,(t[N]!=0)|(w^1u));return r;
 }
 NF_HD static elt sqr(const elt&a){return mul(a,a);}
 NF_HD static elt from_canonical(const elt&a){elt r2;for(int i=0;i<N;i++)r2.v[i]=M::r2(i);return mul(a,r2);}
 NF_HD static elt to_canonical(const elt&a){elt c=zero();c.v[0]=1;return mul(a,c);}
};

typedef NistField<P256Mod> Fp256;
typedef NistField<P384Mod> Fp384;

/* a=-3 Jacobian formulas: no inversion in the hot point path. */
template<class F> struct njac { typename F::elt X,Y,Z; };
template<class F> struct naff { typename F::elt x,y; };
template<class F> NF_HD njac<F> n_double(const njac<F>&P){
 typedef typename F::elt E;
 E delta=F::sqr(P.Z),gamma=F::sqr(P.Y),beta=F::mul(P.X,gamma);
 E xm=F::sub(P.X,delta),xp=F::add(P.X,delta),alpha=F::mul(xm,xp);
 alpha=F::add(alpha,F::dbl(alpha));
 E X3=F::sub(F::sqr(alpha),F::dbl(F::dbl(beta)));X3=F::sub(X3,F::dbl(F::dbl(beta)));
 E gz=F::add(P.Y,P.Z);E Z3=F::sub(F::sub(F::sqr(gz),gamma),delta);
 E g2=F::sqr(gamma),fourb=F::dbl(F::dbl(beta));
 E Y3=F::sub(F::mul(alpha,F::sub(fourb,X3)),F::dbl(F::dbl(F::dbl(g2))));
 return {X3,Y3,Z3};
}
template<class F> NF_HD njac<F> n_madd(const njac<F>&P,const naff<F>&Q){
 typedef typename F::elt E;
 E zz=F::sqr(P.Z),u2=F::mul(Q.x,zz),s2=F::mul(Q.y,F::mul(P.Z,zz));
 E h=F::sub(u2,P.X),hh=F::sqr(h),i=F::dbl(F::dbl(hh)),j=F::mul(h,i);
 E r=F::dbl(F::sub(s2,P.Y)),v=F::mul(P.X,i),X3=F::sub(F::sub(F::sqr(r),j),F::dbl(v));
 E Y3=F::sub(F::mul(r,F::sub(v,X3)),F::dbl(F::mul(P.Y,j)));
 E zh=F::add(P.Z,h),Z3=F::sub(F::sub(F::sqr(zh),zz),hh);
 return {X3,Y3,Z3};
}
#endif
