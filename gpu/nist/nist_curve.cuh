/* Affine short-Weierstrass points and the small bits of curve arithmetic the
 * Pollard-rho walk needs, generic over the field.  The hot path never
 * inverts: the walk batches denominators and calls affine_add_with_inv with
 * a precomputed 1/(x2-x1).  Seeding (rare) uses a self-inverting affine
 * double-and-add, so the whole walk layer needs only the field interface --
 * no Jacobian, no a=-3 assumption.  The curve coefficient a is supplied by
 * an AP policy: AMinus3 for the NIST curves, a literal for the toy.
 */
#ifndef GPU_NIST_CURVE_CUH
#define GPU_NIST_CURVE_CUH
#include "nist_fields.cuh"

template<class F> struct nap { typename F::elt x, y; uint32_t inf; };

/* a = -3, computed in whatever field F is (p256/p384 and any a=-3 toy). */
template<class F> struct AMinus3 {
 NF_HD static typename F::elt a(){
  typename F::elt three=F::add(F::one(),F::add(F::one(),F::one()));
  return F::neg(three);
 }
};

template<class F, class AP> struct EcCurve {
 typedef F Field; typedef typename F::elt elt; typedef nap<F> pt;
 NF_HD static elt a(){return AP::a();}
 NF_HD static pt infinity(){pt r; r.x=F::zero(); r.y=F::zero(); r.inf=1; return r;}
 NF_HD static int is_inf(const pt&P){return P.inf;}
 NF_HD static pt neg(const pt&P){pt r=P; if(!P.inf) r.y=F::neg(P.y); return r;}

 /* P+Q with a precomputed inv = 1/(Q.x-P.x), or 1/(2 P.y) when doubling.
  * The hot step; no field inversion here. */
 NF_HD static pt affine_add_with_inv(const pt&P,const pt&Q,const elt&inv,int doubling){
  elt num;
  if(doubling){ elt x2=F::sqr(P.x); num=F::add(F::add(x2,F::add(x2,x2)),a()); }  /* 3x^2 + a */
  else { num=F::sub(Q.y,P.y); }
  elt lam=F::mul(num,inv);
  pt r; r.inf=0;
  r.x=F::sub(F::sub(F::sqr(lam),P.x),Q.x);
  r.y=F::sub(F::mul(lam,F::sub(P.x,r.x)),P.y);
  return r;
 }

 /* P+Q computing its own inverse: used off the hot path (table build,
  * seeding, the host reference).  Handles identity and the vertical case. */
 NF_HD static pt add(const pt&P,const pt&Q){
  if(P.inf) return Q;
  if(Q.inf) return P;
  if(F::eq(P.x,Q.x)){
   if(F::eq(P.y,Q.y)){
    if(F::is_zero(P.y)) return infinity();
    return affine_add_with_inv(P,Q,F::inv(F::dbl(P.y)),1);
   }
   return infinity();                 /* Q = -P */
  }
  return affine_add_with_inv(P,Q,F::inv(F::sub(Q.x,P.x)),0);
 }

 /* [k]P by affine double-and-add; k is a little-endian limb array of kbits
  * bits.  Only seeding and the table build call it. */
 NF_HD static pt scalar_mul(const pt&P,const uint32_t*k,int kbits){
  pt acc=infinity();
  for(int i=kbits-1;i>=0;i--){
   acc=add(acc,acc);
   if((k[i>>5]>>(i&31))&1u) acc=add(acc,P);
  }
  return acc;
 }
 NF_HD static pt double_scalar_mul(const pt&P,const uint32_t*ka,const pt&Q,const uint32_t*kb,int kbits){
  pt acc=infinity();
  for(int i=kbits-1;i>=0;i--){
   acc=add(acc,acc);
   if((ka[i>>5]>>(i&31))&1u) acc=add(acc,P);
   if((kb[i>>5]>>(i&31))&1u) acc=add(acc,Q);
  }
  return acc;
 }
};

/* The NIST instantiations. */
typedef EcCurve<Fp256, AMinus3<Fp256> > CurveP256;
typedef EcCurve<Fp384, AMinus3<Fp384> > CurveP384;

#endif
