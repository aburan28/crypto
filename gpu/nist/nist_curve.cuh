/* Affine short-Weierstrass points and the small bits of curve arithmetic the
 * Pollard-rho walk needs, generic over the field.  The hot path never
 * inverts: the walk batches denominators and calls affine_add_with_inv with
 * a precomputed 1/(x2-x1).  Seeding and the jump-table build take an
 * inversion-free a=-3 Jacobian ladder (n_double/n_madd) with one final
 * to_affine, instead of the affine double-and-add that paid a field
 * inversion per scalar bit; a general-a curve (the toy) keeps the affine
 * ladder.  The curve coefficient a and which ladder to use are supplied by
 * an AP policy: AMinus3 for the NIST curves, a literal for the toy.
 */
#ifndef GPU_NIST_CURVE_CUH
#define GPU_NIST_CURVE_CUH
#include "nist_fields.cuh"

template<class F> struct nap { typename F::elt x, y; uint32_t inf; };

/* a = -3, computed in whatever field F is (p256/p384 and any a=-3 toy).
 * is_a_minus3() lets the seeding ladder take the inversion-free Jacobian
 * path (n_double/n_madd are a=-3 formulas); a general-a policy returns
 * false and seeds with the affine ladder instead. */
template<class F> struct AMinus3 {
 NF_HD static typename F::elt a(){
  typename F::elt three=F::add(F::one(),F::add(F::one(),F::one()));
  return F::neg(three);
 }
 NF_HD static bool is_a_minus3(){return true;}
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

 /* Jacobian (X:Y:Z) -> affine, one field inversion.  Z=0 is the identity. */
 NF_HD static pt to_affine(const njac<F>&J){
  if(F::is_zero(J.Z)) return infinity();
  elt zi=F::inv(J.Z), zi2=F::sqr(zi), zi3=F::mul(zi2,zi);
  pt r; r.x=F::mul(J.X,zi2); r.y=F::mul(J.Y,zi3); r.inf=0; return r;
 }

 /* Complete mixed addition A(Jacobian, Z!=0) + B(affine).  n_madd alone is
  * the generic-position formula and returns garbage with Z=0 when A == +-B;
  * since Q is a multiple of G, the seeding accumulator coincides with +-P or
  * +-Q often enough to matter, so dispatch the equal (double) and inverse
  * (identity) cases explicitly.  Sets out_inf when the sum is the identity. */
 NF_HD static njac<F> jmadd(const njac<F>&A,const naff<F>&B,bool&out_inf){
  elt zz=F::sqr(A.Z);
  elt u2=F::mul(B.x,zz);                 /* B.x * Z^2 */
  elt s2=F::mul(F::mul(B.y,A.Z),zz);     /* B.y * Z^3 */
  if(F::eq(u2,A.X)){
   if(F::eq(s2,A.Y)){ out_inf=false; return n_double<F>(A); }  /* A == B  -> 2A */
   out_inf=true; return A;                                     /* A == -B -> O */
  }
  out_inf=false;
  return n_madd<F>(A,B);
 }

 /* [ka]P + [kb]Q by an interleaved a=-3 Jacobian ladder: no per-step field
  * inversion, one final to_affine.  `inf` tracks the identity (including a
  * mid-ladder cancellation, where a mixed add yields Z=0) so n_madd is only
  * ever fed a genuine point.  Used only when the curve is a=-3. */
 NF_HD static pt double_scalar_mul_jac(const pt&P,const uint32_t*ka,const pt&Q,const uint32_t*kb,int kbits){
  naff<F> Pa{P.x,P.y}, Qa{Q.x,Q.y};
  njac<F> acc{}; bool inf=true;
  for(int i=kbits-1;i>=0;i--){
   if(!inf) acc=n_double<F>(acc);
   if((ka[i>>5]>>(i&31))&1u){
    if(inf){ acc={P.x,P.y,F::one()}; inf=false; } else acc=jmadd(acc,Pa,inf);
   }
   if((kb[i>>5]>>(i&31))&1u){
    if(inf){ acc={Q.x,Q.y,F::one()}; inf=false; } else acc=jmadd(acc,Qa,inf);
   }
  }
  return inf ? infinity() : to_affine(acc);
 }

 /* [k]P by the Jacobian a=-3 ladder (single scalar). */
 NF_HD static pt scalar_mul_jac(const pt&P,const uint32_t*k,int kbits){
  naff<F> Pa{P.x,P.y};
  njac<F> acc{}; bool inf=true;
  for(int i=kbits-1;i>=0;i--){
   if(!inf) acc=n_double<F>(acc);
   if((k[i>>5]>>(i&31))&1u){
    if(inf){ acc={P.x,P.y,F::one()}; inf=false; } else acc=jmadd(acc,Pa,inf);
   }
  }
  return inf ? infinity() : to_affine(acc);
 }

 /* [k]P.  a=-3 curves take the inversion-free Jacobian ladder; a general-a
  * curve (the toy) falls back to affine double-and-add.  `k` is a
  * little-endian limb array of kbits bits.  Called by seeding and the
  * jump-table build, not the hot step. */
 NF_HD static pt scalar_mul(const pt&P,const uint32_t*k,int kbits){
  if(AP::is_a_minus3()) return scalar_mul_jac(P,k,kbits);
  pt acc=infinity();
  for(int i=kbits-1;i>=0;i--){
   acc=add(acc,acc);
   if((k[i>>5]>>(i&31))&1u) acc=add(acc,P);
  }
  return acc;
 }
 NF_HD static pt double_scalar_mul(const pt&P,const uint32_t*ka,const pt&Q,const uint32_t*kb,int kbits){
  if(AP::is_a_minus3()) return double_scalar_mul_jac(P,ka,Q,kb,kbits);
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
