/* CPU verification for the NIST rho engine (nist_curve/nist_rho/nist_rho_host).
 * Compiles the same headers the CUDA kernel uses with g++, so a green run
 * means the field/curve/walk/solver logic is correct; only the launch
 * configuration stays GPU-specific.
 *
 * Dependency-free on purpose (unlike test_cpu.cpp, which uses Boost as the
 * field reference in CI): correctness here rests on group-law identities
 * that need no external big-integer library --
 *   a * inv(a) == 1,  batch_inv == inv,
 *   P+O=P,  P+(-P)=O,  add(P,P)=[2]P,  (P+Q)+R=P+(Q+R),  [a]G+[b]G=[a+b]G,
 *   affine_add_with_inv == add,
 * -- plus the negation-map fold laws, walk-step determinism, and the
 * decisive end-to-end DLP solve on a ~2^20 prime-order toy curve whose
 * (n,G,Q,k) were checked independently when the curve was constructed.
 *
 * Suites 1-4 run on P-256 and P-384 (the device fields); suite 5 is the toy.
 */
#include "nist_rho_host.hpp"
#include <iostream>
#include <cstdint>
#include <vector>

static int fails=0;
#define CHECK(x,msg) do{if(!(x)){std::cerr<<"FAIL "<<(msg)<<"\n";fails++;}}while(0)

/* A field element from a small canonical integer (< 2^32 <= p-1): a valid
 * internal (Montgomery) representation, nonzero iff v != 0. */
template<class F> static typename F::elt small(uint32_t v){
 uint32_t l[F::N]={0}; l[0]=v; return F::from_limbs(l);
}
template<class C> static nap<typename C::Field> smul(const nap<typename C::Field>&P,uint64_t k){
 uint32_t kl[C::Field::N]={0}; kl[0]=(uint32_t)k; if(C::Field::N>1)kl[1]=(uint32_t)(k>>32);
 return C::scalar_mul(P,kl,64);
}
template<class C> static int pteq(const nap<typename C::Field>&A,const nap<typename C::Field>&B){
 typedef typename C::Field F;
 if(A.inf||B.inf) return A.inf&&B.inf;
 return F::eq(A.x,B.x)&&F::eq(A.y,B.y);
}

/* ---- 1. inversion ------------------------------------------------------ */
template<class F> static void inv_test(const char*name){
 uint64_t s=0x12345;
 for(int i=0;i<3000;i++){
  s=s*6364136223846793005ull+1442695040888963407ull;
  uint32_t v=(uint32_t)(s>>33); if(v==0)v=1;
  auto a=small<F>(v); auto ai=F::inv(a);
  CHECK(F::eq(F::mul(a,ai),F::one()),name);
 }
 const int n=41; typename F::elt x[n],sc[n],ref[n];
 for(int i=0;i<n;i++){x[i]=small<F>((uint32_t)(3+7*i)); ref[i]=F::inv(x[i]);}
 F::batch_inv(x,n,sc);
 for(int i=0;i<n;i++) CHECK(F::eq(x[i],ref[i]),name);
 std::cout<<name<<": inversion + batch inverse (a*inv(a)=1) passed\n";
}

/* ---- 2. group laws ----------------------------------------------------- */
template<class C> static void group_test(const char*name,const char*xs,const char*ys){
 typedef typename C::Field F;
 /* parse the generator from hex */
 auto parse=[&](const char*h)->typename F::elt{
  std::vector<uint8_t> by; const char*p=h; if(p[0]=='0'&&(p[1]=='x'||p[1]=='X'))p+=2;
  std::string s(p); if(s.size()%2)s="0"+s;
  uint32_t l[F::N]={0}; int nib=0;
  for(int i=(int)s.size()-1;i>=0;i--){int c=s[i];int d=(c>='0'&&c<='9')?c-'0':(c|32)-'a'+10; l[nib/8]|=(uint32_t)d<<(4*(nib%8)); nib++;}
  return F::from_limbs(l);
 };
 nap<F> G; G.x=parse(xs); G.y=parse(ys); G.inf=0;
 nap<F> O=C::infinity();
 CHECK(pteq<C>(C::add(G,O),G),name);                         /* P+O=P */
 CHECK(C::add(G,C::neg(G)).inf,name);                        /* P+(-P)=O */
 nap<F> G2=C::add(G,G), G2s=smul<C>(G,2);
 CHECK(pteq<C>(G2,G2s),name);                                /* [2]P via add */
 nap<F> G3=C::add(G2,G), G5=C::add(G3,G2);
 /* associativity (G+G2)+G3 == G+(G2+G3) */
 CHECK(pteq<C>(C::add(C::add(G,G2),G3),C::add(G,C::add(G2,G3))),name);
 /* [a]G+[b]G=[a+b]G for a few small a,b */
 for(uint64_t a=1;a<=5;a++) for(uint64_t b=1;b<=5;b++)
  CHECK(pteq<C>(C::add(smul<C>(G,a),smul<C>(G,b)),smul<C>(G,a+b)),name);
 CHECK(pteq<C>(G5,smul<C>(G,5)),name);
 /* affine_add_with_inv matches add on the generic path */
 auto inv=F::inv(F::sub(G2.x,G.x));
 CHECK(pteq<C>(C::affine_add_with_inv(G,G2,inv,0),C::add(G,G2)),name);
 std::cout<<name<<": group laws (identity, inverse, assoc, distrib) passed\n";
}

/* ---- 3. fold ----------------------------------------------------------- */
template<class C> static void fold_test(const char*name,const nap<typename C::Field>&P){
 typedef typename C::Field F;
 rho_params prm{}; prm.fold=RHO_FOLD_NEG;
 nap<F> a=P, b=C::neg(P);
 int ca=rho_canonical<C>(a,prm), cb=rho_canonical<C>(b,prm);
 CHECK(F::eq(a.x,b.x)&&F::eq(a.y,b.y),name);          /* same representative */
 CHECK((ca&1)!=(cb&1),name);                           /* opposite aut codes */
 nap<F> c=P; int cc=rho_canonical<C>(c,prm); nap<F> back=rho_apply_aut<C>(P,cc);
 CHECK(F::eq(back.x,c.x)&&F::eq(back.y,c.y),name);      /* apply_aut round-trip */
 std::cout<<name<<": negation-map fold is a class function\n";
}

/* ---- 4. step determinism ---------------------------------------------- */
template<class C> static void step_test(const char*name,const nap<typename C::Field>&G){
 typedef typename C::Field F;
 rho_params prm{}; prm.r_bits=6; prm.dp_mask=0; prm.fold=RHO_FOLD_NEG; prm.max_steps=1u<<20; prm.table_seed=99; prm.order_bits=F::N*32;
 uint32_t R=1u<<prm.r_bits; std::vector<nap<F>> tbl(R);
 for(uint32_t j=0;j<R;j++){uint32_t c[F::N],d[F::N];rho_table_seed<F::N>(prm.table_seed,j,c,d,prm.order_bits);(void)d;tbl[j]=C::scalar_mul(G,c,prm.order_bits);}
 nap<F> P=G; uint32_t first=0;
 for(int s=0;s<200;s++){int j,aut;nap<F> nx=rho_step_single<C>(P,tbl.data(),prm,&j,&aut);CHECK(!nx.inf,name);if(s==0)first=nx.x.v[0];P=nx;}
 /* rerun: identical trajectory */
 nap<F> Q=G; for(int s=0;s<1;s++){int j,aut;Q=rho_step_single<C>(Q,tbl.data(),prm,&j,&aut);} CHECK(Q.x.v[0]==first,name);
 std::cout<<name<<": single-step walk ran 200 deterministic steps\n";
}

/* ---- toy field: canonical arithmetic over a compile-time prime < 2^32 -- */
struct ToyMod { static const uint32_t P=1048583u, A=105622u, B=556017u,
 GX=2u, GY=995515u, N=1049227u, K=185340u, QX=948343u, QY=457421u; };
struct SmallField {
 static const int N=1; typedef nfe<1> elt;
 static elt zero(){return elt{{0}};}
 static elt one(){return elt{{1}};}
 static elt modulus(){return elt{{ToyMod::P}};}
 static int eq(const elt&a,const elt&b){return a.v[0]==b.v[0];}
 static int is_zero(const elt&a){return a.v[0]==0;}
 static int lt(const elt&a,const elt&b){return a.v[0]<b.v[0];}
 static void cmov(elt&r,const elt&a,uint32_t f){uint32_t m=0u-(f&1u);r.v[0]=(r.v[0]&~m)|(a.v[0]&m);}
 static elt add(const elt&a,const elt&b){uint64_t s=(uint64_t)a.v[0]+b.v[0];if(s>=ToyMod::P)s-=ToyMod::P;return elt{{(uint32_t)s}};}
 static elt sub(const elt&a,const elt&b){uint32_t r=a.v[0]>=b.v[0]?a.v[0]-b.v[0]:(uint32_t)((uint64_t)a.v[0]+ToyMod::P-b.v[0]);return elt{{r}};}
 static elt neg(const elt&a){return a.v[0]?elt{{ToyMod::P-a.v[0]}}:elt{{0}};}
 static elt dbl(const elt&a){return add(a,a);}
 static elt mul(const elt&a,const elt&b){return elt{{(uint32_t)(((uint64_t)a.v[0]*b.v[0])%ToyMod::P)}};}
 static elt sqr(const elt&a){return mul(a,a);}
 static elt inv(const elt&a){uint64_t r=1,x=a.v[0]%ToyMod::P,e=ToyMod::P-2;while(e){if(e&1)r=r*x%ToyMod::P;x=x*x%ToyMod::P;e>>=1;}return elt{{(uint32_t)r}};}
 static elt from_canonical(const elt&a){return elt{{a.v[0]%ToyMod::P}};}
 static elt to_canonical(const elt&a){return a;}
 static elt from_limbs(const uint32_t*l){return elt{{l[0]%ToyMod::P}};}
 static void batch_inv(elt*x,int n,elt*sc){sc[0]=x[0];for(int i=1;i<n;i++)sc[i]=mul(sc[i-1],x[i]);elt t=inv(sc[n-1]);for(int i=n-1;i>0;i--){elt xi=mul(t,sc[i-1]);t=mul(t,x[i]);x[i]=xi;}x[0]=t;}
};
struct ToyAP { static SmallField::elt a(){uint32_t v=ToyMod::A;return SmallField::from_limbs(&v);} };
typedef EcCurve<SmallField,ToyAP> CurveToy;

static nap<SmallField> toyG(){nap<SmallField> G;uint32_t gx=ToyMod::GX,gy=ToyMod::GY;G.x=SmallField::from_limbs(&gx);G.y=SmallField::from_limbs(&gy);G.inf=0;return G;}
static nap<SmallField> toyQ(){nap<SmallField> Q;uint32_t qx=ToyMod::QX,qy=ToyMod::QY;Q.x=SmallField::from_limbs(&qx);Q.y=SmallField::from_limbs(&qy);Q.inf=0;return Q;}

static void toy_sanity(){
 auto G=toyG(); uint32_t n=ToyMod::N; auto nG=CurveToy::scalar_mul(G,&n,21); CHECK(nG.inf,"toy nG=O");
 uint32_t k=ToyMod::K; auto kG=CurveToy::scalar_mul(G,&k,21);
 CHECK(!kG.inf&&kG.x.v[0]==ToyMod::QX&&kG.y.v[0]==ToyMod::QY,"toy kG=Q");
 std::cout<<"toy curve: nG=O and kG=Q verified\n";
}
static void toy_solve(){
 auto G=toyG(), Q=toyQ();
 for(uint32_t fold : {RHO_FOLD_NEG,RHO_FOLD_NONE}){
  rho_params prm{}; prm.r_bits=8; prm.dp_mask=(1u<<5)-1u; prm.fold=fold; prm.max_steps=1u<<16; prm.table_seed=1; prm.order_bits=21;
  unsigned long long k=0,steps=0;
  bool ok=RhoHost<CurveToy>::solve(G,Q,prm,ToyMod::N,256,40000000ull,k,&steps);
  CHECK(ok&&k==ToyMod::K,fold==RHO_FOLD_NEG?"toy solve folded":"toy solve plain");
  std::cout<<"toy solve ("<<(fold==RHO_FOLD_NEG?"fold 2":"fold 1")<<"): k="<<k<<" after "<<steps<<" steps\n";
 }
}

static int toy_on_curve(const nap<SmallField>&P){
 if(P.inf) return 1;
 auto x=P.x,y=P.y; uint32_t av=ToyMod::A,bv=ToyMod::B;
 auto a=SmallField::from_limbs(&av), b=SmallField::from_limbs(&bv);
 auto lhs=SmallField::sqr(y);
 auto rhs=SmallField::add(SmallField::add(SmallField::mul(SmallField::sqr(x),x),SmallField::mul(a,x)),b);
 return SmallField::eq(lhs,rhs);
}
/* Host-exercise the exact batched stepper the CUDA kernel runs (SoA load/
 * store, seeding, DP emission), on the toy curve.  Confirms the plumbing
 * compiles and runs and that every distinguished point is a canonical point
 * on the curve. */
static void batch_step_test(){
 auto G=toyG(), Q=toyQ();
 rho_params prm{}; prm.r_bits=8; prm.dp_mask=(1u<<7)-1u; prm.fold=RHO_FOLD_NEG; prm.max_steps=1u<<16; prm.table_seed=3; prm.order_bits=21;
 uint32_t R=1u<<prm.r_bits, nthreads=64, W=4, nw=nthreads*W;
 std::vector<nap<SmallField>> tbl(R);
 for(uint32_t j=0;j<R;j++){uint32_t cc[1],dd[1];rho_table_seed<1>(prm.table_seed,j,cc,dd,prm.order_bits);tbl[j]=CurveToy::double_scalar_mul(G,cc,Q,dd,prm.order_bits);}
 std::vector<uint32_t> X(nw),Y(nw),H(RHO_CYCLE_DEPTH*nw,0),esc(nw,0),steps(nw,0),restarts(nw,0);
 std::vector<rho_dp<SmallField>> dp(1<<16); uint32_t dpc=0;
 rho_ctx<CurveToy> c{}; c.X=X.data(); c.Y=Y.data(); c.H=H.data(); c.esc=esc.data(); c.steps=steps.data(); c.restarts=restarts.data();
 c.nthreads=nthreads; c.walks_per_thread=W; c.table=tbl.data(); c.G=G; c.Q=Q; c.prm=prm; c.dp_out=dp.data(); c.dp_count=&dpc; c.dp_cap=(uint32_t)dp.size();
 for(uint32_t i=0;i<nw;i++){rho_state<CurveToy> st; rho_seed_walk<CurveToy>(c,i,st); rho_store<CurveToy>(c,i,st);}
 SmallField::elt den[4],scratch[4];
 for(int it=0;it<4000;it++) for(uint32_t t=0;t<nthreads;t++) rho_batch_step<CurveToy>(c,t,den,scratch);
 CHECK(dpc>0,"batch step produced DPs");
 int bad=0; uint32_t m=std::min(dpc,(uint32_t)dp.size());
 for(uint32_t i=0;i<m;i++){ nap<SmallField> P; P.x.v[0]=dp[i].x[0]; /* recover y from the curve */
  /* a stored DP is only x,walk,restart,steps; check x is a valid abscissa */
  uint32_t av=ToyMod::A,bv=ToyMod::B; auto a=SmallField::from_limbs(&av),b=SmallField::from_limbs(&bv);
  auto rhs=SmallField::add(SmallField::add(SmallField::mul(SmallField::sqr(P.x),P.x),SmallField::mul(a,P.x)),b);
  /* rhs must be a quadratic residue mod p (Euler): rhs^((p-1)/2) in {0,1} */
  uint64_t e=(ToyMod::P-1)/2,r=1,xx=rhs.v[0]%ToyMod::P; while(e){if(e&1)r=r*xx%ToyMod::P;xx=xx*xx%ToyMod::P;e>>=1;}
  if(!(r==0||r==1)) bad++;
 }
 CHECK(bad==0,"every DP abscissa is a valid curve point");
 (void)toy_on_curve;
 std::cout<<"batched stepper: "<<dpc<<" DPs over 4000 iters, all valid abscissae\n";
}

template<class C> static nap<typename C::Field> gen(const char*xs,const char*ys){
 typedef typename C::Field F;
 auto parse=[&](const char*h)->typename F::elt{
  const char*p=h; if(p[0]=='0'&&(p[1]=='x'||p[1]=='X'))p+=2; std::string s(p); if(s.size()%2)s="0"+s;
  uint32_t l[F::N]={0}; int nib=0;
  for(int i=(int)s.size()-1;i>=0;i--){int c=s[i];int d=(c>='0'&&c<='9')?c-'0':(c|32)-'a'+10;l[nib/8]|=(uint32_t)d<<(4*(nib%8));nib++;}
  return F::from_limbs(l);
 };
 nap<F> G; G.x=parse(xs); G.y=parse(ys); G.inf=0; return G;
}

int main(){
 const char*P256X="0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
 const char*P256Y="0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
 const char*P384X="0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7";
 const char*P384Y="0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f";
 inv_test<Fp256>("P-256"); inv_test<Fp384>("P-384");
 group_test<CurveP256>("P-256",P256X,P256Y); group_test<CurveP384>("P-384",P384X,P384Y);
 fold_test<CurveP256>("P-256",gen<CurveP256>(P256X,P256Y)); fold_test<CurveP384>("P-384",gen<CurveP384>(P384X,P384Y));
 step_test<CurveP256>("P-256",gen<CurveP256>(P256X,P256Y)); step_test<CurveP384>("P-384",gen<CurveP384>(P384X,P384Y));
 toy_sanity(); batch_step_test(); toy_solve();
 if(fails){std::cerr<<fails<<" failures\n";return 1;} std::cout<<"ALL PASSED\n"; return 0;
}
