#include "nist_curve.cuh"
#include <cstdio>
#include <chrono>
template<class C> static nap<typename C::Field> parse(const char*xs,const char*ys){using F=typename C::Field;auto p=[&](const char*h){const char*q=h;if(q[0]=='0'&&q[1]=='x')q+=2;int n=0;while(q[n])n++;uint32_t l[F::N]={0};int nib=0;for(int i=n-1;i>=0;i--){int c=q[i];int d=(c>='0'&&c<='9')?c-'0':(c|32)-'a'+10;l[nib/8]|=(uint32_t)d<<(4*(nib%8));nib++;}return F::from_limbs(l);};nap<F> P;P.x=p(xs);P.y=p(ys);P.inf=0;return P;}
template<class C> static nap<typename C::Field> affdsm(const nap<typename C::Field>&P,const uint32_t*ka,const nap<typename C::Field>&Q,const uint32_t*kb,int kb_){nap<typename C::Field> a=C::infinity();for(int i=kb_-1;i>=0;i--){a=C::add(a,a);if((ka[i>>5]>>(i&31))&1u)a=C::add(a,P);if((kb[i>>5]>>(i&31))&1u)a=C::add(a,Q);}return a;}
template<class C> static void run(const char*nm,const nap<typename C::Field>&G,int bits){
 using F=typename C::Field; const int N=2000; int L=F::N;
 uint32_t k7[16]={7}; auto Q=C::scalar_mul(G,k7,8);
 static uint32_t A[16*2000],B[16*2000]; uint64_t s=0x1234;
 for(int i=0;i<L*N;i++){s=s*6364136223846793005ull+1;A[i]=(uint32_t)(s>>33);s=s*6364136223846793005ull+1;B[i]=(uint32_t)(s>>33);}
 volatile uint32_t sink=0; int bad=0;
 auto t0=std::chrono::steady_clock::now();
 for(int i=0;i<N;i++){auto R=C::double_scalar_mul(G,&A[L*i],Q,&B[L*i],bits);sink^=R.x.v[0];}
 auto t1=std::chrono::steady_clock::now();
 for(int i=0;i<N;i++){auto R=affdsm<C>(G,&A[L*i],Q,&B[L*i],bits);sink^=R.x.v[0];}
 auto t2=std::chrono::steady_clock::now();
 for(int i=0;i<N;i++){auto r1=C::double_scalar_mul(G,&A[L*i],Q,&B[L*i],bits);auto r2=affdsm<C>(G,&A[L*i],Q,&B[L*i],bits);if(!((r1.inf&&r2.inf)||(F::eq(r1.x,r2.x)&&F::eq(r1.y,r2.y))))bad++;}
 double jac=std::chrono::duration<double,std::milli>(t1-t0).count()/N, aff=std::chrono::duration<double,std::milli>(t2-t1).count()/N;
 printf("%-6s bits=%d  jacobian %.4f ms/seed  affine %.4f ms/seed  speedup %.2fx  mismatches=%d\n",nm,bits,jac,aff,aff/jac,bad);
 (void)sink;
}
int main(){
 run<CurveP256>("P-256",parse<CurveP256>("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296","0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5"),256);
 run<CurveP384>("P-384",parse<CurveP384>("0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7","0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f"),384);
 return 0;
}
