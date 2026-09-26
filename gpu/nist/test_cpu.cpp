#include "nist_solinas.cuh"
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <random>
using boost::multiprecision::cpp_int;

static int fails=0;
#define CHECK(x,msg) do{if(!(x)){std::cerr<<"FAIL "<<msg<<"\n";fails++;}}while(0)

template<int N> cpp_int limbs(const uint32_t *a){
 cpp_int x=0;for(int i=N-1;i>=0;i--){x<<=32;x+=a[i];}return x;
}
template<class F> typename F::elt raw(cpp_int x){
 typename F::elt r{};for(int i=0;i<F::N;i++){r.v[i]=(uint32_t)(x&0xffffffff);x>>=32;}return r;
}
template<class F> cpp_int canon(const typename F::elt&a){auto c=F::to_canonical(a);return limbs<F::N>(c.v);}
static cpp_int mod(cpp_int x,const cpp_int&p){x%=p;if(x<0)x+=p;return x;}

template<class F,class M> void field_test(const char *name,int rounds){
 cpp_int p=limbs<F::N>(F::modulus().v);
 std::mt19937_64 g(0x50323536);
 for(int k=0;k<rounds;k++){
  cpp_int a=0,b=0;for(int i=0;i<F::N/2;i++){a=(a<<64)+g();b=(b<<64)+g();}a%=p;b%=p;
  auto A=F::from_canonical(raw<F>(a)),B=F::from_canonical(raw<F>(b));
  CHECK(canon<F>(A)==a,name);CHECK(canon<F>(B)==b,name);
  CHECK(canon<F>(F::mul(A,B))==a*b%p,name);
  CHECK(canon<F>(F::sqr(A))==a*a%p,name);
  CHECK(canon<F>(F::add(A,B))==(a+b)%p,name);
  CHECK(canon<F>(F::sub(A,B))==mod(a-b,p),name);
 }
 std::cout<<name<<": "<<rounds<<" random field pairs passed\n";
}
struct AP{cpp_int x,y;bool inf=false;};
static cpp_int invp(cpp_int x,const cpp_int&p){return boost::multiprecision::powm(mod(x,p),p-2,p);}
static AP addref(AP A,AP B,const cpp_int&p,const cpp_int&a){
 if(A.inf)return B;if(B.inf)return A;
 if(A.x==B.x && mod(A.y+B.y,p)==0)return {{},{},true};
 cpp_int m;
 if(A.x==B.x&&A.y==B.y)m=mod((3*A.x*A.x+a)*invp(2*A.y,p),p);
 else m=mod((B.y-A.y)*invp(B.x-A.x,p),p);
 cpp_int x=mod(m*m-A.x-B.x,p),y=mod(m*(A.x-x)-A.y,p);return {x,y,false};
}
template<class F> AP jac_to_ref(const njac<F>&P,const cpp_int&p){
 cpp_int X=canon<F>(P.X),Y=canon<F>(P.Y),Z=canon<F>(P.Z);if(Z==0)return {{},{},true};
 cpp_int zi=invp(Z,p),z2=zi*zi%p,z3=z2*zi%p;return {X*z2%p,Y*z3%p,false};
}
template<class F> naff<F> aff_mont(const AP&A){return {F::from_canonical(raw<F>(A.x)),F::from_canonical(raw<F>(A.y))};}
template<class F> void point_test(const char*name,const char*ps,const char*xs,const char*ys){
 cpp_int p(ps),a=p-3;AP G{cpp_int(xs),cpp_int(ys),false};auto gm=aff_mont<F>(G);
 njac<F> J{gm.x,gm.y,F::one()};auto D=n_double<F>(J);AP want2=addref(G,G,p),got2=jac_to_ref<F>(D,p);
 CHECK(got2.x==want2.x&&got2.y==want2.y,name);
 auto T=n_madd<F>(D,gm);AP want3=addref(want2,G,p),got3=jac_to_ref<F>(T,p);
 CHECK(got3.x==want3.x&&got3.y==want3.y,name);
 std::cout<<name<<": Jacobian double + mixed add passed\n";
}
int main(){
 field_test<Fp256,P256Mod>("P-256",20000);
 field_test<Fp384,P384Mod>("P-384",20000);
 point_test<Fp256>("P-256 point",
  "0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
  "0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296",
  "0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
 point_test<Fp384>("P-384 point",
  "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff",
  "0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7",
  "0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f");
 if(fails){std::cerr<<fails<<" failures\n";return 1;}std::cout<<"ALL PASSED\n";return 0;
}
