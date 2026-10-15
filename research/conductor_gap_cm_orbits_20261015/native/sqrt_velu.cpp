// Codomain-only square-root Velu by kernel power sums and fast multipoint
// evaluation. FLINT supplies product/remainder trees; GMP supplies F_p.
// Does not claim point-map or dual-isogeny certification.
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_vec.h>
#include <gmpxx.h>
#include <array>
#include <vector>
#include <iostream>
#include <chrono>
#include <cmath>
#include <stdexcept>
using Z=mpz_class; using Poly=std::array<Z,7>;
Z p,a,b;
Z red(Z z){mpz_mod(z.get_mpz_t(),z.get_mpz_t(),p.get_mpz_t());return z;}
Z inv(const Z&z){Z v;if(!mpz_invert(v.get_mpz_t(),z.get_mpz_t(),p.get_mpz_t()))throw std::runtime_error("zero denominator");return v;}
struct P{Z x,y;bool inf=true;};
P plus(P u,P v){if(u.inf)return v;if(v.inf)return u;if(u.x==v.x&&red(u.y+v.y)==0)return{};Z m=u.x==v.x?red((3*u.x*u.x+a)*inv(red(2*u.y))):red((v.y-u.y)*inv(red(v.x-u.x)));Z x=red(m*m-u.x-v.x);return{x,red(m*(u.x-x)-u.y),false};}
P times(P q,uint64_t k){P v;while(k){if(k&1)v=plus(v,q);q=plus(q,q);k>>=1;}return v;}
Z rhs(Z x){return red(x*x*x+a*x+b);}
Poly product(const Poly&u,const Poly&v){Poly w{};for(unsigned i=0;i<7;i++)for(unsigned j=0;i+j<7;j++)w[i+j]=red(w[i+j]+u[i]*v[j]);return w;}
void add_scaled(Poly&u,const Poly&v,unsigned c){for(unsigned i=0;i<7;i++)u[i]=red(u[i]+c*v[i]);}
void set(fmpz_t f,const Z&z){fmpz_set_mpz(f,z.get_mpz_t());}
Z get(const fmpz_t f){Z z;fmpz_get_mpz(z.get_mpz_t(),f);return z;}
int main(){try{uint64_t id,ell;Z kx,ky,d,a0,b0,expected_a,expected_b;while(std::cin>>id>>ell>>p>>a>>b>>kx>>ky>>d>>a0>>b0>>expected_a>>expected_b){auto start=std::chrono::steady_clock::now();P K{kx,ky,false};if(red(ky*ky)!=rhs(kx)||times(K,ell).inf==false)throw std::runtime_error("invalid prime-kernel witness");
uint64_t n=(ell-1)/2,B=std::max<uint64_t>(1,std::sqrt(n/2.0)),G=n/(2*B);if(G==0){B=1;G=0;}
fmpz_t mod;fmpz_init(mod);set(mod,p);fmpz_mod_ctx_t ctx;fmpz_mod_ctx_init(ctx,mod);fmpz_mod_poly_t H,D;fmpz_mod_poly_init(H,ctx);fmpz_mod_poly_init(D,ctx);
fmpz* roots=_fmpz_vec_init(B);fmpz* giant=_fmpz_vec_init(G);P step=times(K,2),v=K;for(uint64_t j=0;j<B;j++){if(v.inf)throw std::runtime_error("baby infinity");set(roots+j,v.x);v=plus(v,step);}v=times(K,2*B);P jump=times(K,4*B);for(uint64_t i=0;i<G;i++){if(v.inf)throw std::runtime_error("giant infinity");set(giant+i,v.x);v=plus(v,jump);}
fmpz_mod_poly_product_roots_fmpz_vec(H,roots,B,ctx);_fmpz_vec_clear(roots,B);
auto tree=_fmpz_mod_poly_tree_alloc(G);if(G)_fmpz_mod_poly_tree_build(tree,giant,G,ctx);std::array<fmpz*,7> eval{};for(unsigned k=0;k<7;k++){eval[k]=_fmpz_vec_init(G);if(G)_fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(eval[k],H->coeffs,H->length,tree,G,ctx);fmpz_mod_poly_derivative(D,H,ctx);fmpz_mod_poly_swap(D,H,ctx);}
_fmpz_mod_poly_tree_free(tree,G);std::array<Z,3> moments{};Z fact=1;std::array<Z,7> invfact;for(unsigned j=0;j<7;j++){if(j)fact*=j;invfact[j]=inv(fact);}
for(uint64_t i=0;i<G;i++){Z x=get(giant+i),F=rhs(x),A=red(3*x*x+a);Z reciprocal=inv(get(eval[0]+i));std::array<Z,7> h{},S{};S[0]=B;h[0]=1;for(unsigned k=1;k<=6;k++)h[k]=red(get(eval[k]+i)*invfact[k]*reciprocal);std::array<Z,6>L{};for(unsigned k=0;k<6;k++){L[k]=red((k+1)*h[k+1]);for(unsigned j=1;j<=k;j++)L[k]=red(L[k]-h[j]*L[k-j]);S[k+1]=red(k%2?-L[k]:L[k]);}
Poly u{},delta{};u[0]=x;u[1]=red(-A);u[2]=red(2*F);delta[1]=red(-4*F);delta[2]=red(12*x*F);delta[3]=red(-4*F*A);delta[4]=red(4*F*F);Poly u2=product(u,u),u3=product(u2,u),ud=product(u,delta);add_scaled(u2,delta,1);add_scaled(u3,ud,3);std::array<Poly,3>powers{u,u2,u3};for(unsigned k=0;k<3;k++)for(unsigned j=0;j<7;j++)moments[k]=red(moments[k]+2*powers[k][j]*S[j]);}
v=times(K,4*B*G+1);for(uint64_t index=4*B*G+1;index<ell;index+=2){if(v.inf)throw std::runtime_error("remainder infinity");Z x=v.x,x2=red(x*x);moments[0]=red(moments[0]+x);moments[1]=red(moments[1]+x2);moments[2]=red(moments[2]+x2*x);v=plus(v,step);}
Z dp=d;for(unsigned k=0;k<3;k++){moments[k]=red(moments[k]*dp);dp=red(dp*d);}Z T=red(3*moments[1]+a0*n),W=red(5*moments[2]+3*a0*moments[0]+2*b0*n),a1=red(a0-10*T),b1=red(b0-14*W);bool same=(expected_a<0||(a1==expected_a&&b1==expected_b));double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cout<<"{\"id\":"<<id<<",\"ell\":"<<ell<<",\"p\":\""<<p<<"\",\"a\":\""<<a1<<"\",\"b\":\""<<b1<<"\",\"baby_points\":"<<B<<",\"giant_points\":"<<G<<",\"covered_kernel_pairs\":"<<n<<",\"reference_supplied\":"<<(expected_a>=0?"true":"false")<<",\"matches_direct_codomain\":"<<(expected_a<0?"null":same?"true":"false")<<",\"wall_seconds\":"<<seconds<<"}"<<std::endl;
for(auto e:eval)_fmpz_vec_clear(e,G);_fmpz_vec_clear(giant,G);fmpz_mod_poly_clear(H,ctx);fmpz_mod_poly_clear(D,ctx);fmpz_mod_ctx_clear(ctx);fmpz_clear(mod);if(!same)return 3;}
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}}
