// Memory-bounded version: retain seven baby-polynomial derivatives, then
// evaluate bounded batches of giant points. Same moment identities as v2.
#define CM_CONSTRUCTOR_LIBRARY
#include "sqrt_velu.cpp"
#include <algorithm>

int main(int argc, char** argv) {
  try {
    const uint64_t batch = argc == 2 ? std::stoull(argv[1]) : 262144;
    if (!batch || batch > 1048576) throw std::runtime_error("invalid batch cap");
    uint64_t id, ell;
    Z kx, ky, d, a0, b0, expected_a, expected_b;
    while (std::cin >> id >> ell >> p >> a >> b >> kx >> ky >> d >> a0 >> b0 >> expected_a >> expected_b) {
      auto start = std::chrono::steady_clock::now();
      P K{kx, ky, false};
      if (red(ky*ky) != rhs(kx) || !times(K,ell).inf)
        throw std::runtime_error("invalid prime-kernel witness");
      uint64_t n=(ell-1)/2, B=std::max<uint64_t>(1,std::sqrt(n/2.0)), G=n/(2*B);
      if (G==0) {B=1;G=0;}
      fmpz_t mod;fmpz_init(mod);set(mod,p);
      fmpz_mod_ctx_t ctx;fmpz_mod_ctx_init(ctx,mod);
      fmpz_mod_poly_t derivatives[7];
      for (auto& h:derivatives) fmpz_mod_poly_init(h,ctx);
      fmpz* roots=_fmpz_vec_init(B);
      P step=times(K,2), v=K;
      for (uint64_t j=0;j<B;j++) {
        if(v.inf) throw std::runtime_error("baby infinity");
        set(roots+j,v.x);v=plus(v,step);
      }
      fmpz_mod_poly_product_roots_fmpz_vec(derivatives[0],roots,B,ctx);
      _fmpz_vec_clear(roots,B);
      for(unsigned k=1;k<7;k++) fmpz_mod_poly_derivative(derivatives[k],derivatives[k-1],ctx);
      std::array<Z,3> moments{};
      Z fact=1;std::array<Z,7> invfact;
      for(unsigned j=0;j<7;j++){if(j)fact*=j;invfact[j]=inv(fact);}
      v=times(K,2*B);P jump=times(K,4*B);
      uint64_t batches=0;
      for (uint64_t offset=0;offset<G;offset+=batch) {
        uint64_t count=std::min(batch,G-offset);batches++;
        fmpz* giant=_fmpz_vec_init(count);
        for (uint64_t i=0;i<count;i++) {
          if(v.inf) throw std::runtime_error("giant infinity");
          set(giant+i,v.x);v=plus(v,jump);
        }
        auto tree=_fmpz_mod_poly_tree_alloc(count);
        _fmpz_mod_poly_tree_build(tree,giant,count,ctx);
        std::array<fmpz*,7> eval{};
        for(unsigned k=0;k<7;k++) {
          eval[k]=_fmpz_vec_init(count);
          _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(eval[k],derivatives[k]->coeffs,derivatives[k]->length,tree,count,ctx);
        }
        _fmpz_mod_poly_tree_free(tree,count);
        for(uint64_t i=0;i<count;i++) {
          Z x=get(giant+i),F=rhs(x),A=red(3*x*x+a);
          Z reciprocal=inv(get(eval[0]+i));
          std::array<Z,7> h{},S{};S[0]=B;h[0]=1;
          for(unsigned k=1;k<=6;k++)h[k]=red(get(eval[k]+i)*invfact[k]*reciprocal);
          std::array<Z,6>L{};
          for(unsigned k=0;k<6;k++) {
            L[k]=red((k+1)*h[k+1]);
            for(unsigned j=1;j<=k;j++)L[k]=red(L[k]-h[j]*L[k-j]);
            S[k+1]=red(k%2?-L[k]:L[k]);
          }
          Poly u{},delta{};u[0]=x;u[1]=red(-A);u[2]=red(2*F);
          delta[1]=red(-4*F);delta[2]=red(12*x*F);
          delta[3]=red(-4*F*A);delta[4]=red(4*F*F);
          Poly u2=product(u,u),u3=product(u2,u),ud=product(u,delta);
          add_scaled(u2,delta,1);add_scaled(u3,ud,3);
          std::array<Poly,3> powers{u,u2,u3};
          for(unsigned k=0;k<3;k++)for(unsigned j=0;j<7;j++)moments[k]=red(moments[k]+2*powers[k][j]*S[j]);
        }
        for(auto e:eval)_fmpz_vec_clear(e,count);
        _fmpz_vec_clear(giant,count);
        std::cerr<<"case "<<id<<" giant points "<<offset+count<<" / "<<G<<"\n";
      }
      v=times(K,4*B*G+1);
      for(uint64_t index=4*B*G+1;index<ell;index+=2) {
        if(v.inf)throw std::runtime_error("remainder infinity");
        Z x=v.x,x2=red(x*x);moments[0]=red(moments[0]+x);
        moments[1]=red(moments[1]+x2);moments[2]=red(moments[2]+x2*x);v=plus(v,step);
      }
      Z dp=d;for(unsigned k=0;k<3;k++){moments[k]=red(moments[k]*dp);dp=red(dp*d);}
      Z T=red(3*moments[1]+a0*n),W=red(5*moments[2]+3*a0*moments[0]+2*b0*n);
      Z a1=red(a0-10*T),b1=red(b0-14*W);
      bool same=expected_a>=0&&a1==expected_a&&b1==expected_b;
      double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      std::cout<<"{\"id\":"<<id<<",\"ell\":"<<ell<<",\"p\":\""<<p<<"\",\"a\":\""<<a1<<"\",\"b\":\""<<b1<<"\",\"baby_points\":"<<B<<",\"giant_points\":"<<G<<",\"giant_batch_cap\":"<<batch<<",\"batches\":"<<batches<<",\"covered_kernel_pairs\":"<<n<<",\"reference_supplied\":"<<(expected_a>=0?"true":"false")<<",\"matches_direct_codomain\":"<<(expected_a<0?"null":same?"true":"false")<<",\"wall_seconds\":"<<seconds<<"}"<<std::endl;
      for(auto& h:derivatives) fmpz_mod_poly_clear(h,ctx);
      fmpz_mod_ctx_clear(ctx);fmpz_clear(mod);
      if(expected_a>=0&&!same) return 3;
    }
  } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 2;}
}
