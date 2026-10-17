// Memory-bounded version: retain seven baby-polynomial derivatives, then
// evaluate bounded batches of giant points. Same moment identities as v2.
#define CM_CONSTRUCTOR_LIBRARY
#include "sqrt_velu.cpp"
#include <algorithm>

std::array<Poly,3> paired_powers(const Z& x, const Z& A, const Z& F) {
  Z x2=red(x*x), A2=red(A*A), F2=red(F*F), AF=red(A*F);
  Poly u{},v{},w{};
  u[0]=x;u[1]=red(-A);u[2]=red(2*F);
  v[0]=x2;v[1]=red(-2*x*A-4*F);v[2]=red(A2+16*x*F);
  v[3]=red(-8*AF);v[4]=red(8*F2);
  w[0]=red(x2*x);w[1]=red(-3*x2*A-12*x*F);
  w[2]=red(3*x*A2+42*x2*F+12*AF);
  w[3]=red(-A2*A-60*x*AF-24*F2);
  w[4]=red(18*A2*F+96*x*F2);w[5]=red(-48*A*F2);w[6]=red(32*F2*F);
  return {u,v,w};
}

void polynomial_self_test() {
  uint64_t state=20261015;
  auto next=[&](){state^=state<<13;state^=state>>7;state^=state<<17;return state;};
  for(const char* modulus:{"5947698030651343","78125000378223625457769953627"}) {
    p=Z(modulus);
    auto sample=[&](){return red((Z(next())<<64)+next());};
    for(unsigned i=0;i<2048;i++) {
      Z x=sample(),A=sample(),F=sample();
      Poly u{},delta{};u[0]=x;u[1]=red(-A);u[2]=red(2*F);
      delta[1]=red(-4*F);delta[2]=red(12*x*F);
      delta[3]=red(-4*F*A);delta[4]=red(4*F*F);
      Poly v=product(u,u),w=product(v,u),ud=product(u,delta);
      add_scaled(v,delta,1);add_scaled(w,ud,3);
      if(paired_powers(x,A,F)!=std::array<Poly,3>{u,v,w})
        throw std::runtime_error("closed-form polynomial mismatch");
    }
  }
  std::cout<<"{\"self_test\":\"paired-polynomial coefficients\",\"cases\":4096,\"fields\":2,\"passed\":true}"<<std::endl;
}

int main(int argc, char** argv) {
  try {
    if(argc==2 && std::string(argv[1])=="--self-test") {polynomial_self_test();return 0;}
    const uint64_t batch = argc == 2 ? std::stoull(argv[1]) : 262144;
    if (!batch || batch > 2097152) throw std::runtime_error("invalid batch cap");
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
            Z value=(k+1)*h[k+1];
            for(unsigned j=1;j<=k;j++)value-=h[j]*L[k-j];
            L[k]=red(value);
            S[k+1]=red(k%2?-L[k]:L[k]);
          }
          auto powers=paired_powers(x,A,F);
          for(unsigned k=0;k<3;k++) {
            Z value=moments[k];
            for(unsigned j=0;j<=2*(k+1);j++)value+=2*powers[k][j]*S[j];
            moments[k]=red(value);
          }
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
