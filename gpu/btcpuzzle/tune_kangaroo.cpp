/* tune_kangaroo.cpp -- measure the kangaroo constant C in steps = C*sqrt(W).
 *
 * Kangaroo run lengths are close to exponentially distributed, so a mean
 * over a handful of trials tells you almost nothing: with 20 trials the
 * standard error is over 20% of the mean, and a grid search across twenty
 * parameter combinations will hand you a "winner" that is pure sampling
 * noise.  This reports the standard error alongside the mean so the
 * difference between a real effect and a lucky draw is visible.
 *
 * Built separately from `make test` because it is slow -- 120 trials per
 * configuration takes minutes.
 *
 *   g++ -O2 -std=c++17 -I. -I../ecc -o tune_kangaroo tune_kangaroo.cpp
 *   ./tune_kangaroo
 */
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include "kangaroo_host.hpp"
struct HS {
    std::vector<uint32_t> X,Y,D,st,rs; std::vector<kg_dp> dps; uint32_t n=0; kg_ctx c;
    void init(KangarooHost &h, uint32_t T, uint32_t W, uint32_t cap){
        uint32_t k=T*W; X.assign(8*k,0);Y.assign(8*k,0);D.assign(8*k,0);st.assign(k,0);rs.assign(k,0);
        dps.resize(cap); n=0;
        c.X=X.data();c.Y=Y.data();c.D=D.data();c.steps=st.data();c.restarts=rs.data();
        c.nthreads=T;c.kang_per_thread=W;c.jumps=h.jumps.data();c.Qshift=h.Qshift;c.prm=h.prm;
        c.dp_out=dps.data();c.dp_count=&n;c.dp_cap=cap;
        for(uint32_t t=0;t<T;t++) kg_init_thread(c,t);
    }
};
static double solve(int bits,uint32_t seed,int ms,uint32_t njb,uint32_t dpb,int*ok){
    u256 secret=u256_pow2(bits-1);
    uint64_t s=0x9E3779B97F4A7C15ull*seed+12345; u256 off=u256_zero();
    for(int i=0;i<4;i++){uint64_t z=kg_splitmix64(s);off.v[2*i]=(uint32_t)z;off.v[2*i+1]=(uint32_t)(z>>32);}
    for(int l=0;l<8;l++){int lo=32*l;if(lo>=bits-1)off.v[l]=0;else if(lo+32>bits-1)off.v[l]&=(1u<<(bits-1-lo))-1u;}
    secret=u256_add(secret,off);
    KangarooHost h; h.prm.njump_bits=njb; h.prm.dp_mask=(1u<<dpb)-1u; h.prm.max_steps=1u<<22;
    h.prm.seed=seed; h.prm.reseed_on_dp=1; h.mean_shift=ms;
    h.setup(Curve::to_affine(Curve::scalar_mul(Curve::generator(),secret.v,0)),
            u256_pow2(bits-1),(uint32_t)(bits-1));
    HS st; st.init(h,4,8,1u<<18);
    double total=0; uint32_t consumed=0; u256 found; *ok=0;
    for(int it=0;it<8000000;it++){
        for(uint32_t t=0;t<4;t++) kg_step_batch<8>(st.c,t);
        total+=32;
        while(consumed<st.n&&consumed<st.c.dp_cap)
            if(h.add_dp(st.dps[consumed++],found)){*ok=(u256_cmp(found,secret)==0);return total;}
        if(st.n>=st.c.dp_cap) break;
    }
    return total;
}
int main(){
    const int bits=32, trials=120; const uint32_t dpb=7;
    double rootw=ldexp(1.0,(bits-1)/2.0);
    struct { const char*name; int ms; uint32_t njb; } cases[] = {
        {"apparent best (2^5, shift-2)", -2, 5},
        {"apparent worst (2^7, shift+0)", 0, 7},
        {"current default (2^6, shift+0)", 0, 6},
    };
    printf("interval 2^%d, %d trials per configuration\n", bits-1, trials);
    printf("%-32s %8s %8s %8s\n","configuration","mean","std err","median");
    for (auto &c : cases) {
        std::vector<double> v;
        for(int i=0;i<trials;i++){int ok=0;double x=solve(bits,9000+i,c.ms,c.njb,dpb,&ok);if(ok)v.push_back(x/rootw);}
        double m=0; for(double x:v) m+=x; m/=v.size();
        double s2=0; for(double x:v) s2+=(x-m)*(x-m); double se=sqrt(s2/(v.size()-1)/v.size());
        std::sort(v.begin(),v.end());
        printf("%-32s %8.2f %8.2f %8.2f   (n=%zu)\n",c.name,m,se,v[v.size()/2],v.size());
        fflush(stdout);
    }
    return 0;
}
