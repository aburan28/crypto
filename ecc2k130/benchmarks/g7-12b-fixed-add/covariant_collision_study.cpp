#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include "ecc2k130/include/curveparams.h"
#include "ecc2k130/include/walk.h"

using C=CfgF23; using R=Ref<C>; using E=R::Elem; using P=R::Point;
struct K { unsigned long long v[3]; bool operator==(K const&o) const { return v[0]==o.v[0]&&v[1]==o.v[1]&&v[2]==o.v[2]; } };
struct H { size_t operator()(K const&k) const { return size_t((k.v[0]*0x9e3779b97f4a7c15ULL)^(k.v[1]*0xbf58476d1ce4e5b9ULL)^(k.v[2]*0x94d049bb133111ebULL)); } };
static K key(E x) { x=R::canonical(x); return K{{x.v[0],x.v[1],x.v[2]}}; }
struct W { P p; U192 a,b; unsigned long long it=0,seed=0; };
struct Stats { std::vector<unsigned long long> totals,walks; unsigned solved=0,bad=0,invariantBad=0,exceptional=0,overdue=0; unsigned long long degenerateExact=0; };
static int positions[C::M];

static void initPositions() {
    int gamma=1;
    for(int p=0;p<C::M;p++) {
        positions[gamma-1]=p;
        gamma=(2*gamma)%C::NRING;
        if(gamma>C::M) gamma=C::NRING-gamma;
    }
}
static int phase(E const&x) {
    int w=R::weight(x), sum=0;
    if(w%C::M==0) return -1;
    for(int i=0;i<C::M;i++) if(R::bit(x,i)) sum=(sum+positions[i])%C::M;
    int inv=1; while((w*inv)%C::M!=1) ++inv;
    return (sum*inv)%C::M;
}
static unsigned long long provePhase() {
    unsigned long long checked=0;
    for(unsigned long long bits=1;bits<(1ull<<C::M)-1;bits++) {
        E x=R::zero(); x.v[0]=bits;
        int p=phase(x), q=phase(R::sigma(x,1));
        if(p<0 || q!=(p+1)%C::M) {
            std::fprintf(stderr,"phase proof failed bits=%llu p=%d q=%d\n",bits,p,q);
            std::exit(2);
        }
        checked++;
    }
    return checked;
}
static void show(const char*label,std::vector<unsigned long long> v) {
    std::sort(v.begin(),v.end()); long double sum=0; for(auto x:v) sum+=x;
    std::printf("%s_mean=%.6Lf %s_median=%llu %s_p90=%llu",label,sum/v.size(),label,v[v.size()/2],label,v[v.size()*9/10]);
}
int main(int argc,char**argv) {
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    unsigned trials=argc>1?unsigned(std::strtoul(argv[1],nullptr,10)):100;
    initPositions();
    auto checked=provePhase();
    std::printf("phase_identity_checked=%llu exceptional_coordinates=2 passed=1\n",checked);
    U192 ell=u192_from_dec(eccF23::ELL_DEC),s=u192_from_dec(eccF23::S_DEC),sp[256];
    sp[0]=u192_from(1); for(int i=1;i<256;i++)sp[i]=mod_mul(sp[i-1],s,ell);
    P B=R::make(R::fromLimbs(eccF23::PX),R::fromLimbs(eccF23::PY));
    P Q=R::make(R::fromLimbs(eccF23::QX),R::fromLimbs(eccF23::QY));
    P orbit[C::M],fixed[8]; for(int e=0;e<C::M;e++) orbit[e]=R::frob(B,e); for(int j=0;j<8;j++) fixed[j]=orbit[j+3];
    int firstMode=argc>2?std::atoi(argv[2]):0, lastMode=argc>2?firstMode+1:3;
    for(int mode=firstMode;mode<lastMode;mode++) { Stats z;
        unsigned trialOffset=argc>3?unsigned(std::strtoul(argv[3],nullptr,10)):0;
        for(unsigned localTrial=0;localTrial<trials;localTrial++) { unsigned trial=trialOffset+localTrial; if(mode==2) std::fprintf(stderr,"covariant_trial=%u\n",trial); std::unordered_map<K,W,H> tab; unsigned long long total=0;
            for(unsigned wi=0;wi<20000;wi++) { W w; w.seed=eccSeedFor(20000+trial,wi);
                w.p=R::startPoint(w.seed,B,Q,&w.a,ell,sp); w.b=u192_from(1);
                for(;;w.it++) {
                    if(w.it>=1024) { z.overdue++; std::fprintf(stderr,"overdue trial=%u walk=%u seed=%llu\n",trial,wi,w.seed); w.p.inf=true; break; }
                    int hw=R::weight(w.p.x); if(hw<=eccF23::DP_WEIGHT) break;
                    int j=R::jOf(hw), ph=mode==2?phase(w.p.x):-1;
                    if(mode==0 || (mode==2 && ph<0)) {
                        U192 f=mod_add(u192_from(1),sp[j],ell);
                        w.a=mod_mul(w.a,f,ell); w.b=mod_mul(w.b,f,ell);
                        w.p=R::addPt(w.p,R::frob(w.p,j));
                        if(mode==2 && ph<0) z.exceptional++;
                    } else {
                        int exponent=mode==1?j:(ph+j)%C::M;
                        w.a=mod_add(w.a,sp[exponent],ell);
                        w.p=R::addPt(w.p,mode==1?fixed[j-3]:orbit[exponent]);
                    }
                    if(w.p.inf){z.bad++;break;}
                }
                total+=w.it;
                if(w.p.inf) continue;
                K k=key(w.p.x); auto [pos,fresh]=tab.emplace(k,w); if(fresh) continue;
                W const&a=w; W const&b=pos->second;
                if(R::eq(a.p,b.p) && u192_is_zero(mod_sub(b.b,a.b,ell))) { z.degenerateExact++; continue; }
                P scalarA=R::addPt(R::scalarMul(B,a.a),R::scalarMul(Q,a.b));
                P scalarB=R::addPt(R::scalarMul(B,b.a),R::scalarMul(Q,b.b));
                if(!R::eq(scalarA,a.p) || !R::eq(scalarB,b.p)) { z.invariantBad++; continue; }
                bool ok=false;
                for(int rot=0;rot<C::M&&!ok;rot++) { P q=R::frob(b.p,rot); int eps=R::eq(q,a.p)?1:(R::eq(R::neg(q),a.p)?-1:0); if(!eps) continue;
                    U192 sc=sp[rot]; if(eps<0)sc=mod_neg(sc,ell);
                    U192 num=mod_sub(a.a,mod_mul(sc,b.a,ell),ell), den=mod_sub(mod_mul(sc,b.b,ell),a.b,ell);
                    if(u192_is_zero(den)) continue; U192 x=mod_mul(num,mod_inv(den,ell),ell); ok=R::eq(R::scalarMul(B,x),Q);
                }
                if(ok) { z.totals.push_back(total); z.walks.push_back(wi+1); z.solved++; break; }
            }
        }
        const char*name=mode==0?"selected":(mode==1?"fixed-add":"covariant-fixed-add");
        std::printf("mode=%s trials=%u solved=%u bad=%u invariant_bad=%u exceptional_steps=%u overdue=%u degenerate_exact=%llu ",name,trials,z.solved,z.bad,z.invariantBad,z.exceptional,z.overdue,z.degenerateExact);
        show("iterations",z.totals); std::printf(" "); show("walks",z.walks); std::printf("\n");
    }
}
