#include <algorithm>
#include <cstdio>
#include <unordered_map>
#include <vector>
#include "ecc2k130/include/curveparams.h"
#include "ecc2k130/include/walk.h"
using C=CfgF23; using R=Ref<C>; using E=R::Elem; using P=R::Point;
struct K { unsigned long long v[3]; bool operator==(K const&o) const { return v[0]==o.v[0]&&v[1]==o.v[1]&&v[2]==o.v[2]; } };
struct H { size_t operator()(K const&k) const { return size_t((k.v[0]*0x9e3779b97f4a7c15ULL)^(k.v[1]*0xbf58476d1ce4e5b9ULL)^(k.v[2]*0x94d049bb133111ebULL)); } };
static K key(E x) { x=R::canonical(x); return K{{x.v[0],x.v[1],x.v[2]}}; }
struct W { P p; U192 alpha; unsigned long long n[8]{}; unsigned long long it=0,seed=0; };
struct Stats { std::vector<unsigned long long> totals,walks; unsigned solved=0,bad=0; };
static void show(const char*label,std::vector<unsigned long long> v) {
    std::sort(v.begin(),v.end()); long double sum=0; for(auto x:v) sum+=x;
    std::printf("%s_mean=%.6Lf %s_median=%llu %s_p90=%llu",label,sum/v.size(),label,v[v.size()/2],label,v[v.size()*9/10]);
}
int main() {
 U192 ell=u192_from_dec(eccF23::ELL_DEC),s=u192_from_dec(eccF23::S_DEC),sp[256]; sp[0]=u192_from(1); for(int i=1;i<256;i++)sp[i]=mod_mul(sp[i-1],s,ell);
 P B=R::make(R::fromLimbs(eccF23::PX),R::fromLimbs(eccF23::PY)),Q=R::make(R::fromLimbs(eccF23::QX),R::fromLimbs(eccF23::QY)),addend[8];
 for(int j=0;j<8;j++) addend[j]=R::frob(B,j+3);
 for(int mode=0;mode<2;mode++) { Stats z;
  for(unsigned trial=0;trial<2000;trial++) { std::unordered_map<K,W,H> tab; unsigned long long total=0;
   for(unsigned wi=0;wi<20000;wi++) { W w; w.seed=eccSeedFor(20000+trial,wi); w.p=R::startPoint(w.seed,B,Q,&w.alpha,ell,sp);
    for(;;w.it++) { int hw=R::weight(w.p.x); if(hw<=eccF23::DP_WEIGHT) break; int ix=(hw>>1)&7; w.n[ix]++;
      w.p=mode?R::addPt(w.p,addend[ix]):R::addPt(w.p,R::frob(w.p,ix+3)); if(w.p.inf){z.bad++;break;} }
    if(w.p.inf) continue; total+=w.it; K k=key(w.p.x); auto [pos,fresh]=tab.emplace(k,w); if(fresh) continue;
    W const&a=w; W const&b=pos->second; U192 aa=a.alpha,ab=b.alpha,ba=u192_from(1),bb=u192_from(1);
    if(mode) { for(int j=0;j<8;j++) { U192 ea=u192_zero(),eb=u192_zero(); ea.v[0]=a.n[j]; eb.v[0]=b.n[j]; aa=mod_add(aa,mod_mul(ea,sp[j+3],ell),ell); ab=mod_add(ab,mod_mul(eb,sp[j+3],ell),ell); } }
    else { auto mult=[&](W const&q){ U192 mu=u192_from(1); for(int j=0;j<8;j++) if(q.n[j]) { U192 f=mod_add(u192_from(1),sp[j+3],ell),e=u192_zero(); e.v[0]=q.n[j]; mu=mod_mul(mu,mod_pow(f,e,ell),ell); } return mu; }; U192 ma=mult(a),mb=mult(b); aa=mod_mul(ma,aa,ell); ab=mod_mul(mb,ab,ell); ba=ma; bb=mb; }
    bool ok=false; for(int rot=0;rot<C::M&&!ok;rot++) { P q=R::frob(b.p,rot); int eps=R::eq(q,a.p)?1:(R::eq(R::neg(q),a.p)?-1:0); if(!eps) continue; U192 sc=sp[rot]; if(eps<0)sc=mod_neg(sc,ell);
      U192 num=mod_sub(aa,mod_mul(sc,ab,ell),ell), den=mod_sub(mod_mul(sc,bb,ell),ba,ell);
      if(u192_is_zero(den)) continue; U192 x=mod_mul(num,mod_inv(den,ell),ell); ok=R::eq(R::scalarMul(B,x),Q); }
    if(ok) { z.totals.push_back(total); z.walks.push_back(wi+1); z.solved++; break; }
   }
  }
  std::printf("mode=%s solved=%u bad=%u ",mode?"fixed-add":"selected",z.solved,z.bad); show("iterations",z.totals); std::printf(" "); show("walks",z.walks); std::printf("\n");
 }
}
