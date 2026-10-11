#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include "ecc2k130/include/curveparams.h"
#include "ecc2k130/include/walk.h"
using C=CfgF23; using R=Ref<C>; using E=R::Elem; using P=R::Point;
struct K { unsigned long long x[3],y[3]; bool operator==(K const&o)const{for(int i=0;i<3;i++)if(x[i]!=o.x[i]||y[i]!=o.y[i])return false;return true;} };
struct H { size_t operator()(K const&k)const{size_t h=0;for(int i=0;i<3;i++){h^=size_t(k.x[i]+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2));h^=size_t(k.y[i]+0xbf58476d1ce4e5b9ULL+(h<<6)+(h>>2));}return h;} };
static int pos[C::M];
static void init(){int g=1;for(int p=0;p<C::M;p++){pos[g-1]=p;g=2*g%C::NRING;if(g>C::M)g=C::NRING-g;}}
static int phase(E const&x){int w=R::weight(x),sum=0;if(w%C::M==0)return-1;for(int i=0;i<C::M;i++)if(R::bit(x,i))sum=(sum+pos[i])%C::M;int inv=1;while(w*inv%C::M!=1)++inv;return sum*inv%C::M;}
static K key(P const&p){K k{};for(int i=0;i<3;i++){k.x[i]=p.x.v[i];k.y[i]=p.y.v[i];}return k;}
int main(){
 init(); unsigned trial=94,walk=73; unsigned long long seed=eccSeedFor(20000+trial,walk);
 U192 ell=u192_from_dec(eccF23::ELL_DEC),s=u192_from_dec(eccF23::S_DEC),sp[256];sp[0]=u192_from(1);for(int i=1;i<256;i++)sp[i]=mod_mul(sp[i-1],s,ell);
 P B=R::make(R::fromLimbs(eccF23::PX),R::fromLimbs(eccF23::PY)),Q=R::make(R::fromLimbs(eccF23::QX),R::fromLimbs(eccF23::QY)),orbit[C::M];for(int i=0;i<C::M;i++)orbit[i]=R::frob(B,i);
 U192 a;P p=R::startPoint(seed,B,Q,&a,ell,sp);std::unordered_map<K,unsigned long long,H> seen;int minWeight=C::M+1;unsigned long long exceptional=0;
 for(unsigned long long it=0;it<1000000;it++){
  int hw=R::weight(p.x);if(hw<minWeight)minWeight=hw;
  if(hw<=eccF23::DP_WEIGHT){std::printf("unexpected_dp seed=%llu iteration=%llu weight=%d\n",seed,it,hw);return 2;}
  auto [q,fresh]=seen.emplace(key(p),it);if(!fresh){std::printf("seed=%llu trial=%u walk=%u preperiod=%llu cycle_length=%llu min_weight=%d dp_weight=%d exceptional_steps=%llu exact_repeat=1\n",seed,trial,walk,q->second,it-q->second,minWeight,eccF23::DP_WEIGHT,exceptional);return 0;}
  int j=R::jOf(hw),ph=phase(p.x);if(ph<0){p=R::addPt(p,R::frob(p,j));exceptional++;}else p=R::addPt(p,orbit[(ph+j)%C::M]);
  if(p.inf){std::printf("unexpected_infinity seed=%llu iteration=%llu\n",seed,it+1);return 3;}
 }
 std::printf("no_repeat_within=1000000 seed=%llu min_weight=%d\n",seed,minWeight);return 4;
}
