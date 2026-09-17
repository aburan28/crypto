#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include "ecc2k130/include/curveparams.h"
#include "ecc2k130/include/walk.h"
using C=CfgF23; using R=Ref<C>; using E=R::Elem; using P=R::Point;
struct K{unsigned long long v[3];bool operator==(K const&o)const{return v[0]==o.v[0]&&v[1]==o.v[1]&&v[2]==o.v[2];}};
struct H{size_t operator()(K const&k)const{return size_t((k.v[0]*0x9e3779b97f4a7c15ULL)^(k.v[1]*0xbf58476d1ce4e5b9ULL)^(k.v[2]*0x94d049bb133111ebULL));}};
static K key(E x){x=R::canonical(x);return K{{x.v[0],x.v[1],x.v[2]}};}
struct W{P p;U192 a,b;unsigned long long it=0,seed=0;};
struct Stats{std::vector<unsigned long long> totals,walks;unsigned solved=0,bad=0,invariantBad=0,exceptional=0,overdue=0;unsigned long long zeroDen=0;};
static int pos[C::M];
static void initPositions(){int g=1;for(int p=0;p<C::M;p++){pos[g-1]=p;g=2*g%C::NRING;if(g>C::M)g=C::NRING-g;}}
struct Choice{int phase,index,sign;};
static Choice choose(P const&p){
 int w=R::weight(p.x),sum=0;if(w%C::M==0)return{-1,-1,0};
 for(int i=0;i<C::M;i++)if(R::bit(p.x,i))sum=(sum+pos[i])%C::M;
 int inv=1;while(w*inv%C::M!=1)++inv;int ph=sum*inv%C::M,best=C::M,index=-1;
 for(int i=0;i<C::M;i++)if(R::bit(p.x,i)){int q=(pos[i]-ph+C::M)%C::M;if(q<best){best=q;index=i;}}
 return{ph,index,R::bit(p.y,index)};
}
static void proveSelector(P B){
 unsigned long long checked=0;
 for(unsigned long long bits=1;bits<(1ull<<C::M)-1;bits++){
  P p{};p.x=R::zero();p.y=R::zero();p.x.v[0]=bits;Choice a=choose(p);P q=p;q.x=R::sigma(p.x,1);Choice b=choose(q);
  if(a.phase<0||b.phase!=(a.phase+1)%C::M||pos[b.index]!=(pos[a.index]+1)%C::M){std::fprintf(stderr,"selector proof failed bits=%llu\n",bits);std::exit(2);}checked++;
 }
 P p=B;unsigned valid=0,exceptional=0;
 for(unsigned i=0;i<131072;i++){
  Choice a=choose(p);if(a.phase<0)exceptional++;else{Choice f=choose(R::frob(p,1)),n=choose(R::neg(p));if(f.sign!=a.sign||n.sign==a.sign||f.phase!=(a.phase+1)%C::M){std::fprintf(stderr,"sign proof failed point=%u\n",i);std::exit(3);}valid++;}p=R::addPt(p,B);
 }
 std::printf("selector_identity_checked=%llu valid_point_sign_checks=%u exceptional_points=%u passed=1\n",checked,valid,exceptional);
}
static void show(const char*l,std::vector<unsigned long long>v){std::sort(v.begin(),v.end());long double s=0;for(auto x:v)s+=x;std::printf("%s_mean=%.6Lf %s_median=%llu %s_p90=%llu",l,s/v.size(),l,v[v.size()/2],l,v[v.size()*9/10]);}
int main(int argc,char**argv){std::setvbuf(stdout,nullptr,_IOLBF,0);unsigned trials=argc>1?unsigned(std::strtoul(argv[1],nullptr,10)):100;int first=argc>2?std::atoi(argv[2]):0,last=argc>2?first+1:2;initPositions();
 U192 ell=u192_from_dec(eccF23::ELL_DEC),s=u192_from_dec(eccF23::S_DEC),sp[256];sp[0]=u192_from(1);for(int i=1;i<256;i++)sp[i]=mod_mul(sp[i-1],s,ell);
 P B=R::make(R::fromLimbs(eccF23::PX),R::fromLimbs(eccF23::PY)),Q=R::make(R::fromLimbs(eccF23::QX),R::fromLimbs(eccF23::QY));proveSelector(B);const unsigned rv[8]={1,3,5,7,11,13,17,19},tv[8]={2,5,11,17,23,29,31,37};P orbit[8][C::M];
 for(int jx=0;jx<8;jx++){P A=R::addPt(R::scalarMul(B,u192_from(rv[jx])),R::scalarMul(Q,u192_from(tv[jx])));if(A.inf){std::fprintf(stderr,"decorrelated addend %d is infinity\n",jx);return 4;}for(int e=0;e<C::M;e++)orbit[jx][e]=R::frob(A,e);}
 for(int mode=first;mode<last;mode++){Stats z;for(unsigned trial=0;trial<trials;trial++){if(mode)std::fprintf(stderr,"signed_trial=%u\n",trial);std::unordered_map<K,W,H>tab;unsigned long long total=0;
  for(unsigned wi=0;wi<20000;wi++){W w;w.seed=eccSeedFor(30000+trial,wi);w.p=R::startPoint(w.seed,B,Q,&w.a,ell,sp);w.b=u192_from(1);
   for(;;w.it++){if(w.it>=1024){z.overdue++;w.p.inf=true;break;}int hw=R::weight(w.p.x);if(hw<=eccF23::DP_WEIGHT)break;int j=R::jOf(hw);
    if(!mode){U192 f=mod_add(u192_from(1),sp[j],ell);w.a=mod_mul(w.a,f,ell);w.b=mod_mul(w.b,f,ell);w.p=R::addPt(w.p,R::frob(w.p,j));}
    else{Choice c=choose(w.p);if(c.phase<0){U192 f=mod_add(u192_from(1),sp[j],ell);w.a=mod_mul(w.a,f,ell);w.b=mod_mul(w.b,f,ell);w.p=R::addPt(w.p,R::frob(w.p,j));z.exceptional++;}
     else{int ix=j-3,e=c.phase;U192 da=mod_mul(u192_from(rv[ix]),sp[e],ell),db=mod_mul(u192_from(tv[ix]),sp[e],ell);if(c.sign){w.a=mod_sub(w.a,da,ell);w.b=mod_sub(w.b,db,ell);w.p=R::addPt(w.p,R::neg(orbit[ix][e]));}else{w.a=mod_add(w.a,da,ell);w.b=mod_add(w.b,db,ell);w.p=R::addPt(w.p,orbit[ix][e]);}}}
    if(w.p.inf){z.bad++;break;}}
   total+=w.it;if(w.p.inf)continue;K k=key(w.p.x);auto [at,fresh]=tab.emplace(k,w);if(fresh)continue;W const&a=w;W const&b=at->second;
   P sa=R::addPt(R::scalarMul(B,a.a),R::scalarMul(Q,a.b)),sb=R::addPt(R::scalarMul(B,b.a),R::scalarMul(Q,b.b));if(!R::eq(sa,a.p)||!R::eq(sb,b.p)){z.invariantBad++;continue;}
   bool ok=false;for(int rot=0;rot<C::M&&!ok;rot++){P q=R::frob(b.p,rot);int eps=R::eq(q,a.p)?1:(R::eq(R::neg(q),a.p)?-1:0);if(!eps)continue;U192 sc=sp[rot];if(eps<0)sc=mod_neg(sc,ell);U192 num=mod_sub(a.a,mod_mul(sc,b.a,ell),ell),den=mod_sub(mod_mul(sc,b.b,ell),a.b,ell);if(u192_is_zero(den)){z.zeroDen++;continue;}U192 x=mod_mul(num,mod_inv(den,ell),ell);ok=R::eq(R::scalarMul(B,x),Q);}
   if(ok){z.totals.push_back(total);z.walks.push_back(wi+1);z.solved++;break;}}
  }
  const char*name=mode?"decorrelated-signed-covariant-add":"selected";std::printf("mode=%s trials=%u solved=%u bad=%u invariant_bad=%u exceptional_steps=%u overdue=%u zero_den=%llu ",name,trials,z.solved,z.bad,z.invariantBad,z.exceptional,z.overdue,z.zeroDen);show("iterations",z.totals);std::printf(" ");show("walks",z.walks);std::printf("\n");
 }
}
