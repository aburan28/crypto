#include "group.hpp"
#include "campaign_params.hpp"
#include <thread>
#include <mutex>
#include <atomic>
#include <unordered_set>
#include <map>
#include <iomanip>
#include <time.h>
#include <algorithm>
#include <random>
U power(U x,U k,U m){U v=1;for(;k;k>>=1,x=(__uint128_t)x*x%m)if(k&1)v=(__uint128_t)v*x%m;return v;}
U root(U n){if(!n)return 0;if(power(n,(p-1)/2,p)!=1)return p;if(p%4==3)return power(n,(p+1)/4,p);U q=p-1,z=2;unsigned s=0;while(!(q&1)){q>>=1;s++;}while(power(z,(p-1)/2,p)!=p-1)z++;U c=power(z,q,p),x=power(n,(q+1)/2,p),t=power(n,q,p);unsigned m=s;while(t!=1){unsigned i=1;U v=sq(t);while(v!=1){v=sq(v);i++;}U b=power(c,1ULL<<(m-i-1),p);x=mul(x,b);c=sq(b);t=mul(t,c);m=i;}return x;}
Pt norm(Pt q){if(!q.z)return q;U z=power(q.z,p-2,p),z2=sq(z);return{mul(q.x,z2),mul(q.y,mul(z,z2)),1};}
Pt scal(U k,U x,U y,U a){return k?norm(scalar(k,x,y,a)):Pt{0,1,0};}
bool eq(Pt a,Pt b){return a.z==b.z&&(!a.z||(a.x==b.x&&a.y==b.y));}
Pt endo(Pt P){if(!P.z)return P;U den[2]={sub(P.x,CM_X[0]),sub(P.x,CM_X[1])};if(!den[0]||!den[1])throw std::runtime_error("Kernel intersects prime subgroup");U prodinv=power(mul(den[0],den[1]),p-2,p);U invs[2]={mul(den[1],prodinv),mul(den[0],prodinv)};U x=P.x,derivative=1;for(unsigned i=0;i<2;i++){U inv=invs[i],inv2=sq(inv);x=add(x,add(mul(mul(2,CM_T[i]),inv),mul(mul(4,CM_Y2[i]),inv2)));derivative=sub(derivative,add(mul(mul(2,CM_T[i]),inv2),mul(mul(8,CM_Y2[i]),mul(inv2,inv))));}return {mul(x,CM_INV_C),mul(mul(P.y,derivative),CM_INV_CU),1};}
U hash64(U x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}
double cpu(){timespec t;clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t);return t.tv_sec+t.tv_nsec*1e-9;}
double wall(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
struct Job{U id,pair,role,fp,order,a,b,gx,gy,qx,qy,seed,repeat;};
std::mutex output_mutex;
std::map<U,std::vector<U>> crater_records;
void run(Job w,unsigned B,U limit,double seconds,std::string dir){
 const CMParameters& cp=all_cm_parameters.at(w.pair);CM_P=cp.p;CM_R=cp.r;CM_PAIR=cp.pair;CM_LAMBDA=cp.lambda;CM_INV_C=cp.inv_c;CM_INV_CU=cp.inv_cu;for(unsigned i=0;i<2;i++){CM_X[i]=cp.x[i];CM_T[i]=cp.t[i];CM_Y2[i]=cp.y2[i];}
 p=w.fp;r=w.order;field_muls=field_sqrs=group_doubles=group_mixed=scalar_calls=0;double c0=cpu(),t0=wall();U scans=0,keys_count=0,queries=0,hits=0,duplicates=0,direct=0;unsigned rank=0;U recovered=0;bool verified=false;std::string status="QUERY_LIMIT";
 if(B>16384||limit>=r){throw std::runtime_error("Unsupported base or noninjective budget");}
 constexpr unsigned C=4;if(B%C||w.pair!=CM_PAIR||p!=CM_P||r!=CM_R)throw std::runtime_error("Wrong CM pilot input");double orbit_c0=cpu();
 std::vector<Pt> base;std::vector<unsigned> columns;std::vector<U> coefficients;std::unordered_set<U> xs;unsigned L=B/C;U next_x=1;
 for(unsigned column=0;column<C;column++){
  Pt P;for(;;){U x=next_x++;scans++;if(x==w.gx)continue;U rhs=add(add(mul(sq(x),x),mul(w.a,x)),w.b);U y=root(rhs);if(y<p&&y){P={x,std::min(y,p-y),1};break;}}
  U coefficient=1;
  for(unsigned level=0;level<L;level++){
   if(!xs.insert(P.x).second)throw std::runtime_error("Orbit abscissa collision");
   if(sq(P.y)!=add(add(mul(sq(P.x),P.x),mul(w.a,P.x)),w.b))throw std::runtime_error("Invalid orbit point");
   base.push_back(P);columns.push_back(column);coefficients.push_back(coefficient);
   if(level+1<L){P=(w.role!=0)?scal(CM_LAMBDA,P.x,P.y,w.a):endo(P);coefficient=(__uint128_t)coefficient*CM_LAMBDA%r;if(P.y>p-P.y){P.y=p-P.y;coefficient=r-coefficient;}}
  }
 }
 double orbit_cpu=cpu()-orbit_c0;
 if(seconds<0){std::vector<U>record;for(unsigned i=0;i<B;i++)for(U v:{(U)i,(U)columns[i],coefficients[i],base[i].x,base[i].y})record.push_back(v);bool exact_checked=false;if(w.role!=1){std::lock_guard<std::mutex>lock(output_mutex);auto it=crater_records.find(w.pair);if(it==crater_records.end())crater_records.emplace(w.pair,record);else{if(it->second!=record)throw std::runtime_error("Crater orbit records differ");exact_checked=true;}}U digest=1469598103934665603ULL; for(unsigned i=0;i<B;i++)for(U v:{(U)i,(U)columns[i],coefficients[i],base[i].x,base[i].y})for(unsigned b=0;b<8;b++){digest^=(v>>(8*b))&255;digest*=1099511628211ULL;}std::lock_guard<std::mutex>lock(output_mutex);std::cout<<std::setprecision(12)<<"{\"experiment\":\"orbit-only\",\"job\":"<<w.id<<",\"pair\":"<<w.pair<<",\"role\":"<<w.role<<",\"repeat\":"<<w.repeat<<",\"orbit_cpu_seconds\":"<<orbit_cpu<<",\"same_crater_exact_checked\":"<<(exact_checked?"true":"false")<<",\"point_record_fnv64\":\""<<std::hex<<digest<<std::dec<<"\",\"points\":"<<B<<",\"field_muls\":"<<field_muls<<",\"field_sqrs\":"<<field_sqrs<<",\"group_doubles\":"<<group_doubles<<",\"group_mixed\":"<<group_mixed<<",\"scalar_calls\":"<<scalar_calls<<"}"<<std::endl;return;}

 std::ofstream fb(dir+"/factor-bases/job-"+std::to_string(w.id)+".csv");fb<<"point,column,coefficient,x,y\n";for(unsigned i=0;i<B;i++)fb<<i<<","<<columns[i]<<","<<coefficients[i]<<","<<base[i].x<<","<<base[i].y<<"\n";fb.close();
 size_t cap=1;while(cap<2ULL*B*B)cap<<=1;size_t mask=cap-1;std::vector<U> keys(cap);std::vector<uint32_t> values(cap);
 auto insert=[&](U x,uint32_t v){size_t h=hash64(x)&mask;while(keys[h]&&keys[h]!=x+1)h=(h+1)&mask;if(!keys[h]){keys[h]=x+1;values[h]=v;keys_count++;}};
 std::vector<U> prefix(B),inverses(B);
 for(unsigned i=0;i<B;i++){
  Pt d=norm(dbl(base[i],w.a));if(d.z)insert(d.x,(i<<16)|(i<<1));
  U product=1;for(unsigned j=i+1;j<B;j++){prefix[j]=product;product=mul(product,sub(base[j].x,base[i].x));}
  U inv=power(product,p-2,p);for(unsigned j=B;j-->i+1;){inverses[j]=mul(inv,prefix[j]);inv=mul(inv,sub(base[j].x,base[i].x));}
  for(unsigned j=i+1;j<B;j++)for(unsigned neg=0;neg<2;neg++){U y=neg?p-base[j].y:base[j].y;U slope=mul(sub(y,base[i].y),inverses[j]);U x=sub(sub(sq(slope),base[i].x),base[j].x);insert(x,(i<<16)|(j<<1)|neg);}
 }
 double setup_c=cpu()-c0,setup_t=wall()-t0,collection_c0=cpu(),collection_t0=wall();
 double pdp_wall=0, relation_wall=0, descent_wall=0, recovery_wall=0; std::map<unsigned,std::vector<U>> pivots;std::unordered_set<std::string> column_parts;
 std::ofstream relations(dir+"/relations/job-"+std::to_string(w.id)+".jsonl");
 auto modadd=[](U a,U b){return(a+b)%r;};auto modsub=[](U a,U b){return a>=b?a-b:r-(b-a);};auto modmul=[](U a,U b){return(U)((__uint128_t)a*b%r);};
 auto relation=[&](Pt R,U alpha){
  size_t h=hash64(R.x)&mask;while(keys[h]&&keys[h]!=R.x+1)h=(h+1)&mask;if(!keys[h])return;
  uint32_t code=values[h];unsigned i=code>>16,j=(code&65535)>>1;int si=1,sj=(code&1)?-1:1;
  double check_start=wall(); Pt S=norm(mixed(base[i],base[j].x,sj<0?p-base[j].y:base[j].y,w.a));
  if(!S.z||S.x!=R.x)throw std::runtime_error("Invalid table point");
  if(S.y!=R.y){if(sub(0,S.y)!=R.y)throw std::runtime_error("Invalid decomposition");si=-si;sj=-sj;}
  hits++;std::map<unsigned,U> part;part[columns[i]]=modadd(part[columns[i]],si>0?coefficients[i]:r-coefficients[i]);part[columns[j]]=modadd(part[columns[j]],sj>0?coefficients[j]:r-coefficients[j]);for(auto it=part.begin();it!=part.end();)if(!it->second)it=part.erase(it);else ++it;if(part.empty())throw std::runtime_error("Zero relation part");U reciprocal=power(part.begin()->second,r-2,r);std::string signature;for(auto [col,coef]:part)signature+=std::to_string(col)+":"+std::to_string(modmul(coef,reciprocal))+",";if(!column_parts.insert(signature).second){duplicates++;return;}
  Pt expected=norm(mixed(scal(alpha,w.gx,w.gy,w.a),w.qx,w.qy,w.a));if(!eq(expected,R))throw std::runtime_error("Wrong walk coefficient");
  relation_wall+=wall()-check_start; double descent_start=wall(); std::vector<U> row(C+2);for(auto [column,coef]:part)row[column]=coef;row[C]=r-1;row[C+1]=alpha;
  bool accepted=false;for(unsigned k=0;k<=C;k++)if(row[k]){auto it=pivots.find(k);if(it==pivots.end()){U inv=power(row[k],r-2,r);for(unsigned z=k;z<C+2;z++)row[z]=modmul(row[z],inv);pivots.emplace(k,row);rank++;accepted=true;break;}U factor=row[k];for(unsigned z=k;z<C+2;z++)row[z]=modsub(row[z],modmul(factor,it->second[z]));}
  if(!accepted&&row[C+1])throw std::runtime_error("Inconsistent linear system");
  relations<<"{\"alpha\":"<<alpha<<",\"beta\":1,\"i\":"<<i<<",\"j\":"<<j<<",\"si\":"<<si<<",\"sj\":"<<sj<<",\"Rx\":"<<R.x<<",\"Ry\":"<<R.y<<",\"rank\":"<<rank<<",\"independent\":"<<(accepted?"true":"false")<<"}\n";
  descent_wall+=wall()-descent_start; double recovery_start=wall(); auto pinned=pivots.find(C);if(pinned!=pivots.end()){recovered=pinned->second[C+1];verified=eq(scal(recovered,w.gx,w.gy,w.a),Pt{w.qx,w.qy,1});if(!verified)throw std::runtime_error("Invalid recovered scalar");status="PASS";} recovery_wall+=wall()-recovery_start;
 };
 // An injective shifted-generator stream: R=Q+[alpha]G. No scalar or rho result is an input.
 constexpr unsigned lanes=64;std::array<Pt,lanes> points;std::array<U,lanes> alphas;U start=w.seed%r;
 Pt first=norm(mixed(scal(start,w.gx,w.gy,w.a),w.qx,w.qy,w.a));points[0]=first;alphas[0]=start;
 for(unsigned i=1;i<lanes;i++){points[i]=norm(mixed(points[i-1],w.gx,w.gy,w.a));alphas[i]=(start+i)%r;}
 Pt jump=scal(lanes,w.gx,w.gy,w.a);std::array<U,lanes> denominators,products;
 while(queries<limit&&!verified){
  if(wall()-t0>=seconds){status="TIME_LIMIT";break;}
  unsigned count=std::min<U>(lanes,limit-queries);
  for(unsigned i=0;i<count;i++){if(points[i].z){double pdp_start=wall();relation(points[i],alphas[i]);pdp_wall+=wall()-pdp_start;}else direct++;queries++;if(verified)break;}
  if(verified||queries>=limit)break;
  U product=1;for(unsigned i=0;i<lanes;i++){denominators[i]=points[i].z?sub(jump.x,points[i].x):0;products[i]=product;if(denominators[i])product=mul(product,denominators[i]);}
  U inverse=power(product,p-2,p);for(unsigned i=lanes;i-->0;){U inv=0;if(denominators[i]){inv=mul(inverse,products[i]);inverse=mul(inverse,denominators[i]);}
   if(!denominators[i])points[i]=norm(mixed(points[i],jump.x,jump.y,w.a));else{U slope=mul(sub(jump.y,points[i].y),inv);U x=sub(sub(sq(slope),points[i].x),jump.x);U y=sub(mul(slope,sub(points[i].x,x)),points[i].y);points[i]={x,y,1};}alphas[i]=(alphas[i]+lanes)%r;
  }
 }
 relations.close();double collection_c=cpu()-collection_c0,collection_t=wall()-collection_t0,total_c=cpu()-c0,total_t=wall()-t0;
 std::lock_guard<std::mutex> lock(output_mutex);std::cout<<std::setprecision(12)<<"{\"field_muls\":"<<field_muls<<",\"field_sqrs\":"<<field_sqrs<<",\"group_doubles\":"<<group_doubles<<",\"group_mixed\":"<<group_mixed<<",\"scalar_calls\":"<<scalar_calls<<",\"online_start_monotonic_ns\":"<<(U)(collection_t0*1e9)<<",\"target_query_wall_seconds\":"<<collection_t-pdp_wall<<",\"target_PDP_wall_seconds\":"<<pdp_wall-relation_wall-descent_wall-recovery_wall<<",\"target_relation_check_wall_seconds\":"<<relation_wall<<",\"target_descent_wall_seconds\":"<<descent_wall<<",\"target_recovery_check_wall_seconds\":"<<recovery_wall<<",\"method\":\"CM-orbit-folded-MITM2\",\"columns\":"<<C<<",\"orbit_length\":"<<L<<",\"orbit_cpu_seconds\":"<<orbit_cpu<<",\"job\":"<<w.id<<",\"pair\":"<<w.pair<<",\"role\":"<<w.role<<",\"repeat\":"<<w.repeat<<",\"status\":\""<<status<<"\",\"B\":"<<B<<",\"unknowns\":"<<C+1<<",\"p\":"<<p<<",\"r\":"<<r<<",\"queries\":"<<queries<<",\"query_limit\":"<<limit<<",\"hits\":"<<hits<<",\"rejected_duplicate_columns\":"<<duplicates<<",\"direct_infinity_skipped\":"<<direct<<",\"rank\":"<<rank<<",\"pair_table_keys\":"<<keys_count<<",\"table_bytes\":"<<cap*12<<",\"factor_base_scans\":"<<scans<<",\"setup_cpu_seconds\":"<<setup_c<<",\"setup_wall_seconds\":"<<setup_t<<",\"collection_cpu_seconds\":"<<collection_c<<",\"collection_wall_seconds\":"<<collection_t<<",\"cpu_seconds\":"<<total_c<<",\"wall_seconds\":"<<total_t<<",\"recovered\":"<<recovered<<",\"verified\":"<<(verified?"true":"false")<<"}"<<std::endl;
}
int main(int argc,char**argv){if(argc!=8){std::cerr<<"usage: ic workloads B queries wall_seconds threads output_dir cm_params\n";return 2;}std::ifstream params(argv[7]);CMParameters cp;while(params>>cp.p>>cp.r>>cp.pair>>cp.lambda>>cp.inv_c>>cp.inv_cu>>cp.x[0]>>cp.x[1]>>cp.t[0]>>cp.t[1]>>cp.y2[0]>>cp.y2[1]){if(!all_cm_parameters.emplace(cp.pair,cp).second)throw std::runtime_error("Duplicate CM parameters");}if(!params.eof()||all_cm_parameters.empty())throw std::runtime_error("Bad CM parameters");std::ifstream input(argv[1]);std::string line;getline(input,line);std::vector<Job> jobs;while(getline(input,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);Job w;s>>w.id>>w.pair>>w.role>>w.fp>>w.order>>w.a>>w.b>>w.gx>>w.gy>>w.qx>>w.qy>>w.seed>>w.repeat;if(!s)throw std::runtime_error("bad workload");jobs.push_back(w);}if(std::stod(argv[4])<0){auto originals=jobs;jobs.clear();std::mt19937_64 rng(20261015);for(U rep=0;rep<20;rep++){std::shuffle(originals.begin(),originals.end(),rng);for(auto w:originals){w.repeat=rep;jobs.push_back(w);}}}unsigned B=std::stoul(argv[2]);U limit=std::stoull(argv[3]);double seconds=std::stod(argv[4]);std::atomic<size_t> next{0};std::vector<std::thread> workers;for(unsigned i=0;i<std::stoul(argv[5]);i++)workers.emplace_back([&]{for(;;){size_t j=next++;if(j>=jobs.size())break;run(jobs[j],B,limit,seconds,argv[6]);}});for(auto&t:workers)t.join();}
