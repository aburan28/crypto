#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using U = unsigned __int128;
struct J {
  std::string scalar, type;
  std::vector<J> array;
  std::map<std::string,J> object;
  const J& at(const std::string& k) const { return object.at(k); }
  const J& at(size_t n) const { return array.at(n); }
};
struct Parser {
  std::string s; size_t i=0;
  void ws(){while(i<s.size() && std::isspace((unsigned char)s[i]))++i;}
  std::string str(){assert(s.at(i++)=='"');std::string o;while(s.at(i)!='"'){char c=s.at(i++);if(c=='\\'){c=s.at(i++);if(c=='n')c='\n';else if(c=='t')c='\t';else if(c=='r')c='\r';}o+=c;}++i;return o;}
  J value(){ws();J j;char c=s.at(i);if(c=='{'){j.type="object";++i;ws();if(s.at(i)=='}'){++i;return j;}do{ws();auto k=str();ws();assert(s.at(i++)==':');j.object.emplace(k,value());ws();c=s.at(i++);}while(c==',');assert(c=='}');}else if(c=='['){j.type="array";++i;ws();if(s.at(i)==']'){++i;return j;}do{j.array.push_back(value());ws();c=s.at(i++);}while(c==',');assert(c==']');}else if(c=='"'){j.type="string";j.scalar=str();}else{size_t b=i;while(i<s.size() && s[i]!=',' && s[i]!=']' && s[i]!='}' && !std::isspace((unsigned char)s[i]))++i;j.scalar=s.substr(b,i-b);j.type="literal";}return j;}
};
J read_json(const std::string& path){std::ifstream f(path);if(!f)throw std::runtime_error("Cannot read "+path);Parser p;p.s=std::string(std::istreambuf_iterator<char>(f),{});auto j=p.value();p.ws();assert(p.i==p.s.size());return j;}
U decimal(const std::string& s){U v=0;assert(!s.empty());for(char c:s){assert(c>='0'&&c<='9');v=v*10+(c-'0');}return v;}
U num(const J& j){return decimal(j.scalar);}
double real(const J& j){return std::stod(j.scalar);}
std::string dec(U x){if(!x)return "0";std::string s;while(x){s+=char('0'+x%10);x/=10;}std::reverse(s.begin(),s.end());return s;}
bool truth(const J&j){return j.scalar=="true";}
int degree(U x){if(!x)return -1;int d=0;while(x>>=1)++d;return d;}
constexpr U TOP=U(1)<<83;
constexpr U LOW=(U(1)<<7)|(U(1)<<4)|(U(1)<<2)|1;
constexpr U MOD=TOP|LOW;
U xt(U x){x<<=1;if(x&TOP)x^=MOD;return x;}
U fm(U x,U y){U z=0;while(y){if(y&1)z^=x;y>>=1;x=xt(x);}return z;}
U fs(U x){return fm(x,x);}
U shiftmod(U x,int n){while(n--)x=xt(x);return x;}
U fi(U a){assert(a);U u=a,v=MOD,g1=1,g2=0;while(u!=1){int j=degree(u)-degree(v);if(j<0){std::swap(u,v);std::swap(g1,g2);j=-j;}u^=v<<j;g1^=shiftmod(g2,j);}assert(fm(a,g1)==1);return g1;}
U rem(U a,U b){int db=degree(b);while(a && degree(a)>=db)a^=b<<(degree(a)-db);return a;}
U pgcd(U a,U b){while(b){U t=rem(a,b);a=b;b=t;}return a;}
struct P{U x=0,y=0;bool inf=true;};
bool operator==(P a,P b){return a.inf==b.inf&&(a.inf||(a.x==b.x&&a.y==b.y));}
P point(const J&j){return {num(j.at(0)),num(j.at(1)),false};}
P neg(P p){if(!p.inf)p.y^=p.x;return p;}
P tau(P p){if(!p.inf){p.x=fs(p.x);p.y=fs(p.y);}return p;}
bool on(P p){return p.inf || (p.x<TOP&&p.y<TOP&&(fs(p.y)^fm(p.x,p.y))==(fm(fs(p.x),p.x)^fs(p.x)^1));}
P add(P p,P q){if(p.inf)return q;if(q.inf)return p;U m;if(p.x==q.x){if(p.y!=q.y || p.x==0)return {};m=p.x^fm(p.y,fi(p.x));U x=fs(m)^m^1;return {x,fs(p.x)^fm(m^1,x),false};}m=fm(p.y^q.y,fi(p.x^q.x));U x=fs(m)^m^p.x^q.x^1;return {x,fm(m,p.x^x)^x^p.y,false};}
P mul(U k,P p){P z;while(k){if(k&1)z=add(z,p);p=add(p,p);k>>=1;}return z;}
uint64_t mm(uint64_t x,uint64_t y,uint64_t n){return uint64_t(U(x)*y%n);}
uint64_t power(uint64_t x,uint64_t e,uint64_t n){uint64_t z=1;while(e){if(e&1)z=mm(z,x,n);x=mm(x,x,n);e>>=1;}return z;}
bool prime(uint64_t n){if(n<2)return false;for(uint64_t p:{2,3,5,7,11,13,17,19,23,29,31,37}){if(n==p)return true;if(n%p==0)return false;}uint64_t d=n-1;int s=0;while(!(d&1)){d>>=1;++s;}for(uint64_t a:{2ULL,325ULL,9375ULL,28178ULL,450775ULL,9780504ULL,1795265022ULL}){if(a%n==0)continue;auto x=power(a%n,d,n);if(x==1||x==n-1)continue;bool ok=false;for(int i=1;i<s;++i){x=mm(x,x,n);if(x==n-1){ok=true;break;}}if(!ok)return false;}return true;}
void require(bool b,const char* msg){if(!b)throw std::runtime_error(msg);}
void check_modulus(const J& j){const auto& terms=j.at("field_modulus_low_terms").array;require(terms.size()==4,"Wrong modulus term count");const unsigned expected[4]={0,2,4,7};for(size_t i=0;i<4;++i)require(num(terms[i])==expected[i],"Wrong modulus term");}
int main(int argc,char**argv){try{
 require(argc==3,"usage: check SOURCE_DIR OUT_JSON");std::string source=argv[1];auto fixture=read_json(source+"/frozen/fixture.json");
 require(num(fixture.at("n"))==83 && num(fixture.at("a"))==1,"Wrong fixture curve");
 check_modulus(fixture);
 U x=2;for(int i=0;i<83;++i)x=fs(x);require(x==2,"Rabin final identity failed");require(pgcd(fs(2)^2,MOD)==1,"Rabin gcd failed");
 uint64_t r=uint64_t(num(fixture.at("subgroup_order")));require(prime(r),"r is composite");
 __int128 t0=2,t1=1;for(int n=2;n<=83;++n){auto t=t1-2*t0;t0=t1;t1=t;}require(t1==6151469093347LL,"Recurrence trace mismatch");
 U order=TOP+1-U(t1);require(order==num(fixture.at("curve_order")),"Curve order mismatch");require(order==num(fixture.at("cofactor"))*r,"Order/cofactor mismatch");
 P g=point(fixture.at("generator")),q=point(fixture.at("public_target"));require(on(g)&&on(q),"Fixture off curve");require(!g.inf&&!q.inf,"Identity fixture");require(mul(r,g).inf&&mul(r,q).inf,"Subgroup membership failed");U sidecar=num(fixture.at("fixture_scalar_validation_only"));require(mul(sidecar,g)==q,"Sidecar point replay failed");
 const std::string first="N83A1K600W12We202610061R21R1";auto rho0=read_json(source+"/runs/"+first+"/rho.jsonl");uint64_t lambda=uint64_t(num(rho0.at("lambda")));
 require((mm(lambda,lambda,r)+r-lambda+2)%r==0,"Lambda characteristic equation failed");require(lambda!=1&&power(lambda,83,r)==1,"Lambda order failed");require(mul(lambda,g)==tau(g),"Tau eigenvalue failed");
 auto base=read_json(source+"/base_n83_K600.jsonl");require(num(base.at("subgroup_order"))==r&&num(base.at("a"))==1&&num(base.at("n"))==83,"Base curve mismatch");
 check_modulus(base);
 const auto& raw=base.at("factor_base_point_coordinates").array;const auto& labels=base.at("factor_base_point_labels").array;const auto& reps=base.at("factor_base_representatives").array;
 require(raw.size()==99600&&labels.size()==99600&&reps.size()==600,"Base dimensions failed");require(num(base.at("factor_base_points"))==99600&&num(base.at("orbit_columns"))==600,"Base metadata failed");
 std::vector<P> points;points.reserve(raw.size());std::set<std::pair<U,U>> distinct;
 for(size_t c=0;c<600;++c){P rep=point(reps[c]);require(on(rep)&&!rep.inf&&mul(r,rep).inf,"Base representative subgroup failed");P cur=rep;uint64_t coefficient=1;std::set<std::pair<U,U>> local;for(int k=0;k<83;++k){for(int sign=0;sign<2;++sign){size_t i=166*c+2*k+sign;P expected=sign?neg(cur):cur;P stored=point(raw[i]);require(on(stored)&&stored==expected,"Base signed Frobenius coordinate failed");require(num(labels[i].at(0))==c&&num(labels[i].at(1))==(sign?r-coefficient:coefficient),"Base scalar coefficient failed");require(local.insert({stored.x,stored.y}).second,"Short base orbit");require(distinct.insert({stored.x,stored.y}).second,"Duplicate across base orbits");points.push_back(stored);}cur=tau(cur);coefficient=mm(coefficient,lambda,r);}require(cur==rep&&coefficient==1,"Base orbit closure failed");}
 std::ostringstream rows;rows<<std::setprecision(17);size_t checked_runs=0;
 for(int run=1;run<=3;++run){std::string id="N83A1K600W12We20261006"+std::to_string(run)+"R2"+std::to_string(run)+"R"+std::to_string(run);std::string path=source+"/runs/"+id;auto ic=read_json(path+"/ic.jsonl"),rho=read_json(path+"/rho.jsonl"),sum=read_json(path+"/ic_summary.jsonl"),meta=read_json(path+"/run.json");
 check_modulus(rho);require(num(meta.at("ic_return_code"))==0&&num(meta.at("rho_return_code"))==0&&num(ic.at("exit_code"))==0,"Recorded producer failure");
 require(point(ic.at("generator"))==g&&point(rho.at("generator"))==g,"Generator mismatch between arms");require(point(ic.at("published_q"))==q&&point(rho.at("published_q"))==q&&point(ic.at("target"))==q,"Target mismatch");
 require(ic.at("published_fixture_scalar").scalar=="null"&&rho.at("published_fixture_scalar").scalar=="null","Scalar present in producer input");require(!truth(meta.at("known_answer_sent_to_ic"))&&!truth(meta.at("known_answer_sent_to_rho")),"Metadata scalar input flag");
 require(truth(ic.at("group_verified"))&&truth(rho.at("verified")),"Producer verification flag false");require(num(ic.at("recovered_scalar"))==sidecar&&num(rho.at("recovered_fixture_scalar"))==sidecar,"Producer scalar mismatch");require(mul(num(ic.at("recovered_scalar")),g)==q&&mul(num(rho.at("recovered_fixture_scalar")),g)==q,"Native scalar replay failed");
 const auto& idx=ic.at("point_indices").array;require(idx.size()==4,"Not four summands");P total;for(size_t i=0;i<4;++i){size_t at=size_t(num(idx[i]));require(at<points.size(),"Witness index out of range");require(points[at].x==num(ic.at("x_codes").at(i)),"Witness x code mismatch");total=add(total,points[at]);}require(total==q,"Four-point witness does not sum to Q");P p01=add(points[size_t(num(idx[0]))],points[size_t(num(idx[1]))]);P p23=add(points[size_t(num(idx[2]))],points[size_t(num(idx[3]))]);require(p01.x==num(ic.at("pinned_intermediates").at(0))&&p23.x==num(ic.at("pinned_intermediates").at(1)),"Intermediate roots mismatch");
 require(num(sum.at("rank"))==600&&num(sum.at("rank_relations"))==600&&num(sum.at("rank_failures"))==0&&num(sum.at("targets_solved"))==1&&num(sum.at("targets_failed"))==0,"Summary rank/target status mismatch");require(num(rho.at("lambda"))==lambda&&num(rho.at("automorphism_size"))==166,"Rho symmetry mismatch");
 double phases=real(ic.at("target_query_ms"))+real(ic.at("target_pdp_ms"))+real(ic.at("target_relation_check_ms"))+real(ic.at("target_descent_ms"))+real(ic.at("target_recovery_check_ms"));double ict=real(ic.at("target_ms"));double rhot=real(rho.at("total_ms"));double rhophases=real(rho.at("setup_ms"))+real(rho.at("walk_ms"))+real(rho.at("validation_ms"));
 require(std::abs(phases-real(ic.at("target_phase_sum_ms")))<1e-8,"IC recorded phase sum arithmetic mismatch");require(std::abs(rhot-rhophases)<1e-7,"Rho phase sum arithmetic mismatch");
 if(checked_runs++){rows<<",\n";}
 rows<<"    {\"run_id\":\""<<id<<"\",\"four_point_witness_native_replay\":true,\"scalar_native_replay_both_arms\":true,\"ic_target_ms\":"<<ict<<",\"ic_five_phase_sum_ms\":"<<phases<<",\"ic_unassigned_elapsed_ms\":"<<(ict-phases)<<",\"rho_recorded_total_ms\":"<<rhot<<",\"rho_setup_ms_unresolved_target_dependence\":"<<real(rho.at("setup_ms"))<<",\"rho_internal_total_over_ic_target\":"<<(rhot/ict)<<",\"rho_walk_steps\":"<<rho.at("walk_steps").scalar<<",\"ic_probes\":"<<ic.at("probes").scalar<<"}";
 }
 std::ofstream out(argv[2]);out<<"{\n  \"kind\":\"native_independent_archive_correctness_check\",\n  \"status\":\"PASS\",\n  \"curve_model_slug\":\"icv1-f2m83-t6151469093347-cdcc5432\",\n  \"scope\":\"Standalone C++ GF(2^83) arithmetic and archive witness/metadata replay; no retiming or solver re-execution\",\n  \"field_rabin_irreducibility\":true,\n  \"subgroup_deterministic_miller_rabin_64bit_prime\":true,\n  \"frobenius_recurrence_curve_order\":\""<<dec(order)<<"\",\n  \"subgroup_order\":\""<<r<<"\",\n  \"cofactor\":\""<<fixture.at("cofactor").scalar<<"\",\n  \"generator_target_on_curve_and_in_subgroup\":true,\n  \"validation_sidecar_scalar_replay\":true,\n  \"tau_eigenvalue\":\""<<lambda<<"\",\n  \"tau_eigenvalue_exact_order\":83,\n  \"signed_frobenius_orbit_size\":166,\n  \"factor_base_representatives_subgroup_verified\":600,\n  \"factor_base_distinct_signed_orbit_points_verified\":99600,\n  \"factor_base_labels_verified\":99600,\n  \"runs\":[\n"<<rows.str()<<"\n  ],\n  \"not_verified\":[\"600 independent rank equations and recovered factor logs: witnesses/log database not retained in these summary receipts\",\"Timing isolation and matched enforced memory/resource limits\",\"Target dependence of rho setup\",\"A new unknown-scalar recovery or performance comparison\",\"Exact canonical candidate/workload identity manifests and five-exclusive-phase timing contract\"]\n}\n";out.close();std::cout<<"PASS: field, subgroup, generator/target, 600 reps, 99600 orbit points/labels, 3 four-point witnesses and 6 scalar answer replays\n";return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
