#pragma once
#include "ec.hpp"
#include <tuple>
namespace Algebra {
struct Count{U multiply=0,add_sub=0,inverse=0,reductions=0,pairs_created=0,pairs_processed=0,product_skips=0,basis_insertions=0;int max_actual_degree=0,max_sugar_degree=0;json wire()const{return {{"field_multiply_including_inverse_power",multiply},{"field_add_sub",add_sub},{"field_inverse",inverse},{"leading_reductions",reductions},{"pairs_created",pairs_created},{"pairs_processed",pairs_processed},{"product_criterion_skips",product_skips},{"basis_insertions",basis_insertions},{"max_actual_polynomial_degree",max_actual_degree},{"max_sugar_degree",max_sugar_degree}};}};
inline thread_local Count* counts=nullptr;
inline U plus(U a,U b){if(counts)counts->add_sub++;return add(a,b);}
inline U minus(U a,U b){if(counts)counts->add_sub++;return sub(a,b);}
inline U times(U a,U b){if(counts)counts->multiply++;return mul(a,b);}
inline U inverse(U a){need(a!=0,"polynomial inverse zero");if(counts)counts->inverse++;U k=PRIME-2,r=1;for(;k;k>>=1,a=times(a,a))if(k&1)r=times(r,a);return r;}
struct Mon{int x=0,y=0;int degree()const{return x+y;}bool operator==(Mon b)const{return x==b.x&&y==b.y;}};
struct Order{bool operator()(Mon a,Mon b)const{return a.degree()!=b.degree()?a.degree()>b.degree():a.x>b.x;}};
using Poly=std::map<Mon,U,Order>;
inline int degree(const Poly&p){return p.empty()?-1:p.begin()->first.degree();}
inline void note(const Poly&p){if(counts)counts->max_actual_degree=std::max(counts->max_actual_degree,degree(p));}
inline void term(Poly&p,Mon m,U c){if(!c)return;auto it=p.find(m);if(it==p.end())p.emplace(m,c);else{U v=plus(it->second,c);if(v)it->second=v;else p.erase(it);}note(p);}
inline Poly scale(const Poly&p,Mon shift,U c){Poly r;for(auto&[m,v]:p)term(r,{m.x+shift.x,m.y+shift.y},times(v,c));return r;}
inline Poly monic(Poly p){if(p.empty()||p.begin()->second==1)return p;U v=inverse(p.begin()->second);for(auto&[m,c]:p){(void)m;c=times(c,v);}return p;}
inline bool divides(Mon a,Mon b){return a.x<=b.x&&a.y<=b.y;}
inline Poly reduce(Poly p,const std::vector<Poly>&basis){Poly rest;while(!p.empty()){auto [m,c]=*p.begin();bool did=false;for(const auto&b:basis){if(b.empty())continue;Mon lm=b.begin()->first;if(!divides(lm,m))continue;if(counts)counts->reductions++;U factor=b.begin()->second==1?c:times(c,inverse(b.begin()->second));Mon shift{m.x-lm.x,m.y-lm.y};for(auto&[bm,bc]:b)term(p,{bm.x+shift.x,bm.y+shift.y},minus(0,times(factor,bc)));did=true;break;}if(!did){term(rest,m,c);p.erase(p.begin());}}return rest;}
inline Mon lcm(Mon a,Mon b){return{std::max(a.x,b.x),std::max(a.y,b.y)};}
inline Poly spol(const Poly&a,const Poly&b){Mon aa=a.begin()->first,bb=b.begin()->first,L=lcm(aa,bb);auto r=scale(a,{L.x-aa.x,L.y-aa.y},inverse(a.begin()->second));auto q=scale(b,{L.x-bb.x,L.y-bb.y},inverse(b.begin()->second));for(auto&[m,v]:q)term(r,m,minus(0,v));return r;}
struct Result{std::vector<Poly>basis;json events=json::array();bool unit=false;};
inline Result groebner(const std::vector<Poly>&input,bool trace){Result r;std::vector<int>sugars;struct Pair{int i,j,sugar;Mon L;};std::vector<Pair>pairs;
 auto insert=[&](Poly p,int sugar){int id=int(r.basis.size());p=monic(std::move(p));r.basis.push_back(std::move(p));sugars.push_back(sugar);if(counts)counts->basis_insertions++;for(int j=0;j<id;j++){Mon L=lcm(r.basis[j].begin()->first,r.basis[id].begin()->first);int d=std::max(sugars[j]+L.degree()-degree(r.basis[j]),sugar+L.degree()-degree(r.basis[id]));pairs.push_back({j,id,d,L});if(counts)counts->pairs_created++;}};
 for(auto p:input){note(p);if(!p.empty())insert(p,degree(p));}
 while(!pairs.empty()){auto it=std::min_element(pairs.begin(),pairs.end(),[](const Pair&a,const Pair&b){return std::tie(a.sugar,a.L.x,a.L.y,a.i,a.j)<std::tie(b.sugar,b.L.x,b.L.y,b.i,b.j);});Pair pair=*it;pairs.erase(it);if(counts){counts->pairs_processed++;counts->max_sugar_degree=std::max(counts->max_sugar_degree,pair.sugar);}auto ma=r.basis[pair.i].begin()->first,mb=r.basis[pair.j].begin()->first;if(std::min(ma.x,mb.x)==0&&std::min(ma.y,mb.y)==0){if(counts)counts->product_skips++;if(trace)r.events.push_back({{"i",pair.i},{"j",pair.j},{"sugar",pair.sugar},{"product_skip",true}});continue;}
 auto reduced=reduce(spol(r.basis[pair.i],r.basis[pair.j]),r.basis);if(trace)r.events.push_back({{"i",pair.i},{"j",pair.j},{"sugar",pair.sugar},{"remainder_degree",degree(reduced)},{"remainder_terms",reduced.size()}});if(reduced.empty())continue;if(degree(reduced)==0){r.basis={Poly{{Mon{0,0},1}}};r.unit=true;return r;}insert(std::move(reduced),pair.sugar);
 need(r.basis.size()<256&&pairs.size()<50000,"GB safety cap");}return r;
}
inline Poly factor_poly(const std::vector<Point>&base,bool y){Poly p{{Mon{0,0},1}};for(Point P:base){Poly q;for(auto&[m,c]:p){term(q,{m.x+(y?0:1),m.y+(y?1:0)},c);term(q,m,minus(0,times(c,P.x)));}p=std::move(q);}return p;}
inline Poly f3(U z,U a,U b){Poly f;term(f,{2,2},1);U twoz=times(2,z),z2=times(z,z),twoa=times(2,a);term(f,{2,1},minus(0,twoz));term(f,{1,2},minus(0,twoz));term(f,{2,0},z2);term(f,{0,2},z2);term(f,{1,1},minus(0,plus(times(2,z2),twoa)));U linear=minus(0,plus(times(twoa,z),times(4,b)));term(f,{1,0},linear);term(f,{0,1},linear);term(f,{0,0},minus(times(a,a),times(times(4,b),z)));return f;}
inline U evaluate(const Poly&p,U x,U y){U r=0;for(auto&[m,c]:p)r=add(r,mul(c,mul(power(x,U(m.x)),power(y,U(m.y)))));return r;}
inline json wire(const Poly&p){json j=json::array();for(auto&[m,c]:p)j.push_back({m.x,m.y,std::to_string(c)});return j;}
inline json wire(const std::vector<Poly>&p){json j=json::array();for(const auto&a:p)j.push_back(wire(a));return j;}
inline int dreg_top_forms(int B){for(int d=0;d<2*B+2;d++){int count=0;for(int x=0;x<=d;x++){int y=d-x;if(x<B&&y<B&&!(x>=2&&y>=2))count++;}if(!count)return d;}throw std::runtime_error("Hilbert regularity");}
}
