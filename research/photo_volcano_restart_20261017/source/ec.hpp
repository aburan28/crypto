#pragma once
#include "json.hpp"
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <queue>
#include <stdexcept>
#include <sys/random.h>
#include <unistd.h>
using json=nlohmann::json;using U=uint64_t;using W=__uint128_t;
constexpr U PRIME=72057980207847001ULL,ORDER=18014494963342201ULL,TOTAL=72057979853368804ULL;
inline void need(bool b,const std::string&s){if(!b)throw std::runtime_error(s);}
inline U val(const json&j){return j.is_string()?std::stoull(j.get<std::string>()):j.get<U>();}
inline U add(U a,U b,U p=PRIME){return a>=p-b?a-(p-b):a+b;}
inline U sub(U a,U b,U p=PRIME){return a>=b?a-b:p-(b-a);}
inline U mul(U a,U b,U p=PRIME){return U(W(a)*b%p);}
inline U power(U a,U n,U p=PRIME){U r=1;for(;n;n>>=1,a=mul(a,a,p))if(n&1)r=mul(r,a,p);return r;}
inline U inv(U a,U p=PRIME){need(a!=0,"inverse zero");return power(a,p-2,p);}
inline bool prime(U n){if(n<2)return false;for(U a:{2,3,5,7,11,13,17,19,23,29,31,37})if(n%a==0)return n==a;U d=n-1;int s=0;while(!(d&1)){s++;d>>=1;}for(U a:{2ULL,325ULL,9375ULL,28178ULL,450775ULL,9780504ULL,1795265022ULL}){if(a%n==0)continue;U x=power(a%n,d,n);if(x==1||x==n-1)continue;bool ok=false;for(int i=1;i<s;i++){x=mul(x,x,n);if(x==n-1){ok=true;break;}}if(!ok)return false;}return true;}
inline bool sqrtmod(U a,U&x){if(!a){x=0;return true;}if(power(a,(PRIME-1)/2)!=1)return false;U q=PRIME-1;int s=0;while(!(q&1)){q>>=1;s++;}U z=2;while(power(z,(PRIME-1)/2)!=PRIME-1)z++;U c=power(z,q),t=power(a,q);x=power(a,(q+1)/2);int m=s;while(t!=1){int i=1;U v=mul(t,t);while(v!=1){v=mul(v,v);i++;need(i<m,"Tonelli failed");}U b=power(c,U(1)<<(m-i-1));x=mul(x,b);c=mul(b,b);t=mul(t,c);m=i;}x=std::min(x,PRIME-x);return true;}
struct Point{U x=0,y=0;bool infinity=true;bool operator==(Point p)const{return infinity==p.infinity&&(infinity||(x==p.x&&y==p.y));}};
inline Point point(const json&j){return j.is_null()?Point{}:Point{val(j[0]),val(j[1]),false};}
inline json wire(Point p){return p.infinity?json(nullptr):json::array({std::to_string(p.x),std::to_string(p.y)});}
struct EC{U a,b;bool on(Point p)const{return p.infinity||mul(p.y,p.y)==add(add(mul(mul(p.x,p.x),p.x),mul(a,p.x)),b);}
 Point neg(Point p)const{if(!p.infinity&&p.y)p.y=PRIME-p.y;return p;}
 Point op(Point p,Point q)const{if(p.infinity)return q;if(q.infinity)return p;if(p.x==q.x&&add(p.y,q.y)==0)return{};U t=p.x==q.x?mul(add(mul(3,mul(p.x,p.x)),a),inv(mul(2,p.y))):mul(sub(q.y,p.y),inv(sub(q.x,p.x)));U x=sub(sub(mul(t,t),p.x),q.x);return{x,sub(mul(t,sub(p.x,x)),p.y),false};}
 Point scalar(U k,Point p)const{Point r;while(k){if(k&1)r=op(r,p);p=op(p,p);k>>=1;}return r;}
};
inline Point map_point(Point p,const json&e){if(p.infinity)return{};U d=sub(p.x,val(e["kernel_x"]));if(!d)return{};U di=inv(d),v=val(e["v"]),u=val(e["u"]);U x=add(p.x,mul(v,di)),y=mul(p.y,sub(1,mul(v,mul(di,di))));return{mul(mul(u,u),x),mul(mul(mul(u,u),u),y),false};}
inline json load(const std::string&f){std::ifstream s(f);need(bool(s),"read "+f);json j;s>>j;return j;}
inline void save(const std::string&f,const json&j){std::ofstream s(f);need(bool(s),"write "+f);s<<j.dump(2)<<'\n';}
inline U ns(){return U(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count());}
struct Random{std::ofstream log;U index=0;explicit Random(const std::string&f):log(f){need(bool(log),"random journal");}U below(U bound,const std::string&scope){need(bound>0,"bound zero");U threshold=U(-bound)%bound;for(;;){U x;size_t n=0;while(n<8){ssize_t k=getrandom(reinterpret_cast<char*>(&x)+n,8-n,0);if(k<0&&errno==EINTR)continue;need(k>0,"getrandom failure");n+=size_t(k);}bool accepted=x>=threshold;log<<json({{"index",index++},{"scope",scope},{"bound",std::to_string(bound)},{"raw_word",std::to_string(x)},{"threshold",std::to_string(threshold)},{"accepted",accepted},{"value",accepted?json(std::to_string(x%bound)):json(nullptr)}}).dump()<<'\n';if(accepted)return x%bound;}}};
inline Point sample(const EC&e,Random&r,const std::string&scope){for(;;){U x=r.below(PRIME,scope+"/x"),y;if(!sqrtmod(add(add(mul(mul(x,x),x),mul(e.a,x)),e.b),y))continue;if(r.below(2,scope+"/sign")&&y)y=PRIME-y;Point p=e.scalar(4,{x,y,false});if(p.infinity)continue;need(e.scalar(ORDER,p).infinity,"sample order");return p;}}
