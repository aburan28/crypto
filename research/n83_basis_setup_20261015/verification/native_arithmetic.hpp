#pragma once
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
  J value(){ws();J j;char c=s.at(i);if(c=='{'){j.type="object";++i;ws();if(s.at(i)=='}'){++i;return j;}do{ws();auto k=str();ws();assert(s.at(i++)==':');auto inserted=j.object.emplace(k,value());if(!inserted.second)throw std::runtime_error("Duplicate JSON key");ws();c=s.at(i++);}while(c==',');assert(c=='}');}else if(c=='['){j.type="array";++i;ws();if(s.at(i)==']'){++i;return j;}do{j.array.push_back(value());ws();c=s.at(i++);}while(c==',');assert(c==']');}else if(c=='"'){j.type="string";j.scalar=str();}else{size_t b=i;while(i<s.size() && s[i]!=',' && s[i]!=']' && s[i]!='}' && !std::isspace((unsigned char)s[i]))++i;j.scalar=s.substr(b,i-b);j.type="literal";}return j;}
};
J read_json(const std::string& path){std::ifstream f(path);if(!f)throw std::runtime_error("Cannot read "+path);Parser p;p.s=std::string(std::istreambuf_iterator<char>(f),{});auto j=p.value();p.ws();assert(p.i==p.s.size());return j;}
U decimal(const std::string& s){U v=0;assert(!s.empty());for(char c:s){assert(c>='0'&&c<='9');if(v>(~U(0)-U(c-'0'))/10)throw std::runtime_error("Decimal overflow");v=v*10+(c-'0');}return v;}
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
