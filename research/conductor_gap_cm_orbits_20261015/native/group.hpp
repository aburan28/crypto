#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <chrono>
using U=uint64_t;
thread_local U p=0,r=0;
thread_local U field_muls=0, field_sqrs=0, group_doubles=0, group_mixed=0, scalar_calls=0;
U add(U a,U b){return (a+b)%p;}
U sub(U a,U b){return a>=b?a-b:p-(b-a);}
U mul(U a,U b){++field_muls;return (__uint128_t(a)*b)%p;}
U sq(U a){++field_sqrs;return (__uint128_t(a)*a)%p;}
struct Pt{U x,y,z;};
Pt dbl(Pt q,U a){
 ++group_doubles;
 if(!q.z||!q.y)return {0,1,0};
 U xx=sq(q.x),yy=sq(q.y),yyyy=sq(yy),zz=sq(q.z);
 U s=mul(2,sub(sub(sq(add(q.x,yy)),xx),yyyy));
 U m=add(mul(3,xx),mul(a,sq(zz)));
 U x=sub(sq(m),mul(2,s));
 U y=sub(mul(m,sub(s,x)),mul(8,yyyy));
 U z=sub(sub(sq(add(q.y,q.z)),yy),zz);
 return {x,y,z};
}
Pt mixed(Pt q,U x,U y,U a){
 ++group_mixed;
 if(!q.z)return {x,y,1};
 U zz=sq(q.z),u=mul(x,zz),s=mul(y,mul(q.z,zz));
 U h=sub(u,q.x),rr=mul(2,sub(s,q.y));
 if(!h)return !rr?dbl(q,a):Pt{0,1,0};
 U hh=sq(h),ii=mul(4,hh),jj=mul(h,ii),v=mul(q.x,ii);
 U nx=sub(sub(sq(rr),jj),mul(2,v));
 U ny=sub(mul(rr,sub(v,nx)),mul(2,mul(q.y,jj)));
 U nz=sub(sub(sq(add(q.z,h)),zz),hh);
 return {nx,ny,nz};
}
Pt scalar(U k,U x,U y,U a){
 ++scalar_calls; if(!k)return {0,1,0};
 Pt q{0,1,0};
 for(int i=63-__builtin_clzll(k);i>=0;i--){q=dbl(q,a);if((k>>i)&1)q=mixed(q,x,y,a);}
 return q;
}
