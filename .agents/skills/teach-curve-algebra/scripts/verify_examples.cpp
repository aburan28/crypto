// Exact teaching fixtures; no discrete-log solver or performance benchmark.
#include <cstdlib>
#include <iostream>
#include <stdexcept>
struct Point { int x=0,y=0; bool infinity=true; };
bool operator==(Point p,Point q) { return p.infinity==q.infinity && (p.infinity||(p.x==q.x&&p.y==q.y)); }
int mod(int x,int p) { return (x%p+p)%p; }
void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
int inv(int x,int p) { for(int i=1;i<p;++i) if(mod(x*i,p)==1)return i; throw std::runtime_error("noninvertible denominator"); }
Point add(Point p,Point q,int prime,int a) {
 if(p.infinity)return q;
 if(q.infinity)return p;
 if(p.x==q.x&&mod(p.y+q.y,prime)==0)return {};
 int slope=p==q?mod((3*p.x*p.x+a)*inv(2*p.y,prime),prime):mod((q.y-p.y)*inv(q.x-p.x,prime),prime);
 int x=mod(slope*slope-p.x-q.x,prime);
 return {x,mod(slope*(p.x-x)-p.y,prime),false};
}
Point mul(int k,Point p,int prime,int a) {
 if(k<0){p.y=mod(-p.y,prime);k=-k;}
 Point q;
 while(k){if(k&1)q=add(q,p,prime,a);p=add(p,p,prime,a);k>>=1;}
 return q;
}
Point phi(Point p) { if(!p.infinity)p.x=mod(6*p.x,43);return p; }
int count(int p,int a,int b) {
 int n=1;
 for(int x=0;x<p;++x)for(int y=0;y<p;++y)if(mod(y*y-x*x*x-a*x-b,p)==0)++n;
 return n;
}
int main() {
 try {
  require(count(43,0,7)==31,"positive point count");
  Point p{2,12,false};require(mul(31,p,43,0).infinity,"generator order");
  require(phi(p)==Point{12,12,false},"map coordinates");
  require(phi(p)==mul(5,p,43,0),"subgroup eigenvalue");
  int points=0;
  for(int k=0;k<31;++k){
   Point q=mul(k,p,43,0);require(add(add(q,phi(q),43,0),phi(phi(q)),43,0).infinity,"map polynomial");
   require(phi(q)==mul(5,q,43,0),"all-point eigenvalue");++points;
  }
  require(mul(14,p,43,0)==add(mul(-1,p,43,0),mul(3,phi(p),43,0),43,0),"GLV identity");
  require(-5*6-1==-31,"positive lattice determinant");
  require(mod(15+5*(-3),31)==0&&14-15==-1,"residual relation");
  require(count(167,25,36)==161,"negative point count");
  require(7*7-4*167==-619&&167-3*7+9==155,"negative trace and norm");
  require(mod(21*21-21+155,23)==0,"negative polynomial root");
  require(2*(-9)-5==-23,"negative lattice determinant");
  int subgroup_points=0;
  for(int x=0;x<167;++x)for(int y=0;y<167;++y)if(mod(y*y-x*x*x-25*x-36,167)==0){
   Point q=mul(7,{x,y,false},167,25);
   require(mul(23,q,167,25).infinity,"negative subgroup");
   require(mul(-2,q,167,25)==mul(21,q,167,25),"negative subgroup action");++subgroup_points;
  }
  std::cout<<"PASS positive_order=31 verified_group_points="<<points<<" negative_order=161 checked_cofactor_images="<<subgroup_points<<" performance_claim=none\n";
 } catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return EXIT_FAILURE;}
}
