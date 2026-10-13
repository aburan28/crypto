#pragma once
#define main archived_checker_entry
#include "check.cpp"
#undef main
#include <array>
#include <cerrno>
#include <filesystem>
#include <fcntl.h>
#include <unistd.h>

namespace files = std::filesystem;
inline std::string qstr(const std::string&s){std::string z="\"";for(char c:s){if(c=='"'||c=='\\')z+='\\';if(c=='\n')z+="\\n";else if(c=='\r')z+="\\r";else if(c=='\t')z+="\\t";else z+=c;}return z+'"';}
inline std::string bytes(const files::path&p){std::ifstream f(p,std::ios::binary);require(bool(f),"Cannot read input file");return std::string(std::istreambuf_iterator<char>(f),{});}
inline void create_new(const files::path&p,const std::string&s){int fd=::open(p.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);if(fd<0)throw std::runtime_error("Exclusive create failed: "+p.string());size_t i=0;while(i<s.size()){ssize_t n=::write(fd,s.data()+i,s.size()-i);if(n<0){if(errno==EINTR)continue;::close(fd);throw std::runtime_error("Write failed");}i+=size_t(n);}require(::fsync(fd)==0,"fsync failed");require(::close(fd)==0,"close failed");}
// SHA-256 is used only to bind immutable files and source/executable bytes.
inline uint32_t rr32(uint32_t x,unsigned k){return (x>>k)|(x<<(32-k));}
inline std::string sha256(std::string input){
 static const uint32_t k[64]={0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
 uint64_t bits_length=uint64_t(input.size())*8;input+=char(0x80);while(input.size()%64!=56)input+=char(0);for(int i=7;i>=0;--i)input+=char(bits_length>>(i*8));
 uint32_t h[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
 for(size_t at=0;at<input.size();at+=64){uint32_t w[64];for(int i=0;i<16;++i){w[i]=0;for(int j=0;j<4;++j)w[i]=(w[i]<<8)|uint8_t(input[at+i*4+j]);}for(int i=16;i<64;++i){uint32_t a=w[i-15],b=w[i-2];w[i]=w[i-16]+(rr32(a,7)^rr32(a,18)^(a>>3))+w[i-7]+(rr32(b,17)^rr32(b,19)^(b>>10));}uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];for(int i=0;i<64;++i){uint32_t t1=hh+(rr32(e,6)^rr32(e,11)^rr32(e,25))+((e&f)^((~e)&g))+k[i]+w[i];uint32_t t2=(rr32(a,2)^rr32(a,13)^rr32(a,22))+((a&b)^(a&c)^(b&c));hh=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;}h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=hh;}
 std::ostringstream out;out<<std::hex<<std::setfill('0');for(auto x:h)out<<std::setw(8)<<x;return out.str();
}
inline std::string digest(const files::path&p){return sha256(bytes(p));}
inline U orbit_key(P point){require(!point.inf,"No orbit key for infinity");U key=point.x;U x=point.x;for(int k=1;k<83;++k){x=fs(x);key=std::min(key,x);}return key;}
inline std::string point_json(P p){require(!p.inf,"Affine point required");return "["+qstr(dec(p.x))+","+qstr(dec(p.y))+"]";}
struct Instance {uint64_t r=8569786107849059ULL,lambda=2998004823950457ULL;P g,prior;std::vector<P> base;std::set<U> excluded_orbits;};
inline Instance verify_instance(const files::path&archive){
 auto fixture=read_json((archive/"frozen/fixture.json").string());check_modulus(fixture);require(num(fixture.at("n"))==83&&num(fixture.at("a"))==1,"Wrong curve");
 U x=2;for(int k=0;k<83;++k)x=fs(x);require(x==2&&pgcd(fs(2)^2,MOD)==1,"Rabin irreducibility failure");Instance out;require(prime(out.r)&&num(fixture.at("subgroup_order"))==out.r,"Subgroup prime mismatch");__int128 t0=2,t1=1;for(int n=2;n<=83;++n){auto t=t1-2*t0;t0=t1;t1=t;}require(t1==6151469093347LL&&TOP+1-U(t1)==num(fixture.at("curve_order"))&&num(fixture.at("cofactor"))*out.r==num(fixture.at("curve_order")),"Curve order/cofactor failure");
 out.g=point(fixture.at("generator"));out.prior=point(fixture.at("public_target"));require(out.g.x==decimal("6614244651398184860292015")&&out.g.y==decimal("8868260410029249391539529"),"Exact generator mismatch");require(on(out.g)&&on(out.prior)&&mul(out.r,out.g).inf&&mul(out.r,out.prior).inf,"Generator/prior subgroup failure");require(mul(out.lambda,out.g)==tau(out.g)&&out.lambda!=1&&power(out.lambda,83,out.r)==1&&(mm(out.lambda,out.lambda,out.r)+out.r-out.lambda+2)%out.r==0,"Frobenius eigenvalue failure");
 auto header=read_json((archive/"base_n83_K600.jsonl").string());check_modulus(header);require(num(header.at("n"))==83&&num(header.at("a"))==1&&num(header.at("subgroup_order"))==out.r,"Base curve mismatch");const auto&raw=header.at("factor_base_point_coordinates").array;const auto&labels=header.at("factor_base_point_labels").array;const auto&reps=header.at("factor_base_representatives").array;require(raw.size()==99600&&labels.size()==99600&&reps.size()==600,"Base size mismatch");std::set<std::pair<U,U>> distinct;out.base.reserve(raw.size());
 for(size_t c=0;c<600;++c){P rep=point(reps[c]);require(on(rep)&&mul(out.r,rep).inf,"Base representative subgroup failure");require(out.excluded_orbits.insert(orbit_key(rep)).second,"Duplicate representative orbit");P current=rep;uint64_t coefficient=1;for(int k=0;k<83;++k){for(int sign=0;sign<2;++sign){size_t at=166*c+2*k+sign;P stored=point(raw[at]);require(on(stored)&&stored==(sign?neg(current):current),"Signed orbit point mismatch");require(num(labels[at].at(0))==c&&num(labels[at].at(1))==(sign?out.r-coefficient:coefficient),"Orbit label mismatch");require(distinct.insert({stored.x,stored.y}).second,"Duplicate base point");out.base.push_back(stored);}current=tau(current);coefficient=mm(coefficient,out.lambda,out.r);}require(current==rep&&coefficient==1,"Base orbit closure failure");}
 out.excluded_orbits.insert(orbit_key(out.prior));return out;
}
struct SplitMix64 {uint64_t state;uint64_t draws=0;uint64_t next(){state+=0x9e3779b97f4a7c15ULL;uint64_t z=state;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;++draws;return z^(z>>31);}};
