#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <vector>
#include "ecc2k130/include/curveparams.h"
#include "ecc2k130/include/walk.h"
using C=CfgF131;using R=Ref<C>;using P=R::Point;
static uint32_t r32(const std::vector<unsigned char>&b,size_t o){uint32_t x;std::memcpy(&x,b.data()+o,4);return x;}static uint64_t r64(const std::vector<unsigned char>&b,size_t o){uint64_t x;std::memcpy(&x,b.data()+o,8);return x;}
int main(int argc,char**argv){if(argc!=2)return 2;std::ifstream f(argv[1],std::ios::binary);std::vector<unsigned char>b((std::istreambuf_iterator<char>(f)),{});if(b.size()<40||std::memcmp(b.data(),"ECC2K130",8))return 3;unsigned version=r32(b,8),degree=r32(b,12),threads=r32(b,16),batch=r32(b,20),lanes=r32(b,24),runid=r32(b,28);uint64_t steps=r64(b,32);size_t n=size_t(threads)*batch;if(version!=6||degree!=131||lanes!=1||b.size()!=40+n*60)return 4;
 U192 ell=u192_from_dec(eccF131::ELL_DEC),s=u192_from_dec(eccF131::S_DEC),sp[256];sp[0]=u192_from(1);for(int i=1;i<256;i++)sp[i]=mod_mul(sp[i-1],s,ell);P B=R::make(R::fromLimbs(eccF131::PX),R::fromLimbs(eccF131::PY)),Q=R::make(R::fromLimbs(eccF131::QX),R::fromLimbs(eccF131::QY));size_t bad=0;
 for(unsigned slot=0;slot<batch;slot++)for(unsigned tid=0;tid<threads;tid++){size_t id=size_t(slot)*threads+tid;P p=R::startPoint(eccSeedFor(runid,id),B,Q,nullptr,ell,sp);for(uint64_t k=0;k<steps;k++)p=((R::weight(p.x)%72)==14)?R::addPt(p,R::frob(p,3)):R::addPt(p,R::frob(p,1));uint64_t gx[3]{};for(int w=0;w<5;w++){uint32_t x=r32(b,40+((size_t(slot)*5+w)*threads+tid)*4);gx[w/2]|=uint64_t(x)<<(32*(w&1));}if(std::memcmp(gx,p.x.v,sizeof gx))bad++;}
 std::printf("bridge1-sparse-bridge3-mod72 full-x oracle: points=%zu steps=%llu mismatches=%zu\n",n,(unsigned long long)steps,bad);return bad?1:0;}
