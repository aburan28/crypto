#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "../../build/g7-two-jump-src/include/packed131.h"
namespace eccPacked131 {
#include "../../build/g7-two-jump-src/include/packedtwosigma131.h"
}
using eccPacked131::P131;
struct Pair { P131 a,b; };
static void ck(cudaError_t e){if(e!=cudaSuccess){std::fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e));std::exit(1);}}
static bool same(P131 a,P131 b){return !std::memcmp(a.v,b.v,sizeof a.v);}
static P131 reference(P131 a,int exponent){
 unsigned factor=1;for(int k=0;k<exponent;++k)factor=(2*factor)%263;P131 r{};
 for(int bit=0;bit<131;++bit)if((a.v[bit/32]>>(bit%32))&1u){
  unsigned target=((bit+1)*factor)%263;if(target>131)target=263-target;--target;
  r.v[target/32]^=1u<<(target%32);
 }return r;
}
__global__ void probe(const Pair*in,Pair*out,const int*powers,int n){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
 auto r=eccPacked131::sigmaWalkPairImmediate34(in[i].a,in[i].b,powers[i]-3);
 out[i]={r.first,r.second};
}
int main(){
 std::vector<Pair> in,want;std::vector<int> powers;unsigned state=0x234131u;
 auto rnd=[&](){state^=state<<13;state^=state>>17;state^=state<<5;return state;};
 for(int power=3;power<=4;++power){
  for(int bit=0;bit<131;++bit){P131 x{};x.v[bit/32]=1u<<(bit%32);in.push_back({x,x});powers.push_back(power);}
  for(int row=0;row<8192;++row){P131 x{},y{};for(int w=0;w<5;++w){x.v[w]=rnd();y.v[w]=rnd();}x.v[4]&=7;y.v[4]&=7;in.push_back({x,y});powers.push_back(power);}
 }
 for(size_t i=0;i<in.size();++i)want.push_back({reference(in[i].a,powers[i]),reference(in[i].b,powers[i])});
 Pair *di,*dout;int*dp;ck(cudaMalloc(&di,in.size()*sizeof(Pair)));ck(cudaMalloc(&dout,in.size()*sizeof(Pair)));ck(cudaMalloc(&dp,powers.size()*sizeof(int)));
 ck(cudaMemcpy(di,in.data(),in.size()*sizeof(Pair),cudaMemcpyHostToDevice));ck(cudaMemcpy(dp,powers.data(),powers.size()*sizeof(int),cudaMemcpyHostToDevice));
 probe<<<(in.size()+255)/256,256>>>(di,dout,dp,int(in.size()));ck(cudaGetLastError());ck(cudaDeviceSynchronize());
 std::vector<Pair> got(in.size());ck(cudaMemcpy(got.data(),dout,got.size()*sizeof(Pair),cudaMemcpyDeviceToHost));
 for(size_t i=0;i<got.size();++i)if(!same(got[i].a,want[i].a)||!same(got[i].b,want[i].b)){std::fprintf(stderr,"FAIL %zu power %d\n",i,powers[i]);return 1;}
 std::printf("PASS: %zu fixed sigma3/sigma4 GPU pairs against independent basis routing\n",got.size());
}
