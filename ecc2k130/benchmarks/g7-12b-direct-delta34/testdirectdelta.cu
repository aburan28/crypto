#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#define ECC_PACKED_CLMAD 1
#define ECC_PACKED_DIRECT_REDUCE 1
#define ECC_PACKED_TAIL_LAYOUT 1
#define ECC_PACKED_FAST_CONVERT 1
#include "../../build/g7-direct-delta34-src/include/packed131.h"
using eccPacked131::P131;
struct Pair { P131 x,y; };
static void ck(cudaError_t e){if(e!=cudaSuccess){std::fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e));std::exit(1);}}
static __host__ __device__ bool same(P131 a,P131 b){for(int i=0;i<5;i++)if(a.v[i]!=b.v[i])return false;return true;}
__device__ P131 direct(P131 a,int j){
 P131 s=eccPacked131::squarePolynomial131(a);
 s=eccPacked131::squarePolynomial131(s);
 s=eccPacked131::squarePolynomial131(s);
 if(j==4)s=eccPacked131::squarePolynomial131(s);
 return eccPacked131::add131(a,s);
}
__device__ P131 independent(P131 a,int j){
 P131 n=eccPacked131::fromPolynomial131(a), s=n;
 for(int k=0;k<j;k++)s=eccPacked131::sqr131(s);
 return eccPacked131::toPolynomial131(eccPacked131::add131(n,s));
}
__global__ void probe(const Pair*in,Pair*out,const int*powers,unsigned*bad,int n){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
 Pair got{direct(in[i].x,powers[i]),direct(in[i].y,powers[i])};
 Pair ref{independent(in[i].x,powers[i]),independent(in[i].y,powers[i])};
 if(!same(got.x,ref.x)||!same(got.y,ref.y))atomicAdd(bad,1u);
 out[i]=got;
}
int main(){
 std::vector<Pair> in;std::vector<int> powers;unsigned state=0x34d131u;
 auto rnd=[&](){state^=state<<13;state^=state>>17;state^=state<<5;return state;};
 for(int j=3;j<=4;j++){
  for(int bit=0;bit<131;bit++){P131 a{},b{};a.v[bit/32]=1u<<(bit%32);b.v[(bit*37)%131/32]=1u<<((bit*37)%131%32);in.push_back({a,b});powers.push_back(j);}
  const Pair edges[]={{{0,0,0,0,0},{1,0,0,0,0}},{{~0u,~0u,~0u,~0u,7},{0,0,0,0,0}},{{0xaaaaaaaau,0x55555555u,0x12345678u,0x87654321u,5},{1,2,4,8,3}}};
  for(auto e:edges){in.push_back(e);powers.push_back(j);}
  for(int row=0;row<8192;row++){Pair e{};for(int w=0;w<5;w++){e.x.v[w]=rnd();e.y.v[w]=rnd();}e.x.v[4]&=7;e.y.v[4]&=7;in.push_back(e);powers.push_back(j);}
 }
 Pair *di,*dout;int*dp;unsigned*dbad;ck(cudaMalloc(&di,in.size()*sizeof(Pair)));ck(cudaMalloc(&dout,in.size()*sizeof(Pair)));ck(cudaMalloc(&dp,powers.size()*sizeof(int)));ck(cudaMalloc(&dbad,sizeof(unsigned)));
 ck(cudaMemcpy(di,in.data(),in.size()*sizeof(Pair),cudaMemcpyHostToDevice));ck(cudaMemcpy(dp,powers.data(),powers.size()*sizeof(int),cudaMemcpyHostToDevice));ck(cudaMemset(dout,0xa5,in.size()*sizeof(Pair)));ck(cudaMemset(dbad,0,sizeof(unsigned)));
 probe<<<(in.size()+255)/256,256>>>(di,dout,dp,dbad,int(in.size()));ck(cudaGetLastError());ck(cudaDeviceSynchronize());
 unsigned bad;ck(cudaMemcpy(&bad,dbad,sizeof bad,cudaMemcpyDeviceToHost));if(bad){std::fprintf(stderr,"FAIL: %u mismatches\n",bad);return 1;}
 std::vector<Pair> got(in.size());ck(cudaMemcpy(got.data(),dout,got.size()*sizeof(Pair),cudaMemcpyDeviceToHost));
 std::printf("PASS: %zu polynomial delta pairs, direct squares equal independent normal-basis route\n",got.size());
}
