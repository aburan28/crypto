#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "nist_ptx.cuh"
#define CU(x) do{cudaError_t e=(x);if(e!=cudaSuccess){fprintf(stderr,"%s\n",cudaGetErrorString(e));exit(1);}}while(0)

template<class F,int CHAINS> __global__ void k_mul(typename F::elt*out,const typename F::elt*in,int iters,int n){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
 typename F::elt x[CHAINS];
#pragma unroll
 for(int c=0;c<CHAINS;c++)x[c]=in[(i+c*131)&(n-1)];
 for(int k=0;k<iters;k++){
#pragma unroll
  for(int c=0;c<CHAINS;c++)x[c]=F::mul(x[c],in[(i+17+c*313)&(n-1)]);
 }
 typename F::elt z=x[0];
#pragma unroll
 for(int c=1;c<CHAINS;c++)z=F::add(z,x[c]);
 out[i]=z;
}
template<class F> __global__ void k_sqr(typename F::elt*out,const typename F::elt*in,int iters,int n){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;auto x=in[i];
 for(int k=0;k<iters;k++)x=F::sqr(x);out[i]=x;
}
template<class F> __global__ void k_dbl(njac<F>*out,const njac<F>*in,int iters,int n){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;auto x=in[i];
 for(int k=0;k<iters;k++)x=n_double<F>(x);out[i]=x;
}
template<class F> __global__ void k_madd(njac<F>*out,const njac<F>*in,const naff<F>*q,int iters,int n){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;auto x=in[i];auto a=q[i];
 for(int k=0;k<iters;k++)x=n_madd<F>(x,a);out[i]=x;
}
struct Timer{cudaEvent_t a,b;Timer(){cudaEventCreate(&a);cudaEventCreate(&b);}float runStart(){cudaEventRecord(a);return 0;}double stop(){cudaEventRecord(b);cudaEventSynchronize(b);float ms;cudaEventElapsedTime(&ms,a,b);return ms/1000.;}};

template<class F> typename F::elt seedfe(uint64_t s){
 typename F::elt a{};for(int i=0;i<F::N;i++){s+=0x9e3779b97f4a7c15ull;uint64_t z=s;z=(z^(z>>30))*0xbf58476d1ce4e5b9ull;z=(z^(z>>27))*0x94d049bb133111ebull;z^=z>>31;a.v[i]=(uint32_t)z;}
 typename F::elt p=F::modulus(),d;uint32_t bw=n_sub<F::N>(d.v,a.v,p.v);n_cmov<F::N>(a.v,d.v,bw^1u);return F::from_canonical(a);
}
template<class F> void bench(const char*name,int threads,int n,int iters){
 using E=typename F::elt;std::vector<E> h(n);for(int i=0;i<n;i++)h[i]=seedfe<F>(0x50323536ull+i);
 E *di,*doo;CU(cudaMalloc(&di,n*sizeof(E)));CU(cudaMalloc(&doo,n*sizeof(E)));CU(cudaMemcpy(di,h.data(),n*sizeof(E),cudaMemcpyHostToDevice));
 int blocks=(n+threads-1)/threads;Timer t;
 k_mul<F,1><<<blocks,threads>>>(doo,di,2,n);CU(cudaDeviceSynchronize());
 t.runStart();k_mul<F,1><<<blocks,threads>>>(doo,di,iters,n);double s=t.stop();
 printf("%s mul chain1 threads=%d: %.3f Gmul/s\n",name,threads,(double)n*iters/s/1e9);
 t.runStart();k_mul<F,2><<<blocks,threads>>>(doo,di,iters/2,n);s=t.stop();
 printf("%s mul chain2 threads=%d: %.3f Gmul/s\n",name,threads,(double)n*iters/s/1e9);
 t.runStart();k_sqr<F><<<blocks,threads>>>(doo,di,iters,n);s=t.stop();
 printf("%s sqr threads=%d: %.3f Gsqr/s\n",name,threads,(double)n*iters/s/1e9);
 std::vector<njac<F>> hj(n);std::vector<naff<F>> ha(n);
 for(int i=0;i<n;i++){ha[i]={h[i],h[(i+1)&(n-1)]};hj[i]={ha[i].x,ha[i].y,F::one()};}
 njac<F>*dj,*djo;naff<F>*da;CU(cudaMalloc(&dj,n*sizeof(*dj)));CU(cudaMalloc(&djo,n*sizeof(*djo)));CU(cudaMalloc(&da,n*sizeof(*da)));
 CU(cudaMemcpy(dj,hj.data(),n*sizeof(*dj),cudaMemcpyHostToDevice));CU(cudaMemcpy(da,ha.data(),n*sizeof(*da),cudaMemcpyHostToDevice));
 int pit=iters/16;if(pit<1)pit=1;
 t.runStart();k_dbl<F><<<blocks,threads>>>(djo,dj,pit,n);s=t.stop();
 printf("%s point_double threads=%d: %.3f M/s\n",name,threads,(double)n*pit/s/1e6);
 t.runStart();k_madd<F><<<blocks,threads>>>(djo,dj,da,pit,n);s=t.stop();
 printf("%s mixed_add threads=%d: %.3f M/s\n",name,threads,(double)n*pit/s/1e6);
 cudaFree(di);cudaFree(doo);cudaFree(dj);cudaFree(djo);cudaFree(da);
}
int main(int argc,char**argv){
 int threads=argc>1?atoi(argv[1]):128,n=1<<20,iters=128;
 cudaDeviceProp p{};CU(cudaGetDeviceProperties(&p,0));
 printf("GPU %s cc %d.%d SMs=%d regs/SM=%d NIST_PTX=%d\n",p.name,p.major,p.minor,p.multiProcessorCount,p.regsPerMultiprocessor,(int)NIST_PTX);
 bench<Fp256>("P-256",threads,n,iters);
 bench<Fp384>("P-384",threads,n,iters/2);
 return 0;
}
