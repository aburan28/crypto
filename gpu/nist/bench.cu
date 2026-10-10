#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "nist_solinas.cuh"
#include "nist_rho.cuh"
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
/* ---- Pollard-rho walk throughput --------------------------------------
 * The hot kernel: every resident thread advances its W walks one folded,
 * negation-mapped r-adding step per outer iteration, batching the one field
 * inversion across its W walks.  Same rho_batch_step the CPU test verifies,
 * so the device result is correct by construction; here it is only timed and
 * its distinguished-point output counted.  W is capped at RHO_BENCH_W so the
 * per-thread den/scratch live in registers/local, not dynamic memory. */
#define RHO_BENCH_W 8
template<class C> __global__ void k_rho(rho_ctx<C> c,int iters){
 typedef typename C::Field F;
 uint32_t tid=blockIdx.x*blockDim.x+threadIdx.x; if(tid>=c.nthreads) return;
 typename F::elt den[RHO_BENCH_W],scratch[RHO_BENCH_W];
 for(int it=0;it<iters;it++) rho_batch_step<C>(c,tid,den,scratch);
}
template<class C> __global__ void k_rho_seed(rho_ctx<C> c){
 uint32_t tid=blockIdx.x*blockDim.x+threadIdx.x; uint32_t nw=rho_nwalks(c);
 for(uint32_t w=0;w<c.walks_per_thread;w++){
  uint32_t idx=tid+w*c.nthreads; if(idx>=nw) return;
  c.restarts[idx]=0; c.steps[idx]=0;
  rho_state<C> st; rho_seed_walk<C>(c,idx,st); rho_store<C>(c,idx,st);
 }
}

template<class F> static nap<F> parse_pt(const char*xs,const char*ys){
 auto parse=[&](const char*h)->typename F::elt{
  const char*p=h; if(p[0]=='0'&&(p[1]=='x'||p[1]=='X'))p+=2; int len=0; while(p[len])len++;
  uint32_t l[F::N]={0}; int nib=0;
  for(int i=len-1;i>=0;i--){int ch=p[i];int d=(ch>='0'&&ch<='9')?ch-'0':(ch|32)-'a'+10;l[nib/8]|=(uint32_t)d<<(4*(nib%8));nib++;}
  return F::from_limbs(l);
 };
 nap<F> P; P.x=parse(xs); P.y=parse(ys); P.inf=0; return P;
}

template<class C> static void rho_bench(const char*name,const nap<typename C::Field>&G,
                                        const nap<typename C::Field>&Q,int threads,int W,int iters){
 typedef typename C::Field F;
 if(W>RHO_BENCH_W)W=RHO_BENCH_W;
 int blocks=1; // tune below
 cudaDeviceProp p{}; CU(cudaGetDeviceProperties(&p,0));
 blocks=p.multiProcessorCount*16;
 uint32_t nthreads=(uint32_t)blocks*threads, nw=nthreads*W;
 rho_params prm{}; prm.r_bits=8; prm.dp_mask=(1u<<24)-1u; prm.fold=RHO_FOLD_NEG; prm.max_steps=1u<<28; prm.table_seed=7; prm.order_bits=F::N*32;
 uint32_t R=1u<<prm.r_bits;
 /* jump table on the host, then copy */
 std::vector<nap<F>> tbl(R);
 for(uint32_t j=0;j<R;j++){uint32_t cc[F::N],dd[F::N];rho_table_seed<F::N>(prm.table_seed,j,cc,dd,prm.order_bits);tbl[j]=C::double_scalar_mul(G,cc,Q,dd,prm.order_bits);}
 rho_ctx<C> c{}; c.nthreads=nthreads; c.walks_per_thread=(uint32_t)W; c.prm=prm; c.G=G; c.Q=Q; c.dp_cap=1u<<20;
 nap<F>*dtbl; CU(cudaMalloc(&dtbl,R*sizeof(nap<F>))); CU(cudaMemcpy(dtbl,tbl.data(),R*sizeof(nap<F>),cudaMemcpyHostToDevice)); c.table=dtbl;
 CU(cudaMalloc(&c.X,sizeof(uint32_t)*F::N*nw)); CU(cudaMalloc(&c.Y,sizeof(uint32_t)*F::N*nw));
 CU(cudaMalloc(&c.H,sizeof(uint32_t)*RHO_CYCLE_DEPTH*nw)); CU(cudaMalloc(&c.esc,sizeof(uint32_t)*nw));
 CU(cudaMalloc(&c.steps,sizeof(uint32_t)*nw)); CU(cudaMalloc(&c.restarts,sizeof(uint32_t)*nw));
 CU(cudaMalloc(&c.dp_out,sizeof(rho_dp<F>)*c.dp_cap)); CU(cudaMalloc(&c.dp_count,sizeof(uint32_t)));
 CU(cudaMemset(c.esc,0,sizeof(uint32_t)*nw)); CU(cudaMemset(c.dp_count,0,sizeof(uint32_t)));
 k_rho_seed<C><<<blocks,threads>>>(c); CU(cudaDeviceSynchronize());
 k_rho<C><<<blocks,threads>>>(c,4); CU(cudaDeviceSynchronize());   /* warm-up */
 Timer t; t.runStart(); k_rho<C><<<blocks,threads>>>(c,iters); double s=t.stop();
 uint32_t dps=0; CU(cudaMemcpy(&dps,c.dp_count,sizeof(uint32_t),cudaMemcpyDeviceToHost));
 double steps=(double)nw*iters;
 printf("%s rho: walks=%u threads/blk=%d blocks=%d  %.3f Gstep/s  (%u DPs, dp_bits=24)\n",
        name,nw,threads,blocks,steps/s/1e9,dps);
 cudaFree(dtbl);cudaFree(c.X);cudaFree(c.Y);cudaFree(c.H);cudaFree(c.esc);cudaFree(c.steps);cudaFree(c.restarts);cudaFree(c.dp_out);cudaFree(c.dp_count);
}

int main(int argc,char**argv){
 int threads=argc>1?atoi(argv[1]):128,n=1<<20,iters=128;
 cudaDeviceProp p{};CU(cudaGetDeviceProperties(&p,0));
 printf("GPU %s cc %d.%d SMs=%d regs/SM=%d NIST_PTX=%d\n",p.name,p.major,p.minor,p.multiProcessorCount,p.regsPerMultiprocessor,(int)NIST_PTX);
 bench<Fp256>("P-256 Montgomery",threads,n,iters);
 bench<Fp256Sol>("P-256 Solinas",threads,n,iters);
 bench<Fp384>("P-384 Montgomery",threads,n,iters/2);
 bench<Fp384Sol>("P-384 Solinas",threads,n,iters/2);
 /* rho walk throughput on P-256 and P-384 (Q = [k]G for an arbitrary k). */
 nap<Fp256> G256=parse_pt<Fp256>(
  "0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296",
  "0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
 uint32_t k7[8]={7,0,0,0,0,0,0,0}; nap<Fp256> Q256=CurveP256::scalar_mul(G256,k7,8);
 rho_bench<CurveP256>("P-256",G256,Q256,threads,RHO_BENCH_W,1<<16);
 nap<Fp384> G384=parse_pt<Fp384>(
  "0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7",
  "0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f");
 uint32_t k7b[12]={7,0,0,0,0,0,0,0,0,0,0,0}; nap<Fp384> Q384=CurveP384::scalar_mul(G384,k7b,8);
 rho_bench<CurveP384>("P-384",G384,Q384,threads,RHO_BENCH_W,1<<15);
 return 0;
}
