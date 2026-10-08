#define main correctness_main
#include "probe.cu"
#undef main
#include <cmath>
#include <cstring>

template<bool Hybrid,int Chains>
__device__ __forceinline__ void capacityBody(const DenseFrag *dense,const HybridFrag *hybrid,
                                            uint32_t *out,unsigned rounds) {
    const unsigned tid=blockIdx.x*blockDim.x+threadIdx.x;
    DenseFrag d;HybridFrag h;
    if constexpr(Hybrid) h=hybrid[tid]; else d=dense[tid];
    uint32_t c0[Chains],c1[Chains],c2[Chains],c3[Chains];
#pragma unroll
    for(int chain=0;chain<Chains;++chain) {
        const uint32_t seed=tid*(4*Chains)+4*chain+1;
        c0[chain]=seed;c1[chain]=seed+1;c2[chain]=seed+2;c3[chain]=seed+3;
    }
#pragma unroll 1
    for(unsigned round=0;round<rounds;++round) {
#pragma unroll
        for(int chain=0;chain<Chains;++chain) {
            if constexpr(Hybrid) {
                denseMma(c0[chain],c1[chain],c2[chain],c3[chain],h.da,h.db);
                sparseMma(c0[chain],c1[chain],c2[chain],c3[chain],h.sa,h.sb,h.e);
            } else {
#pragma unroll
                for(int tile=0;tile<3;++tile)
                    denseMma(c0[chain],c1[chain],c2[chain],c3[chain],d.a[tile],d.b[tile]);
            }
        }
    }
#pragma unroll
    for(int chain=0;chain<Chains;++chain) {
        const size_t at=(size_t(tid)*Chains+chain)*4;
        out[at]=c0[chain];out[at+1]=c1[chain];out[at+2]=c2[chain];out[at+3]=c3[chain];
    }
}

extern "C" __global__ __launch_bounds__(128) void dense4(const DenseFrag *d,const HybridFrag *h,uint32_t *o,unsigned r) {capacityBody<false,4>(d,h,o,r);}
extern "C" __global__ __launch_bounds__(128) void hybrid4(const DenseFrag *d,const HybridFrag *h,uint32_t *o,unsigned r) {capacityBody<true,4>(d,h,o,r);}
extern "C" __global__ __launch_bounds__(128) void dense8(const DenseFrag *d,const HybridFrag *h,uint32_t *o,unsigned r) {capacityBody<false,8>(d,h,o,r);}
extern "C" __global__ __launch_bounds__(128) void hybrid8(const DenseFrag *d,const HybridFrag *h,uint32_t *o,unsigned r) {capacityBody<true,8>(d,h,o,r);}

static int benchmark() {
    cudaDeviceProp prop{};CUDA(cudaGetDeviceProperties(&prop,0));
    if(prop.major!=12 || prop.minor!=0 || prop.multiProcessorCount!=188) {
        std::fprintf(stderr,"Unexpected GPU geometry\n");return 1;
    }
    constexpr unsigned Threads=128,BlocksPerSm=8,Rounds=1024;
    const unsigned blocks=unsigned(prop.multiProcessorCount)*BlocksPerSm;
    const unsigned workers=blocks*Threads,warps=workers/32;
    if(workers>1u<<20 || workers!=192512 || warps!=6016) return 1;
    constexpr uint64_t MaxCoefficient=64ull*129*129;
    constexpr uint64_t OverflowBound=(1ull<<20)*32+Rounds*MaxCoefficient;
    static_assert(OverflowBound==1124139008ull && OverflowBound<0x80000000ull);
    std::vector<Input> inputs(warps);
    std::mt19937 rng(131025);
    for(auto &input:inputs) for(int j=0;j<4;++j) {input.a[j]=rng();input.b[j]=rng();}
    inputs[0]={{{~0u,~0u,~0u,~0u}},{{~0u,~0u,~0u,~0u}}};
    inputs[1]={};
    std::vector<DenseFrag> hostDense(workers);std::vector<HybridFrag> hostHybrid(workers);
    std::vector<std::array<uint32_t,128>> expected(warps);
    for(unsigned warp=0;warp<warps;++warp) {
        pack(inputs[warp],hostDense.data()+32*warp,hostHybrid.data()+32*warp);
        expected[warp]=reference(inputs[warp]);
        for(auto c:expected[warp]) if(c>MaxCoefficient) return 1;
    }
    DenseFrag *dd=nullptr;HybridFrag *dh=nullptr;uint32_t *output=nullptr;
    CUDA(cudaMalloc(&dd,hostDense.size()*sizeof(DenseFrag)));
    CUDA(cudaMalloc(&dh,hostHybrid.size()*sizeof(HybridFrag)));
    CUDA(cudaMalloc(&output,size_t(workers)*8*4*sizeof(uint32_t)));
    CUDA(cudaMemcpy(dd,hostDense.data(),hostDense.size()*sizeof(DenseFrag),cudaMemcpyHostToDevice));
    CUDA(cudaMemcpy(dh,hostHybrid.data(),hostHybrid.size()*sizeof(HybridFrag),cudaMemcpyHostToDevice));
    cudaEvent_t start,stop;CUDA(cudaEventCreate(&start));CUDA(cudaEventCreate(&stop));
    std::vector<uint32_t> hostOutput(size_t(workers)*8*4);
    const char *names[] = {"dense4","hybrid4","dense8","hybrid8"};
    for(int mode=0;mode<4;++mode) {
        cudaFuncAttributes a{};int active=0;
        if(mode==0) {CUDA(cudaFuncGetAttributes(&a,dense4));CUDA(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&active,dense4,Threads,0));}
        if(mode==1) {CUDA(cudaFuncGetAttributes(&a,hybrid4));CUDA(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&active,hybrid4,Threads,0));}
        if(mode==2) {CUDA(cudaFuncGetAttributes(&a,dense8));CUDA(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&active,dense8,Threads,0));}
        if(mode==3) {CUDA(cudaFuncGetAttributes(&a,hybrid8));CUDA(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&active,hybrid8,Threads,0));}
        if(a.localSizeBytes || a.sharedSizeBytes) {std::fprintf(stderr,"Unexpected local/shared memory\n");return 1;}
        std::printf("{\"kind\":\"resources\",\"mode\":\"%s\",\"registers\":%d,\"localBytes\":%zu,\"sharedBytes\":%zu,\"occupancyApiBlocksPerSm\":%d}\n",names[mode],a.numRegs,a.localSizeBytes,a.sharedSizeBytes,active);
    }
    auto sample=[&](int mode,const char *phase,int repeat)->int {
        const unsigned chains=mode<2 ? 4 : 8;
        const bool sparse=mode&1;
        CUDA(cudaEventRecord(start));
        if(mode==0) dense4<<<blocks,Threads>>>(dd,dh,output,Rounds);
        if(mode==1) hybrid4<<<blocks,Threads>>>(dd,dh,output,Rounds);
        if(mode==2) dense8<<<blocks,Threads>>>(dd,dh,output,Rounds);
        if(mode==3) hybrid8<<<blocks,Threads>>>(dd,dh,output,Rounds);
        CUDA(cudaGetLastError());CUDA(cudaEventRecord(stop));CUDA(cudaEventSynchronize(stop));
        float ms=0;CUDA(cudaEventElapsedTime(&ms,start,stop));
        if(!std::isfinite(ms) || ms<0.05f) {std::fprintf(stderr,"Invalid or too-short event interval\n");return 1;}
        const size_t words=size_t(workers)*chains*4;
        CUDA(cudaMemcpy(hostOutput.data(),output,words*sizeof(uint32_t),cudaMemcpyDeviceToHost));
        uint32_t observedMax=0;
        for(unsigned tid=0;tid<workers;++tid) {
            const unsigned lane=tid&31,warp=tid/32;
            for(unsigned chain=0;chain<chains;++chain) for(unsigned q=0;q<4;++q) {
                const unsigned row=(lane>>2)+8*(q>>1),col=2*(lane&3)+(q&1);
                const uint64_t seed=uint64_t(tid)*(4*chains)+4*chain+q+1;
                const uint64_t want=seed+uint64_t(Rounds)*expected[warp][16*col+row];
                const auto actual=hostOutput[(size_t(tid)*chains+chain)*4+q];
                if(want>=0x80000000ull || actual!=want) {
                    std::printf("MISMATCH capacity mode=%s tid=%u chain=%u q=%u want=%llu got=%u\n",names[mode],tid,chain,q,static_cast<unsigned long long>(want),actual);return 1;
                }
                if(actual>observedMax) observedMax=actual;
            }
        }
        const uint64_t cores=uint64_t(warps)*Rounds*chains;
        const uint64_t denseInstructions=cores*(sparse?1:3),sparseInstructions=sparse?cores:0;
        std::printf("{\"kind\":\"sample\",\"phase\":\"%s\",\"repeat\":%d,\"mode\":\"%s\",\"workers\":%u,\"warps\":%u,\"chains\":%u,\"rounds\":%u,\"milliseconds\":%.9f,\"rawCoreEquivalents\":%llu,\"denseMatrixInstructions\":%llu,\"sparseMatrixInstructions\":%llu,\"outputsChecked\":%zu,\"maximumOutput\":%u,\"valid\":true}\n",
                    phase,repeat,names[mode],workers,warps,chains,Rounds,double(ms),
                    static_cast<unsigned long long>(cores),static_cast<unsigned long long>(denseInstructions),
                    static_cast<unsigned long long>(sparseInstructions),words,observedMax);
        std::fflush(stdout);
        return 0;
    };
    for(int mode=0;mode<4;++mode) for(int repeat=0;repeat<2;++repeat) if(sample(mode,"warmup",repeat)) return 1;
    for(int pair=0;pair<2;++pair) for(int repeat=0;repeat<5;++repeat) {
        const int a=2*pair+(repeat&1),b=2*pair+1-(repeat&1);
        if(sample(a,"measurement",repeat) || sample(b,"measurement",repeat)) return 1;
    }
    CUDA(cudaEventDestroy(start));CUDA(cudaEventDestroy(stop));
    CUDA(cudaFree(dd));CUDA(cudaFree(dh));CUDA(cudaFree(output));
    std::printf("{\"kind\":\"capacitySummary\",\"valid\":true,\"warmups\":8,\"measurements\":20,\"workers\":%u,\"rounds\":%u,\"scope\":\"raw matrix capacity including device fragment load and output store; excludes host packing, field reduction, inversion and walk\"}\n",workers,Rounds);
    return 0;
}

int main(int argc,char **argv) {
    if(argc==1 || (argc==2 && std::strcmp(argv[1],"--check")==0)) return correctness_main();
    if(argc==2 && std::strcmp(argv[1],"--bench")==0) {
        const int result=correctness_main();
        if(result) return result;
        return benchmark();
    }
    std::fprintf(stderr,"usage: probe [--check|--bench]\n");return 2;
}
