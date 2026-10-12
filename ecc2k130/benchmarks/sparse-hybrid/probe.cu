#include <cuda_runtime.h>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#define CUDA(call) do { const cudaError_t error=(call); if(error!=cudaSuccess) { \
    std::fprintf(stderr,"CUDA error: %s\n",cudaGetErrorString(error)); return 2; } } while(0)

struct DenseFrag { uint32_t a[3][4], b[3][2]; };
struct HybridFrag { uint32_t da[4], db[2], sa[4], sb[4], e; };
struct Input { std::array<uint32_t,4> a,b; };
static_assert(sizeof(DenseFrag)*32==2304);
static_assert(sizeof(HybridFrag)*32==1920);

__device__ __forceinline__ void denseMma(uint32_t &c0,uint32_t &c1,uint32_t &c2,uint32_t &c3,
                                        const uint32_t *a,const uint32_t *b) {
    asm volatile("mma.sync.aligned.m16n8k32.row.col.s32.u8.u8.s32 "
                 "{%0,%1,%2,%3},{%4,%5,%6,%7},{%8,%9},{%0,%1,%2,%3};"
                 : "+r"(c0),"+r"(c1),"+r"(c2),"+r"(c3)
                 : "r"(a[0]),"r"(a[1]),"r"(a[2]),"r"(a[3]),"r"(b[0]),"r"(b[1]));
}

__device__ __forceinline__ void sparseMma(uint32_t &c0,uint32_t &c1,uint32_t &c2,uint32_t &c3,
                                         const uint32_t *a,const uint32_t *b,uint32_t e) {
    asm volatile("mma.sp::ordered_metadata.sync.aligned.m16n8k64.row.col.s32.u8.u8.s32 "
                 "{%0,%1,%2,%3},{%4,%5,%6,%7},{%8,%9,%10,%11},{%0,%1,%2,%3},%12,0;"
                 : "+r"(c0),"+r"(c1),"+r"(c2),"+r"(c3)
                 : "r"(a[0]),"r"(a[1]),"r"(a[2]),"r"(a[3]),
                   "r"(b[0]),"r"(b[1]),"r"(b[2]),"r"(b[3]),"r"(e));
}

extern "C" __global__ __launch_bounds__(32) void denseKernel(const DenseFrag *input,uint32_t *out) {
    const size_t id=size_t(blockIdx.x)*32+threadIdx.x;
    const DenseFrag f=input[id];
    uint32_t c0=0,c1=0,c2=0,c3=0;
#pragma unroll
    for(int tile=0;tile<3;++tile) denseMma(c0,c1,c2,c3,f.a[tile],f.b[tile]);
    out[4*id]=c0;out[4*id+1]=c1;out[4*id+2]=c2;out[4*id+3]=c3;
}

extern "C" __global__ __launch_bounds__(32) void hybridKernel(const HybridFrag *input,uint32_t *out) {
    const size_t id=size_t(blockIdx.x)*32+threadIdx.x;
    const HybridFrag f=input[id];
    uint32_t c0=0,c1=0,c2=0,c3=0;
    denseMma(c0,c1,c2,c3,f.da,f.db);
    sparseMma(c0,c1,c2,c3,f.sa,f.sb,f.e);
    out[4*id]=c0;out[4*id+1]=c1;out[4*id+2]=c2;out[4*id+3]=c3;
}

static std::array<uint32_t,64> encode(const std::array<uint32_t,4> &a) {
    std::array<uint32_t,64> out{};
    for(int i=0;i<64;++i) {
        const uint32_t bits=a[i/16]>>(2*(i%16));
        out[i]=(bits&1)+128*((bits>>1)&1);
    }
    return out;
}

static uint32_t aval(const std::array<uint32_t,64>& a,int row,int k) {
    const int i=112+row-k;
    return i>=0 && i<64 ? a[i] : 0;
}
static uint32_t bval(const std::array<uint32_t,64>& b,int k,int col) {
    const int j=k+16*col-112;
    return j>=0 && j<64 ? b[j] : 0;
}

static void pack(const Input &input,DenseFrag *dense,HybridFrag *hybrid) {
    const auto aa=encode(input.a),bb=encode(input.b);
    int columns[64],kept[16][32];uint32_t meta[16][16];
    for(int p=0;p<32;++p) {
        columns[2*p]=p<17 ? 96+p : 49+(p-17);
        columns[2*p+1]=p<17 ? 32+p : 113+(p-17);
    }
    for(int row=0;row<16;++row) {
        int count=0;
        for(int group=0;group<16;++group) {
            int indices[2],n=0;
            for(int j=0;j<4;++j) {
                const int k=columns[4*group+j],ai=112+row-k;
                if(ai>=0 && ai<64) {
                    if(n==2) std::abort();
                    indices[n++]=j;kept[row][count++]=k;
                }
            }
            if(n!=2 || indices[0]>=indices[1]) std::abort();
            meta[row][group]=unsigned(indices[0]+4*indices[1]);
        }
        if(count!=32) std::abort();
    }
    for(int lane=0;lane<32;++lane) {
        const int g=lane>>2,t=lane&3;
        dense[lane]={};hybrid[lane]={};
        for(int q=0;q<4;++q) for(int byte=0;byte<4;++byte) {
            const int row=g+8*(q&1),col=4*t+byte+16*(q>>1);
            for(int tile=0;tile<3;++tile)
                dense[lane].a[tile][q] |= aval(aa,row,32+32*tile+col)<<(8*byte);
            hybrid[lane].da[q] |= aval(aa,row,64+col)<<(8*byte);
            hybrid[lane].sa[q] |= aval(aa,row,kept[row][col])<<(8*byte);
            const int k=16*q+4*t+byte;
            hybrid[lane].sb[q] |= bval(bb,columns[k],g)<<(8*byte);
        }
        for(int q=0;q<2;++q) for(int byte=0;byte<4;++byte) {
            const int col=16*q+4*t+byte;
            for(int tile=0;tile<3;++tile)
                dense[lane].b[tile][q] |= bval(bb,32+32*tile+col,g)<<(8*byte);
            hybrid[lane].db[q] |= bval(bb,64+col,g)<<(8*byte);
        }
        // E ownership differs from A ownership: a lane supplies eight entire
        // ordered metadata nibbles for one row and one half of the K dimension.
        const int erow=g+8*(t&1),firstGroup=8*(t>>1);
        for(int j=0;j<8;++j) hybrid[lane].e |= meta[erow][firstGroup+j]<<(4*j);
    }
}

static std::array<uint32_t,128> reference(const Input &input) {
    const auto a=encode(input.a),b=encode(input.b);
    std::array<uint32_t,128> out{};
    for(int i=0;i<64;++i) for(int j=0;j<64;++j) out[i+j]+=a[i]*b[j];
    return out;
}
static std::array<uint32_t,8> serial(const Input &input) {
    std::array<uint32_t,8> out{};
    for(int i=0;i<128;++i) if((input.a[i/32]>>(i%32))&1)
        for(int j=0;j<128;++j) if((input.b[j/32]>>(j%32))&1)
            out[(i+j)/32]^=uint32_t(1)<<((i+j)%32);
    return out;
}
static std::array<uint32_t,8> reconstruct(const std::array<uint32_t,128>& c,const Input &input) {
    std::array<uint32_t,8> out{};
    auto bit=[&](int i,uint32_t value) { if(i<256 && value) out[i/32]^=uint32_t(1)<<(i%32); };
    for(int i=0;i<128;++i) { bit(2*i,c[i]&1);bit(2*i+1,(c[i]>>7)&1);bit(2*i+2,(c[i]>>14)&1); }
    uint32_t all=~uint32_t(0);
    for(int i=0;i<4;++i) all &= input.a[i]&input.b[i];
    if(all==~uint32_t(0)) out[4]^=1;
    return out;
}

int main() {
    std::vector<Input> cases;
    for(int i=0;i<128;++i) for(int j=0;j<128;++j) {
        Input input{};input.a[i/32]=uint32_t(1)<<(i%32);input.b[j/32]=uint32_t(1)<<(j%32);
        cases.push_back(input);
    }
    const std::array<std::array<uint32_t,4>,8> edges={{{0,0,0,0},{1,0,0,0},
        {~0u,~0u,~0u,~0u},{~1u,~0u,~0u,~0u},{0x55555555u,0x55555555u,0x55555555u,0x55555555u},
        {0xaaaaaaaau,0xaaaaaaaau,0xaaaaaaaau,0xaaaaaaaau},{0,0,0,0x80000000u},{~0u,~0u,~0u,0x7fffffffu}}};
    for(auto a:edges) for(auto b:edges) cases.push_back({a,b});
    std::mt19937 rng(131024);
    for(int i=0;i<64;++i) {
        Input input;for(int j=0;j<4;++j) { input.a[j]=rng();input.b[j]=rng(); }cases.push_back(input);
    }
    const size_t n=cases.size();
    std::vector<DenseFrag> dense(n*32);std::vector<HybridFrag> hybrid(n*32);
    for(size_t i=0;i<n;++i) pack(cases[i],dense.data()+32*i,hybrid.data()+32*i);
    CUDA(cudaSetDevice(0));cudaDeviceProp prop{};CUDA(cudaGetDeviceProperties(&prop,0));
    std::printf("device: %s, sm_%d%d, %d SMs\n",prop.name,prop.major,prop.minor,prop.multiProcessorCount);
    DenseFrag *dd=nullptr;HybridFrag *dh=nullptr;uint32_t *od=nullptr,*oh=nullptr;
    CUDA(cudaMalloc(&dd,dense.size()*sizeof(DenseFrag)));CUDA(cudaMalloc(&dh,hybrid.size()*sizeof(HybridFrag)));
    CUDA(cudaMalloc(&od,n*128*sizeof(uint32_t)));CUDA(cudaMalloc(&oh,n*128*sizeof(uint32_t)));
    CUDA(cudaMemcpy(dd,dense.data(),dense.size()*sizeof(DenseFrag),cudaMemcpyHostToDevice));
    CUDA(cudaMemcpy(dh,hybrid.data(),hybrid.size()*sizeof(HybridFrag),cudaMemcpyHostToDevice));
    denseKernel<<<unsigned(n),32>>>(dd,od);CUDA(cudaGetLastError());CUDA(cudaDeviceSynchronize());
    hybridKernel<<<unsigned(n),32>>>(dh,oh);CUDA(cudaGetLastError());CUDA(cudaDeviceSynchronize());
    std::vector<uint32_t> denseOut(n*128),hybridOut(n*128);
    CUDA(cudaMemcpy(denseOut.data(),od,denseOut.size()*sizeof(uint32_t),cudaMemcpyDeviceToHost));
    CUDA(cudaMemcpy(hybridOut.data(),oh,hybridOut.size()*sizeof(uint32_t),cudaMemcpyDeviceToHost));
    for(size_t i=0;i<n;++i) {
        const auto expect=reference(cases[i]);std::array<uint32_t,128> dc{},hc{};
        for(int lane=0;lane<32;++lane) for(int q=0;q<4;++q) {
            const int row=(lane>>2)+8*(q>>1),col=2*(lane&3)+(q&1),l=16*col+row;
            dc[l]=denseOut[(i*32+lane)*4+q];hc[l]=hybridOut[(i*32+lane)*4+q];
        }
        for(int l=0;l<128;++l) if(dc[l]!=expect[l] || hc[l]!=expect[l]) {
            std::printf("MISMATCH case=%zu coefficient=%d expected=%u dense=%u hybrid=%u\n",i,l,expect[l],dc[l],hc[l]);return 1;
        }
        if(reconstruct(dc,cases[i])!=serial(cases[i]) || reconstruct(hc,cases[i])!=serial(cases[i])) {
            std::printf("MISMATCH raw reconstruction case=%zu\n",i);return 1;
        }
    }
    CUDA(cudaFree(dd));CUDA(cudaFree(dh));CUDA(cudaFree(od));CUDA(cudaFree(oh));
    std::printf("{\"valid\":true,\"basisPairs\":16384,\"edgePairs\":64,\"densePairs\":64,\"cases\":%zu,\"coefficientComparisons\":%zu,\"rawProductComparisons\":%zu,\"timed\":false}\n",n,n*128*2,n*2);
    return 0;
}
