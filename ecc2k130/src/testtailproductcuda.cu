#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include "../include/packed131.h"
using eccPacked131::P131;
static void checked(cudaError_t e){if(e!=cudaSuccess){std::fprintf(stderr,"CUDA raw product: %s\n",cudaGetErrorString(e));std::exit(2);}}
struct Raw { uint32_t v[9]; };
static Raw reference(P131 a,P131 b){
    Raw r{};
    for(int i=0;i<131;++i)if((a.v[i/32]>>(i%32))&1u)
        for(int j=0;j<131;++j)if((b.v[j/32]>>(j%32))&1u)
            r.v[(i+j)/32]^=1u<<((i+j)%32);
    return r;
}
extern "C" __global__ __launch_bounds__(128)
void tailProductKernel(const P131* a,const P131* b,uint32_t* output,unsigned n){
    unsigned i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n){uint32_t c[9];eccPacked131::product131(a[i],b[i],c);
#pragma unroll
        for(int w=0;w<9;++w)output[9*i+w]=c[w];
    }
}
int main(){
    std::vector<P131>a,b;
    for(int i=0;i<131;++i)for(int j=0;j<131;++j){P131 x{},y{};x.v[i/32]=1u<<(i%32);y.v[j/32]=1u<<(j%32);a.push_back(x);b.push_back(y);}
    const uint32_t patterns[]={0u,~0u,0xaaaaaaaau,0x55555555u,0x80000000u,1u,0x80000001u,0x01234567u};
    for(unsigned ah=0;ah<8;++ah)for(unsigned bh=0;bh<8;++bh)
        for(uint32_t ap:patterns)for(uint32_t bp:patterns){a.push_back(P131{{ap,ap,ap,ap,ah}});b.push_back(P131{{bp,bp,bp,bp,bh}});}
    uint32_t seed=0x2193131u;
    auto random=[&](){seed^=seed<<13;seed^=seed>>17;seed^=seed<<5;return seed;};
    for(int i=0;i<8192;++i){P131 x,y;for(int w=0;w<5;++w){x.v[w]=random();y.v[w]=random();}x.v[4]&=7u;y.v[4]&=7u;a.push_back(x);b.push_back(y);}
    const unsigned n=unsigned(a.size()),guard=256;const size_t words=size_t(n)*9;
    if(n!=29449||b.size()!=n)return 3;
    std::vector<Raw>expected;expected.reserve(n);
    for(unsigned i=0;i<n;++i)expected.push_back(reference(a[i],b[i]));
    uint64_t digest=1469598103934665603ull;
    for(const auto& row:expected)for(uint32_t word:row.v)digest=(digest^word)*1099511628211ull;
    P131 *da,*db;uint32_t *storage;
    checked(cudaMalloc(&da,n*sizeof(P131)));checked(cudaMalloc(&db,n*sizeof(P131)));checked(cudaMalloc(&storage,(words+guard*2)*4));
    checked(cudaMemcpy(da,a.data(),n*sizeof(P131),cudaMemcpyHostToDevice));checked(cudaMemcpy(db,b.data(),n*sizeof(P131),cudaMemcpyHostToDevice));
    for(uint32_t poison:{0xa5a5a5a5u,0x5a5a5a5au}){
        std::vector<uint32_t>output(words+guard*2,poison);checked(cudaMemcpy(storage,output.data(),output.size()*4,cudaMemcpyHostToDevice));
        tailProductKernel<<<(n+127)/128,128>>>(da,db,storage+guard,n);checked(cudaGetLastError());checked(cudaDeviceSynchronize());
        checked(cudaMemcpy(output.data(),storage,output.size()*4,cudaMemcpyDeviceToHost));
        for(unsigned i=0;i<n;++i)for(int w=0;w<9;++w)if(output[guard+9*i+w]!=expected[i].v[w]){std::fprintf(stderr,"MISMATCH pair=%u word=%d poison=%08x\n",i,w,poison);return 4;}
        for(unsigned i=0;i<guard;++i)if(output[i]!=poison||output[guard+words+i]!=poison){std::fprintf(stderr,"GUARD mismatch %u\n",i);return 5;}
    }
    checked(cudaFree(da));checked(cudaFree(db));checked(cudaFree(storage));
    std::printf("TAIL_PASS mode=%d pairs=29449 poison_passes=2 output_words=530082 guard_words=1024 expected_digest=%016llx\n",ECC_PACKED_TAIL_LAYOUT,(unsigned long long)digest);
    return 0;
}
