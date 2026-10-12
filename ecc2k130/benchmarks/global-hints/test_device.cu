// Exercise the production selector and resolver with forced empty/full/sparse
// histories. Expected tags come from the scalar reference, not another queue.
#include <cuda_runtime.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../../include/curveparams.h"
#include "../../include/solver.h"
#include "../../include/packedkernels.cuh"

#if !ECC_TABLE_GLOBAL_HINTS || !ECC_PACKED_COMPACT_STATE || ECC_WITNESS
#error "device queue control requires the frozen compact, witness-free global preset"
#endif
using namespace eccPacked131;
using R = Ref<CfgF131>;
static void checked(cudaError_t e) {
    if(e!=cudaSuccess) { std::fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e)); std::exit(1); }
}
static void need(bool ok) { if(!ok) { std::fprintf(stderr,"FAIL queue device control\n"); std::exit(1); } }
template<class T> static T* alloc(size_t n) { T* p; checked(cudaMalloc(&p,n*sizeof(T))); return p; }
static P131 pack(const R::Elem &e) {
    P131 p; for(int i=0;i<5;++i) p.v[i]=unsigned(e.v[i/2]>>(32*(i&1))); return p;
}
int main() {
    Solver<CfgF131> sol;
    sol.setup(eccF131::PX,eccF131::PY,eccF131::QX,eccF131::QY,
              eccF131::ELL_DEC,eccF131::S_DEC,48,1ull<<40);
    std::vector<uint32_t> table(TW_WORDS); twFillConsts(sol.walk,table.data());
    auto *dt=alloc<uint32_t>(table.size());
    checked(cudaMemcpy(dt,table.data(),table.size()*4,cudaMemcpyHostToDevice));
    if(TW_SHARED_BYTES>48*1024) {
        checked(cudaFuncSetAttribute(selectGlobalHints,cudaFuncAttributeMaxDynamicSharedMemorySize,int(TW_SHARED_BYTES)));
        checked(cudaFuncSetAttribute(resolveGlobalHints,cudaFuncAttributeMaxDynamicSharedMemorySize,int(TW_SHARED_BYTES)));
    }
    cudaFuncAttributes resolverAttrs{};
    checked(cudaFuncGetAttributes(&resolverAttrs,resolveGlobalHints));
    need(resolverAttrs.maxThreadsPerBlock==ECC_TABLE_GLOBAL_HINT_THREADS);
    need(resolverAttrs.sharedSizeBytes==0);
    int resolverBlocksPerSm=0;
    checked(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &resolverBlocksPerSm,resolveGlobalHints,ECC_TABLE_GLOBAL_HINT_THREADS,TW_SHARED_BYTES));
    need(resolverBlocksPerSm>=1);
    std::vector<R::Point> points;
    for(unsigned i=0;i<64;++i)
        points.push_back(R::startPoint(0x5eed0000ull+7919ull*i,sol.basis,sol.target,0,sol.ell,sol.spow));
    points[0]=R::scalarMul(sol.basis,u192_from(1184));
    points[1]=R::addPt(points[0],sol.walk.addend(sol.walk.rawTag(points[0],R::weight(points[0].x))));
    std::vector<unsigned> rawTags(points.size()),resolved(points.size());
    unsigned changed=0;
    for(size_t i=0;i<points.size();++i) {
        rawTags[i]=sol.walk.rawTag(points[i],R::weight(points[i].x));
        resolved[i]=sol.walk.resolveTag(points[i],rawTags[i],eccHistPush(ECC_HIST_EMPTY,rawTags[i]^ECC_TAG_EPS));
        changed+=resolved[i]!=rawTags[i];
    }
    need(changed>0);
    unsigned cases=0;
    for(int workers : {1,127,128,511,512,513,1537}) {
        const size_t n=size_t(workers)*ECC_BATCH;
        const size_t fw=compactPhysicalFieldWords(workers);
        std::vector<unsigned> x(fw),y(fw),q(n+2,0xdecafbad),seen(n);
        std::vector<unsigned long long> hist(n),expected(n),out(n);
        WalkParams<unsigned> p{}; p.threads=workers; p.steps=1; p.dpWeight=48;
        p.x=alloc<unsigned>(fw); p.y=alloc<unsigned>(fw); p.hist=alloc<unsigned long long>(n);
        p.dead=alloc<unsigned>(n); p.twConsts=dt;
        // Dead lanes still select/arithmetic exactly as production. Suppress DP
        // reporting here so this test isolates queue ownership and resolution.
        checked(cudaMemset(p.dead,1,n*4));
        auto *queue=alloc<unsigned>(n+2), *count=alloc<unsigned>(1);
        for(size_t i=0;i<n;++i) {
            const auto &pt=points[i%points.size()];
            compactStore131(x.data(),int(i/workers),int(i%workers),toPolynomial131(pack(pt.x)));
            compactStore131(y.data(),int(i/workers),int(i%workers),toPolynomial131(pack(pt.y)));
        }
        checked(cudaMemcpy(p.x,x.data(),fw*4,cudaMemcpyHostToDevice));
        checked(cudaMemcpy(p.y,y.data(),fw*4,cudaMemcpyHostToDevice));
        // Reuse the same allocation and reset after every pattern. In particular
        // full -> empty must not resolve any stale owner from the preceding run.
        for(int pattern : {0,1,0,2,3,1,0}) {
            unsigned expectedCount=0;
            for(size_t i=0;i<n;++i) {
                const unsigned raw=rawTags[i%points.size()];
                const bool hint=pattern==1 || (pattern==2 && i%401==0) ||
                    (pattern==3 && (i==0 || i+1==n));
                hist[i]=hint ? eccHistPush(ECC_HIST_EMPTY,raw^ECC_TAG_EPS) : ECC_HIST_EMPTY;
                need(eccTagFruitless(raw,hist[i],131)==hint);
                const unsigned tag=hint?resolved[i%points.size()]:raw;
                expected[i]=eccHistPush(hist[i],tag); expectedCount+=hint;
            }
            std::fill(q.begin(),q.end(),0xdecafbad);
            checked(cudaMemcpy(queue,q.data(),q.size()*4,cudaMemcpyHostToDevice));
            checked(cudaMemcpy(p.hist,hist.data(),n*8,cudaMemcpyHostToDevice));
            checked(cudaMemsetAsync(count,0,4));
            selectGlobalHints<<<(workers+ECC_THREADS-1)/ECC_THREADS,ECC_THREADS,TW_SHARED_BYTES>>>(p,queue+1,count);
            checked(cudaGetLastError());
            resolveGlobalHints<<<3,ECC_TABLE_GLOBAL_HINT_THREADS,TW_SHARED_BYTES>>>(p,queue+1,count);
            checked(cudaGetLastError()); checked(cudaDeviceSynchronize());
            unsigned actualCount; checked(cudaMemcpy(&actualCount,count,4,cudaMemcpyDeviceToHost));
            checked(cudaMemcpy(q.data(),queue,q.size()*4,cudaMemcpyDeviceToHost));
            checked(cudaMemcpy(out.data(),p.hist,n*8,cudaMemcpyDeviceToHost));
            need(actualCount==expectedCount && q.front()==0xdecafbad && q.back()==0xdecafbad);
            std::fill(seen.begin(),seen.end(),0);
            for(unsigned j=0;j<actualCount;++j) {
                const unsigned owner=q[j+1]; need(owner<n && !seen[owner]++);
            }
            for(size_t i=0;i<n;++i) {
                const unsigned raw=rawTags[i%points.size()];
                need(seen[i]==unsigned(eccTagFruitless(raw,hist[i],131)) && out[i]==expected[i]);
            }
            ++cases;
        }
        checked(cudaFree(p.x)); checked(cudaFree(p.y)); checked(cudaFree(p.hist));
        checked(cudaFree(p.dead)); checked(cudaFree(queue)); checked(cudaFree(count));
    }
    checked(cudaFree(dt));
    std::printf("PASS: %u production selector/resolver device cases; empty/full/sparse, canaries, partial blocks, repeated resets and scalar-reference histories; resolver %d threads, launch max %d, %d active block(s)/SM, %d dynamic shared bytes\n",
                cases,ECC_TABLE_GLOBAL_HINT_THREADS,resolverAttrs.maxThreadsPerBlock,
                resolverBlocksPerSm,TW_SHARED_BYTES);
}
