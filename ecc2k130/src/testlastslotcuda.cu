// Storage-contract checks using the actual walk kernel; no timing.
#include <cuda_runtime.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <utility>
#include "../include/packedkernels.cuh"
#if ECC_PACKED_LAST_SLOT_CACHE != 1 && ECC_PACKED_LAST_SLOT_CACHE != 2
#error "Last-slot storage probe requires cache mode 1 or 2"
#endif
#if !ECC_PACKED_COMPACT_STATE || ECC_BATCH != 16 || ECC_THREADS != 256
#error "Last-slot storage probe requires the selected compact layout"
#endif
static void need(bool ok,const char*why) { if(!ok) { std::fprintf(stderr,"FAIL: %s\n",why);std::exit(1); } }
static void checked(cudaError_t e) { if(e!=cudaSuccess) { std::fprintf(stderr,"CUDA cache test: %s\n",cudaGetErrorString(e));std::exit(1); } }
struct Counts { size_t launches=0,guardBytes=0,retainedBytes=0,canonicalFields=0,writtenScratchFields=0,comparedBytes=0; };
static Counts totals;
struct Buffer {
    static constexpr size_t guard=64;
    unsigned char*raw=nullptr;size_t size;
    std::vector<unsigned char> initial;
    explicit Buffer(std::vector<unsigned char> data):size(data.size()),initial(std::move(data)) {
        std::vector<unsigned char> all(size+2*guard,0xd3);
        std::copy(initial.begin(),initial.end(),all.begin()+guard);
        checked(cudaMalloc(&raw,all.size()));checked(cudaMemcpy(raw,all.data(),all.size(),cudaMemcpyHostToDevice));
    }
    ~Buffer() { checked(cudaFree(raw)); }
    Buffer(const Buffer&)=delete;Buffer&operator=(const Buffer&)=delete;
    unsigned char*data() { return raw+guard; }
    std::vector<unsigned char> read() {
        std::vector<unsigned char> all(size+2*guard);checked(cudaMemcpy(all.data(),raw,all.size(),cudaMemcpyDeviceToHost));
        for(size_t i=0;i<guard;++i)need(all[i]==0xd3&&all[guard+size+i]==0xd3,"buffer guard overwritten");
        totals.guardBytes+=2*guard;
        return {all.begin()+guard,all.begin()+guard+size};
    }
    void unchanged() { need(read()==initial,"walk changed suppressed-report metadata"); }
};
static size_t fieldBytes(int workers) { return size_t((workers+255)/256)*16*4352; }
// Independent byte addressing, rather than the compact load/store helper.
static size_t lowOffset(int slot,int tid) { return (size_t(tid/256)*16+slot)*4352+size_t(tid%256)*16; }
static size_t tailOffset(int slot,int tid) { return (size_t(tid/256)*16+slot)*4352+4096+tid%256; }
static unsigned randomWord(unsigned &s) { s^=s<<13;s^=s>>17;s^=s<<5;return s; }
static std::vector<unsigned char> coordinates(int workers,int pattern,unsigned seed) {
    std::vector<unsigned char> result(fieldBytes(workers),0x42);
    for(int slot=0;slot<16;++slot)for(int tid=0;tid<workers;++tid) {
        for(int word=0;word<4;++word) {
            const unsigned value=pattern?0:randomWord(seed);
            std::memcpy(result.data()+lowOffset(slot,tid)+4*word,&value,4);
        }
        result[tailOffset(slot,tid)]=pattern?0:randomWord(seed)&7;
    }
    return result;
}
// kind 0/1: x/y; 2: weighted prefix; 3: denominator.
static void inspectField(const std::vector<unsigned char>&data,int workers,int kind,
                         unsigned char poison,std::vector<unsigned char>&normalized) {
    const int padded=((workers+255)/256)*256;
    for(int slot=0;slot<16;++slot)for(int tid=0;tid<padded;++tid) {
        const bool active=tid<workers;
        const bool cached=active&&slot>=14&&(kind==3||(kind==2&&ECC_PACKED_LAST_SLOT_CACHE==2));
        const size_t low=lowOffset(slot,tid),tail=tailOffset(slot,tid);
        if(!active||cached) {
            const unsigned char expected=kind<2?0x42:poison;
            for(int i=0;i<16;++i)need(data[low+i]==expected,"padding or retained scratch field was written");
            need(data[tail]==expected,"padding or retained scratch tail was written");
            if(cached)totals.retainedBytes+=17;
        } else {
            const unsigned char mask=kind==3?63:7;
            need((data[tail]&~mask)==0,"noncanonical coordinate or unwritten scratch field");
            if(kind<2)++totals.canonicalFields;else ++totals.writtenScratchFields;
            normalized.insert(normalized.end(),data.begin()+low,data.begin()+low+16);
            normalized.push_back(data[tail]);
        }
    }
}
static std::vector<unsigned char> run(int workers,int steps,int pattern,unsigned char poison) {
    const size_t n=size_t(workers)*16,bytes=fieldBytes(workers);
    Buffer x(coordinates(workers,pattern,0x5a131001u)),y(coordinates(workers,pattern,0x6b131002u));
    Buffer chain(std::vector<unsigned char>(bytes,poison)),denom(std::vector<unsigned char>(bytes,poison));
    std::vector<unsigned char> deadData(n*sizeof(unsigned));
    for(size_t i=0;i<n;++i) { const unsigned one=1;std::memcpy(deadData.data()+i*sizeof(unsigned),&one,sizeof(one)); }
    Buffer dead(std::move(deadData)),seed(std::vector<unsigned char>(n*sizeof(unsigned long long),0));
    Buffer start(std::vector<unsigned char>(n*sizeof(unsigned long long),0));
    Buffer count(std::vector<unsigned char>(3*sizeof(unsigned),0)),dp(std::vector<unsigned char>(sizeof(DpRecord),0xe1));
    WalkParams<unsigned> p{};p.threads=workers;p.steps=steps;p.dpWeight=-1;p.dpCap=1;p.runId=1567;
    p.x=reinterpret_cast<unsigned*>(x.data());p.y=reinterpret_cast<unsigned*>(y.data());
    p.pchain=reinterpret_cast<unsigned*>(chain.data());p.dead=reinterpret_cast<unsigned*>(dead.data());
    p.seed=reinterpret_cast<unsigned long long*>(seed.data());p.startIter=reinterpret_cast<unsigned long long*>(start.data());
    p.dpCount=reinterpret_cast<unsigned*>(count.data());p.dp=reinterpret_cast<DpRecord*>(dp.data());
    const unsigned blocks=(workers+127)/128+1; // Always include a fully inactive block.
    std::vector<unsigned char> normalized;
    for(int launch=0;launch<2;++launch) {
        eccPacked131::walk<<<blocks,256>>>(p,reinterpret_cast<unsigned*>(denom.data()));
        checked(cudaGetLastError());checked(cudaDeviceSynchronize());++totals.launches;
        inspectField(x.read(),workers,0,poison,normalized);inspectField(y.read(),workers,1,poison,normalized);
        inspectField(chain.read(),workers,2,poison,normalized);inspectField(denom.read(),workers,3,poison,normalized);
        dead.unchanged();seed.unchanged();start.unchanged();count.unchanged();dp.unchanged();
        p.iterBase+=steps;
    }
    return normalized;
}
int main() {
    int scenarios=0;
    for(int workers:{0,1,31,127,128,129,255,256,257,511,512,513,1025})
        for(int steps:{1,2,7})for(int pattern=0;pattern<2;++pattern) {
            const auto first=run(workers,steps,pattern,0xa5),second=run(workers,steps,pattern,0x5a);
            need(first==second,"walk outputs depend on old scratch contents");
            totals.comparedBytes+=first.size();++scenarios;
        }
    std::printf("PASS: %d last-slot storage scenarios, %zu walk launches, two scratch poisons and two successive launches\n",scenarios,totals.launches);
    std::printf("PASS: %zu guard bytes, %zu retained scratch bytes, %zu canonical coordinate fields, %zu written scratch fields, %zu compared output bytes\n",totals.guardBytes,totals.retainedBytes,totals.canonicalFields,totals.writtenScratchFields,totals.comparedBytes);
    std::printf("PASS: dense/zero coordinates, ragged/fully inactive blocks, canonical tails, unchanged metadata and all padding\n");
}
