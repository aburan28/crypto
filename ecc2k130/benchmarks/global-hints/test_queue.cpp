#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

// A phase/ownership model. It proves neither field arithmetic nor GPU memory
// ordering; those need the source audit and matched device replay. The cold
// transform depends only on its owner's selected history, so changing which
// CUDA lane evaluates it may change order but cannot change its value.
static void need(bool ok) { if (!ok) throw std::runtime_error("queue invariant"); }
static uint32_t cold(uint32_t selected, uint32_t point) {
    return (selected & 0xffff0000u) | ((selected ^ point ^ 0x55aa) & 0xffffu);
}
static void check(int workers, int batch, const std::vector<bool> &hints, bool reverse) {
    const size_t count = size_t(workers) * batch;
    need(count <= std::numeric_limits<uint32_t>::max() && hints.size() == count);
    std::vector<uint32_t> selected(count), expected(count), actual(count), queue;
    std::vector<unsigned char> visits(count);
    for (size_t i=0;i<count;++i) {
        selected[i] = uint32_t(i * 31337 + 7);
        expected[i] = hints[i] ? cold(selected[i],uint32_t(i*17)) : selected[i];
    }
    // Emulate arbitrary worker-range reservation order, including a reversed
    // block scheduling order. Within a worker, selected slots stay ascending.
    for (int w=0;w<workers;++w) {
        const int owner = reverse ? workers-1-w : w;
        for (int slot=0;slot<batch;++slot) {
            const uint32_t key = uint32_t(size_t(slot)*workers + owner);
            if(hints[key]) queue.push_back(key);
        }
    }
    need(queue.size() <= count);
    actual=selected;
    const size_t stride=128*188;
    for(size_t lane=0;lane<std::min(stride,queue.size());++lane)
        for(size_t pos=lane;pos<queue.size();pos+=stride) {
            const uint32_t key=queue[pos];
            need(key<count && !visits[key]++ && hints[key]);
            const unsigned slot=key/unsigned(workers), tid=key%unsigned(workers);
            need(size_t(slot)*workers+tid==key);
            actual[key]=cold(selected[key],uint32_t(key*17));
        }
    need(actual==expected);
    for(size_t i=0;i<count;++i) need(visits[i]==unsigned(hints[i]));
    // The hot pass runs only after all selected histories are resolved, then
    // the next launch reinitializes the counter/queue. This must hold even for
    // empty and full queues.
    for(size_t i=0;i<count;++i) need((actual[i]+uint32_t(i))==(expected[i]+uint32_t(i)));
}
int main() {
    try {
        uint64_t cases=0;
        for(int workers=1;workers<=4;++workers)
            for(int batch=1;batch<=4;++batch) {
                const int n=workers*batch;
                std::vector<bool> hints(n);
                for(uint32_t mask=0;mask<(1u<<n);++mask) {
                    for(int i=0;i<n;++i) hints[i]=(mask>>i)&1u;
                    check(workers,batch,hints,false); check(workers,batch,hints,true);
                    cases+=2;
                }
            }
        for(int workers : {511,512,513,96256}) {
            std::vector<bool> hints(size_t(workers)*16);
            for(int pattern=0;pattern<4;++pattern) {
                for(size_t i=0;i<hints.size();++i)
                    hints[i]=pattern==1 || (pattern==2 && i%401==0) ||
                             (pattern==3 && (i==0 || i+1==hints.size()));
                check(workers,16,hints,false); check(workers,16,hints,true); cases+=2;
            }
        }
        // The resolver must use size_t induction. An unsigned pos+stride can
        // wrap near UINT32_MAX even when key/count setup is valid.
        const size_t last=std::numeric_limits<uint32_t>::max()-10ull;
        need(last+128*188 > std::numeric_limits<uint32_t>::max());
        std::cout << "PASS: " << cases << " ownership/phase cases; empty/full/sparse,"
                     " partial-block owners, production population and 32-bit stride boundary\n";
        return 0;
    } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
