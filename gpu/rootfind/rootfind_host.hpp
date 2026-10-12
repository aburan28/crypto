/* rootfind_host.hpp -- host-side glue shared by the emulation backend,
 * the CUDA backend and the test harness: lookup-table construction,
 * canonical <-> Montgomery conversion, and the emulation driver that runs
 * rf_decompose_thread for every thread index on the CPU.
 *
 * The emulation is the *same* function the kernel runs (rootfind.cuh is
 * compiled by g++/clang++ here and by nvcc in rootfind_cuda.cu); only
 * the launch loop differs.  This mirrors the convention of gpu/ecc:
 * header-only `__host__ __device__` arithmetic, verified on the host.
 */
#ifndef GPU_ROOTFIND_HOST_HPP
#define GPU_ROOTFIND_HOST_HPP

#include <algorithm>
#include <atomic>
#include <chrono>
#include <thread>
#include <vector>

#include "rootfind.cuh"

namespace rf {

/* --- lookup tables -------------------------------------------------- */

struct SortedTableHost {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> idx;
    SortedTable view() const { SortedTable t; t.keys = keys.data(); t.idx = idx.data(); t.n = (uint32_t)keys.size(); return t; }
};

struct HashTableHost {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> vals;
    uint32_t mask = 0;
    uint32_t capacity() const { return mask + 1; }
    HashTable view() const { HashTable t; t.keys = keys.data(); t.vals = vals.data(); t.mask = mask; return t; }
};

template <class F>
SortedTableHost build_sorted(const typename F::ctx &f, const typename F::elt *fb_x, uint32_t n) {
    std::vector<std::pair<uint64_t, uint32_t>> e(n);
    for (uint32_t i = 0; i < n; i++) e[i] = std::make_pair(F::key(f, fb_x[i]), i);
    std::sort(e.begin(), e.end());
    SortedTableHost t;
    t.keys.resize(n);
    t.idx.resize(n);
    for (uint32_t i = 0; i < n; i++) { t.keys[i] = e[i].first; t.idx[i] = e[i].second; }
    return t;
}

inline uint32_t hash_capacity_for(uint32_t n) {
    uint32_t cap = 16;
    while (cap < 2 * n) cap <<= 1;   /* load factor <= 1/2 */
    return cap;
}

template <class F>
HashTableHost build_hash(const typename F::ctx &f, const typename F::elt *fb_x, uint32_t n, uint32_t cap = 0) {
    if (cap == 0) cap = hash_capacity_for(n);
    HashTableHost t;
    t.mask = cap - 1;
    t.keys.assign(cap, 0);
    t.vals.assign(cap, RF_NO_INDEX);
    for (uint32_t i = 0; i < n; i++) {
        uint64_t key = F::key(f, fb_x[i]);
        uint32_t slot = rf_hash64(key, t.mask);
        while (t.vals[slot] != RF_NO_INDEX) slot = (slot + 1) & t.mask;
        t.keys[slot] = key;
        t.vals[slot] = i;
    }
    return t;
}

/* --- conversions ---------------------------------------------------- */

inline std::vector<uint64_t> to_mont64(const Fp64Ctx &f, const uint64_t *v, uint32_t n) {
    std::vector<uint64_t> out(n);
    for (uint32_t i = 0; i < n; i++) out[i] = fp64_to_mont(f, v[i]);
    return out;
}

inline std::vector<U256x64> to_mont256(const Fp256Ctx &f, const uint64_t *limbs, uint32_t n) {
    std::vector<U256x64> out(n);
    for (uint32_t i = 0; i < n; i++) {
        U256x64 c;
        for (int l = 0; l < 4; l++) c.v[l] = limbs[4 * i + l];
        out[i] = fp256_to_mont(f, c);
    }
    return out;
}

/* --- emulation driver ----------------------------------------------- */

/* Run every thread index on the host.  Threads are split across hardware
 * threads only for speed; each tid's work is independent and the output
 * order is irrelevant (the Rust harness compares sets). */
template <class F, class T, int CHUNK>
void emulate(const DecompArgs<F> &a, const T &table, unsigned hw_threads = 0) {
    const uint32_t total = rf_thread_count(a, CHUNK);
    if (hw_threads == 0) hw_threads = std::max(1u, std::thread::hardware_concurrency());
    if (hw_threads > 1 && total >= 4096) {
        std::vector<std::thread> pool;
        std::atomic<uint32_t> next{0};
        const uint32_t grain = 1024;
        for (unsigned h = 0; h < hw_threads; h++) {
            pool.emplace_back([&]() {
                for (;;) {
                    uint32_t start = next.fetch_add(grain);
                    if (start >= total) break;
                    uint32_t end = std::min(total, start + grain);
                    for (uint32_t tid = start; tid < end; tid++) rf_decompose_thread<F, T, CHUNK>(a, table, tid);
                }
            });
        }
        for (auto &th : pool) th.join();
    } else {
        for (uint32_t tid = 0; tid < total; tid++) rf_decompose_thread<F, T, CHUNK>(a, table, tid);
    }
}

template <class F, class T>
void emulate_chunk(const DecompArgs<F> &a, const T &table, uint32_t chunk) {
    switch (chunk) {
        case 16: emulate<F, T, 16>(a, table); break;
        case 32: emulate<F, T, 32>(a, table); break;
        case 64: emulate<F, T, 64>(a, table); break;
        default: break;   /* validated by the caller */
    }
}

inline bool chunk_ok(uint32_t c) { return c == 16 || c == 32 || c == 64; }

/* Choose the lookup structure when the caller says AUTO: hash whenever
 * it would fit the device's shared-memory budget (the same rule the CUDA
 * backend applies, so emulation and device agree on the table used). */
inline int resolve_lookup(int requested, uint32_t n_fb) {
    if (requested != RF_LOOKUP_AUTO) return requested;
    return hash_capacity_for(n_fb) <= RF_SHARED_HASH_CAP ? RF_LOOKUP_HASH : RF_LOOKUP_SORTED;
}

inline double now_seconds() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

} /* namespace rf */

#endif /* GPU_ROOTFIND_HOST_HPP */
