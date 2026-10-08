/* Bounded local CPU/CUDA pipeline. Callbacks own symbolic construction and
 * result processing; only regular dense reductions are sent to the device.
 * Completion callbacks receive stable job IDs and may arrive out of order.
 * This is pinned PCIe staging, not a network/GPUDirect RDMA transport. */
#pragma once
#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "kernels.cuh"

namespace mac_pipeline {
using Clock = std::chrono::steady_clock;
inline double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
inline void check(cudaError_t result) {
    if (result != cudaSuccess) throw std::runtime_error(cudaGetErrorString(result));
}
struct Metrics {
    double setup = 0, prepare = 0, consume = 0, cpu_reduce = 0;
    double h2d = 0, kernel = 0, d2h = 0, total = 0;
    size_t cpu_jobs = 0, gpu_jobs = 0, transfer_bytes = 0;
};
struct Slot {
    cudaStream_t stream = nullptr;
    cudaEvent_t events[4]{};
    uint32_t *host = nullptr, *device = nullptr;
    int *piv = nullptr, *rank = nullptr, *device_piv = nullptr, *device_rank = nullptr;
    size_t first = 0, count = 0;
    Slot() = default;
    Slot(const Slot&) = delete;
    Slot& operator=(const Slot&) = delete;
    ~Slot() {
        // Callbacks can throw while a transfer is in flight. Drain before
        // releasing pinned memory, including partially initialized slots.
        if (stream) cudaStreamSynchronize(stream);
        for (auto event : events) if (event) cudaEventDestroy(event);
        if (device) cudaFree(device);
        if (device_piv) cudaFree(device_piv);
        if (device_rank) cudaFree(device_rank);
        if (host) cudaFreeHost(host);
        if (piv) cudaFreeHost(piv);
        if (rank) cudaFreeHost(rank);
        if (stream) cudaStreamDestroy(stream);
    }
    void allocate(size_t batch, size_t cells, int rows) {
        check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
        for (auto& event : events) check(cudaEventCreate(&event));
        check(cudaMallocHost(&host, batch*cells*sizeof(uint32_t)));
        check(cudaMallocHost(&piv, batch*rows*sizeof(int)));
        check(cudaMallocHost(&rank, batch*sizeof(int)));
        check(cudaMalloc(&device, batch*cells*sizeof(uint32_t)));
        check(cudaMalloc(&device_piv, batch*rows*sizeof(int)));
        check(cudaMalloc(&device_rank, batch*sizeof(int)));
    }
};

// prepare(id, canonical_matrix); consume(id, canonical_rref, pivots, rank).
// Matrix/pivot pointers are borrowed only for the duration of the callback.
// Jobs below min_gpu_batch run on CPU, including a short final batch. No
// CUDA context or allocations are required for an entirely CPU-sized run.
// One slot gives a serialized reference; two overlap CPU callbacks and DMA
// with GPU work without allocating memory proportional to the total job count.
template<class Prepare, class Consume>
Metrics run(size_t jobs, int rows, int cols, size_t batch, unsigned slots,
            size_t min_gpu_batch, Prepare prepare, Consume consume) {
    if (rows <= 0 || rows > MAC_MAX_ROWS || cols <= 0 || cols > 65536 ||
        batch == 0 || batch > 65536 || slots == 0 || slots > 2 ||
        min_gpu_batch == 0 || min_gpu_batch > batch)
        throw std::invalid_argument("invalid shape, batch, slot count or CPU cutoff");
    const auto start = Clock::now();
    Metrics m;
    const size_t cells = size_t(rows)*cols;
    std::vector<std::unique_ptr<Slot>> pool;
    if (jobs >= min_gpu_batch) {
        const auto t = Clock::now();
        for (unsigned i = 0; i < std::min<size_t>(slots, (jobs+batch-1)/batch); ++i) {
            auto slot = std::make_unique<Slot>();
            slot->allocate(std::min(batch, jobs), cells, rows);
            pool.push_back(std::move(slot));
        }
        m.setup = seconds(t);
    }
    auto finish = [&](Slot& s) {
        if (!s.count) return;
        check(cudaEventSynchronize(s.events[3]));
        float elapsed;
        check(cudaEventElapsedTime(&elapsed, s.events[0], s.events[1])); m.h2d += elapsed/1000.;
        check(cudaEventElapsedTime(&elapsed, s.events[1], s.events[2])); m.kernel += elapsed/1000.;
        check(cudaEventElapsedTime(&elapsed, s.events[2], s.events[3])); m.d2h += elapsed/1000.;
        auto t = Clock::now();
        for (size_t i = 0; i < s.count; ++i)
            consume(s.first+i, s.host+i*cells, s.piv+i*rows, s.rank[i]);
        m.consume += seconds(t);
        s.count = 0;
    };
    size_t next_slot = 0;
    for (size_t first = 0; first < jobs; first += batch) {
        const size_t count = std::min(batch, jobs-first);
        if (count < min_gpu_batch) {
            std::vector<uint32_t> a(cells);
            std::vector<int> piv(rows);
            for (size_t i = 0; i < count; ++i) {
                auto t = Clock::now(); prepare(first+i, a.data()); m.prepare += seconds(t);
                std::fill(piv.begin(), piv.end(), -1);
                t = Clock::now();
                mac_to_mont(a.data(), cells);
                int rank = mac_rref_serial(a.data(), rows, cols, piv.data());
                mac_from_mont(a.data(), cells);
                m.cpu_reduce += seconds(t);
                t = Clock::now(); consume(first+i, a.data(), piv.data(), rank); m.consume += seconds(t);
                ++m.cpu_jobs;
            }
            continue;
        }
        Slot& s = *pool[next_slot++ % pool.size()];
        finish(s); // Backpressure: a slot cannot be reused before completion.
        auto t = Clock::now();
        for (size_t i = 0; i < count; ++i) prepare(first+i, s.host+i*cells);
        m.prepare += seconds(t);
        s.first = first; s.count = count;
        const size_t bytes = count*cells*sizeof(uint32_t);
        check(cudaEventRecord(s.events[0], s.stream));
        check(cudaMemcpyAsync(s.device, s.host, bytes, cudaMemcpyHostToDevice, s.stream));
        check(cudaMemsetAsync(s.device_piv, 0xff, count*rows*sizeof(int), s.stream));
        check(cudaEventRecord(s.events[1], s.stream));
        to_mont_kernel<<<256,256,0,s.stream>>>(s.device, count*cells);
        check(cudaGetLastError());
        rref_batch_kernel<<<count,MAC_THREADS,(MAC_MAX_ROWS+4)*sizeof(uint32_t),s.stream>>>(
            s.device, rows, cols, s.device_piv, s.device_rank);
        check(cudaGetLastError());
        from_mont_kernel<<<256,256,0,s.stream>>>(s.device, count*cells);
        check(cudaGetLastError());
        check(cudaEventRecord(s.events[2], s.stream));
        check(cudaMemcpyAsync(s.host, s.device, bytes, cudaMemcpyDeviceToHost, s.stream));
        check(cudaMemcpyAsync(s.piv, s.device_piv, count*rows*sizeof(int), cudaMemcpyDeviceToHost, s.stream));
        check(cudaMemcpyAsync(s.rank, s.device_rank, count*sizeof(int), cudaMemcpyDeviceToHost, s.stream));
        check(cudaEventRecord(s.events[3], s.stream));
        m.gpu_jobs += count;
        m.transfer_bytes += 2*bytes + count*(rows+1)*sizeof(int);
    }
    for (auto& slot : pool) finish(*slot);
    pool.clear(); // Include resource teardown in total wall time.
    m.total = seconds(start);
    return m;
}
} // namespace mac_pipeline
