/* Matched serialized/pipelined stage check. The CPU oracle is deliberately
 * outside timed work; generation, transfer, reduction, result checking and
 * allocation are inside. No end-to-end ECDLP operation claim is made. */
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <limits>
#include "params.h"
#include "pipeline.cuh"

static size_t number(const char* text) {
    if (!*text || *text == '-') throw std::invalid_argument("expected unsigned integer");
    char* end = nullptr; errno = 0;
    auto value = std::strtoull(text, &end, 0);
    if (errno || *end || value > std::numeric_limits<size_t>::max())
        throw std::invalid_argument("invalid integer");
    return size_t(value);
}
int main(int argc, char** argv) {
    try {
        if (argc != 6) {
            std::fprintf(stderr, "usage: %s jobs batch slots seed min_gpu_batch\n", argv[0]);
            return 2;
        }
        const size_t jobs = number(argv[1]), batch = number(argv[2]);
        const size_t slots = number(argv[3]), seed = number(argv[4]), cutoff = number(argv[5]);
        if (jobs > 65536 || slots > 2) throw std::invalid_argument("jobs/slots out of range");
        const int rows = MAC_ROWS, cols = MAC_COLS;
        const size_t cells = size_t(rows)*cols;
        struct Expected { uint64_t digest; int rank; std::vector<int> piv; };
        std::vector<Expected> expected(jobs);
        std::vector<bool> seen(jobs, false);
        std::vector<uint32_t> one(cells);
        auto prepare = [&](size_t id, uint32_t* a) {
            // Exercise zero, full-rank, rank-deficient and pivot-skip cases.
            if (id % 11 == 0) std::fill(a, a+cells, 0u);
            else {
                mac_gen_matrix(a, rows, cols, seed+1000ull*id, id%3 == 0 ? 0 : std::min(4,rows-1));
                if (id % 7 == 0) for (int r = 0; r < rows; ++r) a[size_t(r)*cols] = 0;
            }
        };
        for (size_t id = 0; id < jobs; ++id) {
            prepare(id, one.data());
            auto& e = expected[id]; e.piv.assign(rows, -1);
            mac_to_mont(one.data(), cells);
            e.rank = mac_rref_serial(one.data(), rows, cols, e.piv.data());
            mac_from_mont(one.data(), cells);
            e.digest = mac_digest(one.data(), cells);
        }
        auto consume = [&](size_t id, const uint32_t* a, const int* piv, int rank) {
            if (id >= jobs || seen[id]) throw std::runtime_error("duplicate or invalid job ID");
            const auto& e = expected[id];
            if (rank != e.rank || !std::equal(e.piv.begin(), e.piv.end(), piv) ||
                mac_digest(a,cells) != e.digest) throw std::runtime_error("CPU/GPU result mismatch");
            seen[id] = true;
        };
        auto m = mac_pipeline::run(jobs,rows,cols,batch,unsigned(slots),cutoff,prepare,consume);
        if (std::count(seen.begin(), seen.end(), true) != long(jobs))
            throw std::runtime_error("missing results");
        std::printf("{\"jobs\":%zu,\"batch\":%zu,\"slots\":%zu,\"seed\":%zu,\"cutoff\":%zu,"
                    "\"p\":%u,\"rows\":%d,\"cols\":%d,\"verified\":%zu,\"cpu_jobs\":%zu,\"gpu_jobs\":%zu,"
                    "\"transfer_bytes\":%zu,\"setup_s\":%.9f,\"prepare_s\":%.9f,\"consume_s\":%.9f,"
                    "\"cpu_reduce_s\":%.9f,\"h2d_s\":%.9f,\"kernel_s\":%.9f,\"d2h_s\":%.9f,\"total_s\":%.9f}\n",
                    jobs,batch,slots,seed,cutoff,MAC_P,rows,cols,jobs,m.cpu_jobs,m.gpu_jobs,
                    m.transfer_bytes,m.setup,m.prepare,m.consume,m.cpu_reduce,m.h2d,m.kernel,m.d2h,m.total);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr,"pipeline: %s\n",e.what()); return 1;
    }
}
