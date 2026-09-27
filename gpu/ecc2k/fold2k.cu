/* fold2k.cu -- build a folded pair table on the GPU.
 *
 *   fold2k PLAN TABLE [--chunk C] [--block B] [--max-threads T]
 *
 * Reads a plan (`examples/dump_fold_plan.rs`), builds the folded table
 * `PairSumTable::build_folded_within` would build from it -- same buckets,
 * same tagged words, same presence filter -- and writes it where
 * `examples/load_fold_table.rs` loads it and checks it against the CPU's
 * own.  `make device-roundtrip` does all three.
 *
 * Everything here but the CUDA calls is shared with the host emulation
 * (`test_pairtable_emu.cpp`), which runs the same kernels on concurrent
 * CPU threads under ThreadSanitizer and compares the result with the
 * CPU's table: the plan reader and table writer (`fold_io.hpp`), the
 * chunking and the scan (`pt_fold_chunks`, `pt_fold_scan`), the geometry
 * (`pt_fold_geometry`, `pt_filter_bits`), and the kernels themselves.
 * What only a device can say is whether it runs, and how fast.
 *
 * The build is two launches with a host step between, as on the CPU:
 *
 *   1. `pairtable_fold_count_kernel` counts every entry's bucket;
 *   2. the host scans the counts into offsets, which gives the total the
 *      presence filter is sized from;
 *   3. `pairtable_fold_fill_kernel` writes every tagged word at the next
 *      slot of its bucket and sets its presence bit.
 *
 * Both launches recompute every row.  The work is split into chunks of
 * `--chunk` entries of one row (default 128), each a batched inversion
 * over its chunk, grid-striding over all chunks; scratch is
 * `pt_scratch_elems(chunk)` field elements per thread, so a launch's
 * scratch is `threads * 2 * chunk` elements whatever the base.
 *
 * Timing is by CUDA events around each launch, and by the clock around
 * the host step and the transfers.  It is wall-clock on one device: a
 * practicality figure, not the repository's metric (AGENTS.md §6), and
 * the curve additions it performs -- twice the stored entries -- are
 * printed beside it so that it can be put in operation units.
 */
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include "fold_io.hpp"
#include "pairtable.cuh"

#define CUDA_OK(call)                                                                   \
    do {                                                                                \
        const cudaError_t e_ = (call);                                                  \
        if (e_ != cudaSuccess) {                                                        \
            fprintf(stderr, "%s:%d: %s: %s\n", __FILE__, __LINE__, #call,               \
                    cudaGetErrorString(e_));                                            \
            exit(1);                                                                    \
        }                                                                               \
    } while (0)

template <class T>
static T *upload(const std::vector<T> &host) {
    T *dev = nullptr;
    CUDA_OK(cudaMalloc(&dev, (host.empty() ? 1 : host.size()) * sizeof(T)));
    if (!host.empty()) {
        CUDA_OK(cudaMemcpy(dev, host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice));
    }
    return dev;
}

template <class T>
static std::vector<T> download(const T *dev, size_t count) {
    std::vector<T> host(count);
    if (count) CUDA_OK(cudaMemcpy(host.data(), dev, count * sizeof(T), cudaMemcpyDeviceToHost));
    return host;
}

static double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0)
        .count();
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s PLAN TABLE [--chunk C] [--block B] [--max-threads T]\n",
                argv[0]);
        return 2;
    }
    const char *plan_path = argv[1], *table_path = argv[2];
    int chunk = 128, block = 128;
    long max_threads = 1l << 17;
    for (int i = 3; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "--chunk")) chunk = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--block")) block = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--max-threads")) max_threads = atol(argv[i + 1]);
        else {
            fprintf(stderr, "unknown option %s\n", argv[i]);
            return 2;
        }
    }
    if (chunk < 1 || block < 1 || max_threads < 1) {
        fprintf(stderr, "chunk, block and max-threads must be positive\n");
        return 2;
    }

    PtFoldPlan plan;
    std::string err;
    if (!pt_read_plan(plan_path, plan, err)) {
        fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    if (plan.degree != (uint32_t)F2M_M) {
        fprintf(stderr, "%s is a plan for n = %u; this binary is built for n = %d\n", plan_path,
                plan.degree, F2M_M);
        return 1;
    }
    const int n_points = (int)plan.by_orbit.size(), n_reps = (int)plan.rep_pts.size();
    PtFoldGeometry g;
    if (!pt_fold_geometry(plan.n_orbits, n_reps, n_points, &g)) {
        fprintf(stderr, "the CPU would not build this table in the tagged form this builds\n");
        return 1;
    }
    std::vector<uint32_t> chunk_start(n_reps + 1);
    const uint64_t items = pt_fold_chunks(n_points, plan.suffix.data(), plan.rep_orbit.data(),
                                          n_reps, chunk, chunk_start.data());

    cudaDeviceProp prop;
    int device = 0;
    CUDA_OK(cudaGetDevice(&device));
    CUDA_OK(cudaGetDeviceProperties(&prop, device));

    const auto t_setup = std::chrono::steady_clock::now();
    PtFoldRows rows;
    rows.by_orbit = upload(plan.by_orbit);
    rows.n_points = n_points;
    rows.rep_pts = upload(plan.rep_pts);
    rows.rep_orbit = upload(plan.rep_orbit);
    rows.n_reps = n_reps;
    rows.suffix = upload(plan.suffix);
    rows.chunk_start = upload(chunk_start);
    rows.chunk = chunk;
    rows.items = items;
    rows.canon_tables = upload(plan.canon_tables);
    rows.canon_bytes = plan.canon_bytes;

    const uint64_t want = items < (uint64_t)max_threads ? items : (uint64_t)max_threads;
    const unsigned blocks = (unsigned)((want + block - 1) / block);
    const uint64_t threads = (uint64_t)blocks * block;
    const int stride = (int)pt_scratch_elems((size_t)chunk);
    f2e *scratch = nullptr;
    CUDA_OK(cudaMalloc(&scratch, threads * stride * sizeof(f2e)));
    uint32_t *bucket_start = nullptr;
    CUDA_OK(cudaMalloc(&bucket_start, ((size_t)g.buckets + 1) * sizeof(uint32_t)));
    CUDA_OK(cudaMemset(bucket_start, 0, ((size_t)g.buckets + 1) * sizeof(uint32_t)));
    CUDA_OK(cudaDeviceSynchronize());
    const double setup_ms = ms_since(t_setup);

    cudaEvent_t e0, e1, e2, e3;
    CUDA_OK(cudaEventCreate(&e0));
    CUDA_OK(cudaEventCreate(&e1));
    CUDA_OK(cudaEventCreate(&e2));
    CUDA_OK(cudaEventCreate(&e3));

    /* 1. Count. */
    CUDA_OK(cudaEventRecord(e0));
    pairtable_fold_count_kernel<<<blocks, block>>>(rows, g.bucket_shift, bucket_start, scratch,
                                                   stride);
    CUDA_OK(cudaGetLastError());
    CUDA_OK(cudaEventRecord(e1));
    CUDA_OK(cudaEventSynchronize(e1));

    /* 2. Scan, on the host, and size the filter from the total. */
    const auto t_scan = std::chrono::steady_clock::now();
    PtFoldTable table;
    table.bucket_shift = g.bucket_shift;
    table.bucket_start = download(bucket_start, (size_t)g.buckets + 1);
    const uint64_t total = pt_fold_scan(table.bucket_start.data(), g.buckets);
    const int filter_bits = pt_filter_bits(total);
    table.present_mask = (1ull << filter_bits) - 1;
    const size_t present_words = (size_t)((1ull << filter_bits) / 64);
    uint32_t *cursor = nullptr, *words = nullptr;
    uint64_t *present = nullptr;
    CUDA_OK(cudaMalloc(&cursor, (size_t)g.buckets * sizeof(uint32_t)));
    CUDA_OK(cudaMemcpy(cursor, table.bucket_start.data(), (size_t)g.buckets * sizeof(uint32_t),
                       cudaMemcpyHostToDevice));
    CUDA_OK(cudaMalloc(&words, (total ? total : 1) * sizeof(uint32_t)));
    CUDA_OK(cudaMalloc(&present, present_words * sizeof(uint64_t)));
    CUDA_OK(cudaMemset(present, 0, present_words * sizeof(uint64_t)));
    const double scan_ms = ms_since(t_scan);

    /* 3. Fill. */
    CUDA_OK(cudaEventRecord(e2));
    pairtable_fold_fill_kernel<<<blocks, block>>>(rows, g.bucket_shift, cursor, words, present,
                                                  table.present_mask, scratch, stride);
    CUDA_OK(cudaGetLastError());
    CUDA_OK(cudaEventRecord(e3));
    CUDA_OK(cudaEventSynchronize(e3));

    const auto t_down = std::chrono::steady_clock::now();
    table.words = download(words, total);
    table.present = download(present, present_words);
    const std::vector<uint32_t> cursor_end = download(cursor, g.buckets);
    const double download_ms = ms_since(t_down);

    float count_ms = 0, fill_ms = 0;
    CUDA_OK(cudaEventElapsedTime(&count_ms, e0, e1));
    CUDA_OK(cudaEventElapsedTime(&fill_ms, e2, e3));

    /* Every slot handed out exactly once: each cursor ends where the
     * next bucket begins.  Cheap, and the first thing a lost atomic
     * would break. */
    for (uint32_t b = 0; b < g.buckets; b++) {
        if (cursor_end[b] != table.bucket_start[b + 1]) {
            fprintf(stderr, "bucket %u: cursor ended at %u, the next bucket starts at %u\n", b,
                    cursor_end[b], table.bucket_start[b + 1]);
            return 1;
        }
    }
    if (!pt_write_table(table_path, plan, table, err)) {
        fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }

    const double adds = 2.0 * (double)total;
    printf("%s on %s (sm_%d%d): n = %d, %d points, %d rows, %llu stored words\n", plan_path,
           prop.name, prop.major, prop.minor, F2M_M, n_points, n_reps,
           (unsigned long long)total);
    printf("  launch: %llu chunks of <= %d, %u blocks x %d threads, scratch %.1f MiB\n",
           (unsigned long long)items, chunk, blocks, block,
           (double)threads * stride * sizeof(f2e) / (1 << 20));
    printf("  count %.3f ms, host scan %.3f ms, fill %.3f ms (setup %.3f ms, download %.3f ms)\n",
           count_ms, scan_ms, fill_ms, setup_ms, download_ms);
    printf("  %.0f curve additions in the two kernels, %.3g per second\n", adds,
           adds / ((count_ms + fill_ms) / 1e3));
    printf("  written to %s\n", table_path);

    cudaFree(scratch);
    cudaFree(bucket_start);
    cudaFree(cursor);
    cudaFree(words);
    cudaFree(present);
    cudaFree((void *)rows.by_orbit);
    cudaFree((void *)rows.rep_pts);
    cudaFree((void *)rows.rep_orbit);
    cudaFree((void *)rows.suffix);
    cudaFree((void *)rows.chunk_start);
    cudaFree((void *)rows.canon_tables);
    return 0;
}
