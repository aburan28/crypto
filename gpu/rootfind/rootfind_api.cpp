/* rootfind_api.cpp -- the C ABI (rootfind_api.h): argument validation,
 * canonical -> Montgomery conversion, table construction, and dispatch to
 * the emulation backend or (when built with -DRF_HAVE_CUDA=1 and linked
 * against rootfind_cuda.o) the CUDA backend.
 *
 * Built two ways:
 *   make lib        -> librootfind_oracle.a   (emulation only; what cargo
 *                      compiles via build.rs for the `gpu-rootfind` feature)
 *   make lib-cuda   -> librootfind_cuda.a     (emulation + CUDA; needs nvcc)
 */
#include <cstring>
#include <string>

#include "rootfind_host.hpp"

#ifndef RF_HAVE_CUDA
#define RF_HAVE_CUDA 0
#endif

#if RF_HAVE_CUDA
/* implemented in rootfind_cuda.cu */
extern "C" int rf_cuda_decompose64(const Fp64Ctx *f,
                                   const uint64_t *fb_x, const uint64_t *fb_y, uint32_t n_fb,
                                   const uint64_t *tx, const uint64_t *ty, uint32_t n_targets,
                                   int lookup, uint32_t chunk,
                                   rf_relation *out, uint32_t cap, uint32_t *n_found, double *seconds);
extern "C" int rf_cuda_decompose256(const Fp256Ctx *f,
                                    const U256x64 *fb_x, const U256x64 *fb_y, uint32_t n_fb,
                                    const U256x64 *tx, const U256x64 *ty, uint32_t n_targets,
                                    int lookup, uint32_t chunk,
                                    rf_relation *out, uint32_t cap, uint32_t *n_found, double *seconds);
extern "C" int rf_cuda_device_present(char *name, int name_len);
#endif

namespace {

template <class F>
int run_emulated(const typename F::ctx &f,
                 const std::vector<typename F::elt> &fbx, const std::vector<typename F::elt> &fby,
                 const std::vector<typename F::elt> &tx, const std::vector<typename F::elt> &ty,
                 int lookup, uint32_t chunk, rf_relation *out, uint32_t cap, rf_result *res) {
    DecompArgs<F> a;
    a.f = f;
    a.fb_x = fbx.data(); a.fb_y = fby.data(); a.n_fb = (uint32_t)fbx.size();
    a.t_x = tx.data(); a.t_y = ty.data(); a.n_targets = (uint32_t)tx.size();
    a.out = out; a.out_cap = cap;
    uint32_t count = 0;
    a.out_count = &count;

    double t0 = rf::now_seconds();
    if (lookup == RF_LOOKUP_HASH) {
        rf::HashTableHost h = rf::build_hash<F>(f, fbx.data(), a.n_fb);
        HashTable view = h.view();
        t0 = rf::now_seconds();
        rf::emulate_chunk<F, HashTable>(a, view, chunk);
    } else {
        rf::SortedTableHost s = rf::build_sorted<F>(f, fbx.data(), a.n_fb);
        SortedTable view = s.view();
        t0 = rf::now_seconds();
        rf::emulate_chunk<F, SortedTable>(a, view, chunk);
    }
    res->seconds = rf::now_seconds() - t0;
    res->n_found = count;
    res->n_written = count < cap ? count : cap;
    res->backend_used = RF_BACKEND_EMULATE;
    res->lookup_used = lookup;
    res->chunk_used = chunk;
    return count > cap ? RF_ERR_OVERFLOW : RF_OK;
}

int resolve_backend(int requested) {
    if (requested == RF_BACKEND_AUTO) return rf_cuda_available() ? RF_BACKEND_CUDA : RF_BACKEND_EMULATE;
    return requested;
}

int check_common(const void *fb_x, const void *fb_y, uint32_t n_fb, const void *tx, const void *ty,
                 uint32_t n_targets, rf_relation *out, uint32_t cap, rf_options *opt, rf_result *res) {
    if (!fb_x || !fb_y || !tx || !ty || !res || n_fb == 0 || n_targets == 0) return RF_ERR_BAD_ARG;
    if (cap > 0 && !out) return RF_ERR_BAD_ARG;
    std::memset(res, 0, sizeof(*res));
    if (opt->chunk == 0) opt->chunk = 32;
    if (!rf::chunk_ok(opt->chunk)) return RF_ERR_BAD_CHUNK;
    if (opt->lookup != RF_LOOKUP_AUTO && opt->lookup != RF_LOOKUP_SORTED && opt->lookup != RF_LOOKUP_HASH) return RF_ERR_BAD_ARG;
    if (opt->backend != RF_BACKEND_AUTO && opt->backend != RF_BACKEND_EMULATE && opt->backend != RF_BACKEND_CUDA) return RF_ERR_BAD_ARG;
    return RF_OK;
}

} /* namespace */

extern "C" int rf64_decompose(uint64_t p,
                              const uint64_t *fb_x, const uint64_t *fb_y, uint32_t n_fb,
                              const uint64_t *tx, const uint64_t *ty, uint32_t n_targets,
                              rf_relation *out, uint32_t cap,
                              const rf_options *opt_in, rf_result *res) {
    rf_options opt = opt_in ? *opt_in : rf_options{RF_BACKEND_AUTO, RF_LOOKUP_AUTO, 0};
    int rc = check_common(fb_x, fb_y, n_fb, tx, ty, n_targets, out, cap, &opt, res);
    if (rc != RF_OK) return rc;
    if ((p & 1) == 0 || p <= 3 || p >= (1ull << 62)) return RF_ERR_BAD_ARG;
    for (uint32_t i = 0; i < n_fb; i++) if (fb_x[i] >= p || fb_y[i] >= p) return RF_ERR_BAD_ARG;
    for (uint32_t k = 0; k < n_targets; k++) if (tx[k] >= p || ty[k] >= p) return RF_ERR_BAD_ARG;

    Fp64Ctx f = fp64_ctx_new(p);
    std::vector<uint64_t> fbx = rf::to_mont64(f, fb_x, n_fb), fby = rf::to_mont64(f, fb_y, n_fb);
    std::vector<uint64_t> mx = rf::to_mont64(f, tx, n_targets), my = rf::to_mont64(f, ty, n_targets);
    int lookup = rf::resolve_lookup(opt.lookup, n_fb);
    int backend = resolve_backend(opt.backend);

    if (backend == RF_BACKEND_CUDA) {
#if RF_HAVE_CUDA
        uint32_t found = 0;
        double secs = 0;
        rc = rf_cuda_decompose64(&f, fbx.data(), fby.data(), n_fb, mx.data(), my.data(), n_targets,
                                 lookup, opt.chunk, out, cap, &found, &secs);
        res->n_found = found;
        res->n_written = found < cap ? found : cap;
        res->backend_used = RF_BACKEND_CUDA;
        res->lookup_used = lookup;
        res->chunk_used = opt.chunk;
        res->seconds = secs;
        if (rc != RF_OK) return rc;
        return found > cap ? RF_ERR_OVERFLOW : RF_OK;
#else
        return RF_ERR_NO_CUDA;
#endif
    }
    return run_emulated<FieldFp64>(f, fbx, fby, mx, my, lookup, opt.chunk, out, cap, res);
}

extern "C" int rf256_decompose(const uint64_t *fb_x, const uint64_t *fb_y, uint32_t n_fb,
                               const uint64_t *tx, const uint64_t *ty, uint32_t n_targets,
                               rf_relation *out, uint32_t cap,
                               const rf_options *opt_in, rf_result *res) {
    rf_options opt = opt_in ? *opt_in : rf_options{RF_BACKEND_AUTO, RF_LOOKUP_AUTO, 0};
    int rc = check_common(fb_x, fb_y, n_fb, tx, ty, n_targets, out, cap, &opt, res);
    if (rc != RF_OK) return rc;

    Fp256Ctx f = fp256_ctx_p256();
    /* canonical bound check: every coordinate < p */
    auto below_p = [&](const uint64_t *l) {
        for (int i = 3; i >= 0; i--) {
            if (l[i] < f.p.v[i]) return true;
            if (l[i] > f.p.v[i]) return false;
        }
        return false;
    };
    for (uint32_t i = 0; i < n_fb; i++) if (!below_p(fb_x + 4 * i) || !below_p(fb_y + 4 * i)) return RF_ERR_BAD_ARG;
    for (uint32_t k = 0; k < n_targets; k++) if (!below_p(tx + 4 * k) || !below_p(ty + 4 * k)) return RF_ERR_BAD_ARG;

    std::vector<U256x64> fbx = rf::to_mont256(f, fb_x, n_fb), fby = rf::to_mont256(f, fb_y, n_fb);
    std::vector<U256x64> mx = rf::to_mont256(f, tx, n_targets), my = rf::to_mont256(f, ty, n_targets);
    int lookup = rf::resolve_lookup(opt.lookup, n_fb);
    int backend = resolve_backend(opt.backend);

    if (backend == RF_BACKEND_CUDA) {
#if RF_HAVE_CUDA
        uint32_t found = 0;
        double secs = 0;
        rc = rf_cuda_decompose256(&f, fbx.data(), fby.data(), n_fb, mx.data(), my.data(), n_targets,
                                  lookup, opt.chunk, out, cap, &found, &secs);
        res->n_found = found;
        res->n_written = found < cap ? found : cap;
        res->backend_used = RF_BACKEND_CUDA;
        res->lookup_used = lookup;
        res->chunk_used = opt.chunk;
        res->seconds = secs;
        if (rc != RF_OK) return rc;
        return found > cap ? RF_ERR_OVERFLOW : RF_OK;
#else
        return RF_ERR_NO_CUDA;
#endif
    }
    return run_emulated<FieldP256>(f, fbx, fby, mx, my, lookup, opt.chunk, out, cap, res);
}

extern "C" int rf_cuda_available(void) {
#if RF_HAVE_CUDA
    return rf_cuda_device_present(nullptr, 0);
#else
    return 0;
#endif
}

extern "C" const char *rf_backend_info(void) {
    static std::string info;
#if RF_HAVE_CUDA
    char name[256] = {0};
    if (rf_cuda_device_present(name, (int)sizeof(name))) {
        info = std::string("cuda: ") + name + " (+ host emulation)";
    } else {
        info = "cuda backend linked, no device present; host emulation only";
    }
#else
    info = "host emulation only (built without RF_HAVE_CUDA)";
#endif
    return info.c_str();
}

extern "C" const char *rf_strerror(int status) {
    switch (status) {
        case RF_OK: return "ok";
        case RF_ERR_OVERFLOW: return "output buffer overflow (n_found > cap)";
        case RF_ERR_BAD_ARG: return "bad argument (null pointer, empty input, coordinate >= p, or p out of range)";
        case RF_ERR_NO_CUDA: return "CUDA backend not available";
        case RF_ERR_CUDA: return "CUDA runtime error";
        case RF_ERR_BAD_CHUNK: return "chunk must be 16, 32 or 64";
        case RF_ERR_NOMEM: return "out of memory";
        default: return "unknown status";
    }
}
