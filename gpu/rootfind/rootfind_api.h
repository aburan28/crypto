/* rootfind_api.h -- C ABI for the GPU root-finding 2-decomposition oracle.
 *
 * Linked from Rust through src/cryptanalysis/gpu_rootfind.rs (cargo
 * feature `gpu-rootfind`), from the CUDA driver bench.cu, and from the
 * CPU test harness.  Two backends implement the same contract:
 *
 *   emulate -- the kernel's per-thread function (rootfind.cuh) run on the
 *              host for every (target, slice) thread index.  Always
 *              available; what this repository's validation ran.
 *   cuda    -- rootfind_cuda.cu, compiled only when the library is built
 *              with RF_HAVE_CUDA (make lib-cuda).  rf_cuda_available()
 *              reports whether it is linked in *and* a device is present.
 *
 * Coordinates cross this boundary in canonical form (integers < p,
 * little-endian limbs for P-256); the library converts to Montgomery
 * form internally.  Results are index-level relations:
 *
 *   target k:  R_k = eps * F_i + tau * F_j          (tau = +-1, j >= i)
 *              R_k = eps * F_i                      (tau = 0, j == i)
 *
 * Every result must still be re-added in the group by the caller (the
 * Rust harness does); the library itself only promises to report what
 * the kernel computed.
 */
#ifndef GPU_ROOTFIND_API_H
#define GPU_ROOTFIND_API_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rf_relation {
    uint32_t target;
    uint32_t i;
    uint32_t j;
    int8_t eps;
    int8_t tau;
    uint8_t pad[2];
} rf_relation;   /* 16 bytes */

enum rf_backend { RF_BACKEND_AUTO = 0, RF_BACKEND_EMULATE = 1, RF_BACKEND_CUDA = 2 };
enum rf_lookup  { RF_LOOKUP_AUTO = 0, RF_LOOKUP_SORTED = 1, RF_LOOKUP_HASH = 2 };

typedef struct rf_options {
    int backend;        /* enum rf_backend */
    int lookup;         /* enum rf_lookup */
    uint32_t chunk;     /* factor-base slice per thread: 16, 32 or 64; 0 = 32 */
} rf_options;

typedef struct rf_result {
    uint32_t n_written;   /* relations stored in out[] (<= cap) */
    uint32_t n_found;     /* relations the kernel found; > cap means overflow */
    int backend_used;     /* enum rf_backend actually run */
    int lookup_used;      /* enum rf_lookup actually run */
    uint32_t chunk_used;
    double seconds;       /* wall time of the oracle call, excluding table build */
} rf_result;

enum rf_status {
    RF_OK = 0,
    RF_ERR_OVERFLOW = 1,          /* out[] too small; n_found tells how many */
    RF_ERR_BAD_ARG = -1,
    RF_ERR_NO_CUDA = -2,          /* CUDA backend requested but not linked / no device */
    RF_ERR_CUDA = -3,             /* a CUDA runtime call failed */
    RF_ERR_BAD_CHUNK = -4,
    RF_ERR_NOMEM = -5
};

/* 64-bit field: odd prime 3 < p < 2^62 (same bound as prime_fast::Fp64).
 * fb_x/fb_y: n_fb canonical coordinates; tx/ty: n_targets canonical
 * coordinates.  Writes up to `cap` relations to out[]. */
int rf64_decompose(uint64_t p,
                   const uint64_t *fb_x, const uint64_t *fb_y, uint32_t n_fb,
                   const uint64_t *tx, const uint64_t *ty, uint32_t n_targets,
                   rf_relation *out, uint32_t cap,
                   const rf_options *opt, rf_result *res);

/* NIST P-256 field (fixed modulus).  Each coordinate is 4 little-endian
 * 64-bit limbs, so fb_x has 4*n_fb words etc. */
int rf256_decompose(const uint64_t *fb_x, const uint64_t *fb_y, uint32_t n_fb,
                    const uint64_t *tx, const uint64_t *ty, uint32_t n_targets,
                    rf_relation *out, uint32_t cap,
                    const rf_options *opt, rf_result *res);

/* 1 if the CUDA backend is compiled in and at least one device answers. */
int rf_cuda_available(void);

/* Human-readable backend description, e.g. "emulate (host)" or the
 * device name when CUDA is available. */
const char *rf_backend_info(void);

const char *rf_strerror(int status);

#ifdef __cplusplus
}
#endif

#endif /* GPU_ROOTFIND_API_H */
