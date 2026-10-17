#ifndef CA_PROBE_H
#define CA_PROBE_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct ca_probe_counts {
 uint64_t add[3],sub[3],mul_plain[3],mont_mul[3],mont_square[3],mont_redc[3],inverse[3],euclid_divisions[3],complement[3];
} ca_probe_counts;
#ifdef __cplusplus
extern thread_local ca_probe_counts ca_probe_state;
extern thread_local int ca_probe_enabled;
#else
extern _Thread_local ca_probe_counts ca_probe_state;
extern _Thread_local int ca_probe_enabled;
#endif
static inline unsigned ca_probe_slot(uint64_t m){return m==72057980207847001ULL?0U:m==18014494963342201ULL?1U:2U;}
#define CA_PROBE(field,m) do { if(ca_probe_enabled) ca_probe_state.field[ca_probe_slot(m)]++; } while(0)
static inline uint64_t ca_probe_complement(uint64_t m,uint64_t x){CA_PROBE(complement,m);return m-x;}
CA_API void ca_probe_begin(void);
CA_API ca_probe_counts ca_probe_end(void);
#ifdef __cplusplus
}
#endif
#endif
