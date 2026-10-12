#include "../../include/packed131.h"
#include "direct_sigma_half5.generated.h"

using eccPacked131::P131;

#if defined(__clang__) || defined(__GNUC__)
#define DIRECT_NOINLINE __attribute__((noinline, flatten))
#else
#define DIRECT_NOINLINE
#endif

extern "C" {
DIRECT_NOINLINE P131 half5_j3(P131 p) { return eccDirectSigmaHalf5::apply5_j3(p); }
DIRECT_NOINLINE P131 half5_j4(P131 p) { return eccDirectSigmaHalf5::apply5_j4(p); }
DIRECT_NOINLINE P131 half5_j5(P131 p) { return eccDirectSigmaHalf5::apply5_j5(p); }
DIRECT_NOINLINE P131 half5_j6(P131 p) { return eccDirectSigmaHalf5::apply5_j6(p); }
DIRECT_NOINLINE P131 half5_j7(P131 p) { return eccDirectSigmaHalf5::apply5_j7(p); }
DIRECT_NOINLINE P131 half5_j8(P131 p) { return eccDirectSigmaHalf5::apply5_j8(p); }
DIRECT_NOINLINE P131 half5_j9(P131 p) { return eccDirectSigmaHalf5::apply5_j9(p); }
DIRECT_NOINLINE P131 half5_j10(P131 p) { return eccDirectSigmaHalf5::apply5_j10(p); }
}
