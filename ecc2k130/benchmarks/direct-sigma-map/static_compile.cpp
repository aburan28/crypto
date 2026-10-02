#include "../../include/packed131.h"
#include "direct_sigma_table3.generated.h"

using eccPacked131::P131;

#if defined(__clang__) || defined(__GNUC__)
#define DIRECT_NOINLINE __attribute__((noinline, flatten))
#define DIRECT_INLINE inline __attribute__((always_inline))
#else
#define DIRECT_NOINLINE
#define DIRECT_INLINE inline
#endif

template <int J>
DIRECT_INLINE P131 composedFixed(P131 p) {
    const P131 normal = eccPacked131::fromPolynomial131(p);
    return eccPacked131::toPolynomial131(
        eccPacked131::add131(normal, eccPacked131::sigma131(normal, J)));
}

extern "C" {
DIRECT_NOINLINE P131 composed_j3(P131 p) { return composedFixed<3>(p); }
DIRECT_NOINLINE P131 composed_j4(P131 p) { return composedFixed<4>(p); }
DIRECT_NOINLINE P131 composed_j5(P131 p) { return composedFixed<5>(p); }
DIRECT_NOINLINE P131 composed_j6(P131 p) { return composedFixed<6>(p); }
DIRECT_NOINLINE P131 composed_j7(P131 p) { return composedFixed<7>(p); }
DIRECT_NOINLINE P131 composed_j8(P131 p) { return composedFixed<8>(p); }
DIRECT_NOINLINE P131 composed_j9(P131 p) { return composedFixed<9>(p); }
DIRECT_NOINLINE P131 composed_j10(P131 p) { return composedFixed<10>(p); }

DIRECT_NOINLINE P131 table3_j3(P131 p) { return eccDirectSigma131::apply3_j3(p); }
DIRECT_NOINLINE P131 table3_j4(P131 p) { return eccDirectSigma131::apply3_j4(p); }
DIRECT_NOINLINE P131 table3_j5(P131 p) { return eccDirectSigma131::apply3_j5(p); }
DIRECT_NOINLINE P131 table3_j6(P131 p) { return eccDirectSigma131::apply3_j6(p); }
DIRECT_NOINLINE P131 table3_j7(P131 p) { return eccDirectSigma131::apply3_j7(p); }
DIRECT_NOINLINE P131 table3_j8(P131 p) { return eccDirectSigma131::apply3_j8(p); }
DIRECT_NOINLINE P131 table3_j9(P131 p) { return eccDirectSigma131::apply3_j9(p); }
DIRECT_NOINLINE P131 table3_j10(P131 p) { return eccDirectSigma131::apply3_j10(p); }
}
