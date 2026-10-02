#include "../../include/packed131.h"
#include "direct_sigma_table3.generated.h"
#include "direct_sigma_half5.generated.h"

#include <cstdint>
#include <cstdio>

using eccPacked131::P131;

namespace {

bool same(const P131 &a, const P131 &b) {
    for (int i = 0; i < 5; ++i)
        if (a.v[i] != b.v[i]) return false;
    return (a.v[4] & ~7u) == 0;
}

P131 basis(int bit) {
    P131 out{{0, 0, 0, 0, 0}};
    out.v[bit / 32] = 1u << (bit % 32);
    return out;
}

P131 oracle(P131 p, int index) {
    const P131 normal = eccPacked131::fromPolynomial131(p);
    return eccPacked131::toPolynomial131(
        eccPacked131::add131(normal, eccPacked131::sigma131(normal, index + 3)));
}

P131 half5(P131 p, int index) {
    switch (index) {
        case 0: return eccDirectSigmaHalf5::apply5_j3(p);
        case 1: return eccDirectSigmaHalf5::apply5_j4(p);
        case 2: return eccDirectSigmaHalf5::apply5_j5(p);
        case 3: return eccDirectSigmaHalf5::apply5_j6(p);
        case 4: return eccDirectSigmaHalf5::apply5_j7(p);
        case 5: return eccDirectSigmaHalf5::apply5_j8(p);
        case 6: return eccDirectSigmaHalf5::apply5_j9(p);
        default: return eccDirectSigmaHalf5::apply5_j10(p);
    }
}

uint32_t next(uint32_t &state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

}  // namespace

int main() {
    for (int index = 0; index < 8; ++index) {
        for (int bit = 0; bit < 131; ++bit) {
            const P131 input = basis(bit);
            if (!same(eccDirectSigma131::apply3(input, index), oracle(input, index))) {
                std::printf("basis mismatch: j=%d bit=%d\n", index + 3, bit);
                return 1;
            }
            if (!same(half5(input, index), oracle(input, index))) {
                std::printf("half5 basis mismatch: j=%d bit=%d\n", index + 3, bit);
                return 1;
            }
        }
    }
    uint32_t state = 0xd1eec7u;
    for (int index = 0; index < 8; ++index) {
        for (int test = 0; test < 4096; ++test) {
            P131 input{{next(state), next(state), next(state), next(state), next(state) & 7u}};
            if (!same(eccDirectSigma131::apply3(input, index), oracle(input, index))) {
                std::printf("dense mismatch: j=%d test=%d\n", index + 3, test);
                return 1;
            }
            if (!same(half5(input, index), oracle(input, index))) {
                std::printf("half5 dense mismatch: j=%d test=%d\n", index + 3, test);
                return 1;
            }
        }
    }
    std::puts("PASS: generated table3 and half5 direct maps match all 1,048 basis cases and 32,768 dense cases");
    return 0;
}
