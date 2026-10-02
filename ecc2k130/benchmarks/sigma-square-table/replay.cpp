#include "../../include/curveparams.h"
#include "../../include/packed131.h"

#include <cstdint>
#include <cstdio>

using R = Ref<CfgF131>;
using P = eccPacked131::P131;

static uint64_t rngState = 0x131ab123456789ULL;

static uint64_t randomWord() {
    rngState ^= rngState << 13;
    rngState ^= rngState >> 7;
    rngState ^= rngState << 17;
    return rngState;
}

static R::Elem unpack(P p) {
    const uint64_t limbs[3] = {
        p.v[0] | (uint64_t(p.v[1]) << 32),
        p.v[2] | (uint64_t(p.v[3]) << 32),
        p.v[4]
    };
    return R::fromLimbs(limbs);
}

static bool same(P a, P b) {
    for (int i = 0; i < 5; ++i)
        if (a.v[i] != b.v[i]) return false;
    return (a.v[4] & ~7u) == 0;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: replay OUTPUT\n");
        return 2;
    }
    std::FILE *out = std::fopen(argv[1], "wb");
    if (!out) {
        std::perror(argv[1]);
        return 2;
    }

    uint32_t table[eccPacked131::SQ_TAB_WORDS];
    eccPacked131::fillSquareTable131(table);
    const P edges[] = {P{}, P{{1u, 0u, 0u, 0u, 0u}},
                       P{{~0u, ~0u, ~0u, ~0u, 7u}}};
    constexpr int cases = 131 + 3 + 20000;
    for (int test = 0; test < cases; ++test) {
        P input{};
        if (test < 131) {
            input.v[test / 32] = 1u << (test % 32);
        } else if (test < 134) {
            input = edges[test - 131];
        } else {
            for (uint32_t &word : input.v) word = uint32_t(randomWord());
            input.v[4] &= 7u;
        }

        const P want = eccPacked131::squarePolynomial131(input);
        const P got = eccPacked131::sigmaLambdaSquare131(input, table);
        const R::Elem reference = R::sqr(unpack(eccPacked131::fromPolynomial131(input)));
        const bool referenceMatch =
            unpack(eccPacked131::fromPolynomial131(got)) == reference;
        if (!same(got, want) || !referenceMatch) {
            std::fprintf(stderr,
                         "sigma lambda square mismatch at case %d, flag %d, "
                         "got %08x:%08x:%08x:%08x:%08x, "
                         "want %08x:%08x:%08x:%08x:%08x, ref %d\n",
                         test, ECC_SIGMA_SQUARE_TABLE,
                         got.v[4], got.v[3], got.v[2], got.v[1], got.v[0],
                         want.v[4], want.v[3], want.v[2], want.v[1], want.v[0],
                         referenceMatch ? 1 : 0);
            std::fclose(out);
            return 1;
        }
        if (std::fwrite(got.v, sizeof(got.v), 1, out) != 1) {
            std::perror("write replay");
            std::fclose(out);
            return 2;
        }
    }
    if (std::fclose(out) != 0) {
        std::perror("close replay");
        return 2;
    }
    std::printf("{\"valid\":true,\"sigmaSquareTable\":%d,\"cases\":%d,"
                "\"streamBytes\":%zu}\n",
                ECC_SIGMA_SQUARE_TABLE, cases, size_t(cases) * sizeof(P));
    return 0;
}
