// Native ownership, compact-transfer and tagged-top-word checks for the
// full-word fused-sigma shared scratch helpers. No CUDA runtime or GPU.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#ifndef ECC_BATCH
#define ECC_BATCH 16
#endif
#ifndef ECC_PACKED_COMPACT_STATE
#define ECC_PACKED_COMPACT_STATE 1
#endif

#include "../../../include/packed131.h"
#include "../../../include/packedcompactstate.cuh"
#include "../../../include/packedsigmascratch.h"

using eccPacked131::P131;

namespace {

void need(bool condition, const char *message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        std::exit(1);
    }
}

bool same(P131 a, P131 b) {
    return std::memcmp(a.v, b.v, sizeof a.v) == 0;
}

P131 value(int field, int slot, int tid, int round) {
    P131 out;
    std::uint32_t state = 0x9e3779b9u ^ std::uint32_t(field * 0x85ebca6bu) ^
        std::uint32_t(slot * 0xc2b2ae35u) ^ std::uint32_t(tid * 0x27d4eb2du) ^
        std::uint32_t(round * 0x165667b1u);
    for (int word = 0; word < 4; ++word) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        out.v[word] = state;
    }
    // Prefixes are canonical field values. Denominators retain all six bits
    // used by the compact tail, including the three jump-tag bits.
    out.v[4] = field == eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN
        ? unsigned((slot + tid + round) & 7)
        : unsigned((slot * 13 + tid * 17 + round * 7) & 63);
    return out;
}

std::size_t compactBytes(int threads) {
    return ((std::size_t(threads) + 255) / 256) * ECC_BATCH * 4352;
}

void checkCompactGuards(const std::vector<unsigned char> &bytes, std::size_t guard) {
    need(std::all_of(bytes.begin(), bytes.begin() + guard,
                     [](unsigned char value) { return value == 0xa5; }),
         "compact leading canary");
    need(std::all_of(bytes.end() - guard, bytes.end(),
                     [](unsigned char value) { return value == 0xa5; }),
         "compact trailing canary");
}

} // namespace

int main() {
    static_assert(sizeof(P131) == 20, "full-word scratch needs five words");
    static_assert(ECC_BATCH == 16, "frozen shared-scratch batch");
    constexpr std::size_t byteGuard = 64;
    constexpr std::size_t wordGuard = 16;
    const int workerCounts[] = {1, 3, 31, 32, 255, 256};
    std::size_t cases = 0, records = 0, tagged = 0;

    for (int cachedSlots : {2, 3, 4}) {
        const std::size_t scratchWords = std::size_t(
            eccPacked131::SIGMA_FUSED_SCRATCH_FIELDS) * cachedSlots * 5 * 256;
        for (int workers : workerCounts) {
            const std::size_t fieldBytes = compactBytes(workers);
            for (int round = 0; round < 8; ++round) {
                std::vector<unsigned char> chainBytes(fieldBytes + 2 * byteGuard, 0xa5);
                std::vector<unsigned char> denominatorBytes(fieldBytes + 2 * byteGuard, 0xa5);
                auto *chain = reinterpret_cast<unsigned *>(chainBytes.data() + byteGuard);
                auto *denominator = reinterpret_cast<unsigned *>(denominatorBytes.data() + byteGuard);
                need(reinterpret_cast<std::uintptr_t>(chain) % 16 == 0 &&
                     reinterpret_cast<std::uintptr_t>(denominator) % 16 == 0,
                     "compact fixture alignment");

                for (int slot = 0; slot < ECC_BATCH; ++slot) {
                    for (int tid = 0; tid < workers; ++tid) {
                        eccPacked131::compactStore131(
                            chain, slot, tid,
                            value(eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN,
                                  slot, tid, round));
                        eccPacked131::compactStore131(
                            denominator, slot, tid,
                            value(eccPacked131::SIGMA_FUSED_SCRATCH_DENOMINATOR,
                                  slot, tid, round));
                    }
                }
                checkCompactGuards(chainBytes, byteGuard);
                checkCompactGuards(denominatorBytes, byteGuard);

                std::vector<unsigned> scratch(scratchWords + 2 * wordGuard, 0xa5a5a5a5u);
                std::vector<unsigned> expected = scratch;
                unsigned *body = scratch.data() + wordGuard;
                unsigned *expectedBody = expected.data() + wordGuard;
                for (int field = 0; field < eccPacked131::SIGMA_FUSED_SCRATCH_FIELDS; ++field) {
                    const unsigned *source = field == eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN
                        ? chain : denominator;
                    for (int slot = 0; slot < cachedSlots; ++slot) {
                        for (int tid = 0; tid < workers; ++tid) {
                            const P131 stored = eccPacked131::compactLoad131(source, slot, tid);
                            eccPacked131::sigmaFusedScratchStore131(
                                body, field, slot, tid, cachedSlots, 256, stored);
                            for (int word = 0; word < 5; ++word) {
                                const std::size_t index =
                                    (((std::size_t(field) * cachedSlots + slot) * 5 + word) * 256) + tid;
                                expectedBody[index] = stored.v[word];
                            }
                        }
                    }
                }
                need(scratch == expected,
                     "shared SoA byte image, inactive lanes and canaries");

                std::vector<unsigned char> chainRoundTrip(fieldBytes + 2 * byteGuard, 0xa5);
                std::vector<unsigned char> denominatorRoundTrip(fieldBytes + 2 * byteGuard, 0xa5);
                auto *chainOut = reinterpret_cast<unsigned *>(chainRoundTrip.data() + byteGuard);
                auto *denominatorOut = reinterpret_cast<unsigned *>(
                    denominatorRoundTrip.data() + byteGuard);
                for (int field = 0; field < eccPacked131::SIGMA_FUSED_SCRATCH_FIELDS; ++field) {
                    unsigned *destination = field == eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN
                        ? chainOut : denominatorOut;
                    for (int slot = 0; slot < cachedSlots; ++slot) {
                        for (int tid = 0; tid < workers; ++tid) {
                            const P131 actual = eccPacked131::sigmaFusedScratchLoad131(
                                body, field, slot, tid, cachedSlots, 256);
                            const P131 wanted = value(field, slot, tid, round);
                            need(same(actual, wanted), "shared helper tagged-value round trip");
                            eccPacked131::compactStore131(destination, slot, tid, actual);
                            ++records;
                            if (field == eccPacked131::SIGMA_FUSED_SCRATCH_DENOMINATOR &&
                                actual.v[4] > 7)
                                ++tagged;
                        }
                    }
                }
                for (int slot = 0; slot < cachedSlots; ++slot) {
                    for (int tid = 0; tid < workers; ++tid) {
                        need(same(eccPacked131::compactLoad131(chainOut, slot, tid),
                                  value(eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN,
                                        slot, tid, round)),
                             "shared-to-compact chain conversion");
                        need(same(eccPacked131::compactLoad131(denominatorOut, slot, tid),
                                  value(eccPacked131::SIGMA_FUSED_SCRATCH_DENOMINATOR,
                                        slot, tid, round)),
                             "shared-to-compact tagged denominator conversion");
                    }
                }
                checkCompactGuards(chainRoundTrip, byteGuard);
                checkCompactGuards(denominatorRoundTrip, byteGuard);
                ++cases;
            }
        }
    }

    need(cases == 144, "native case inventory");
    need(records == 83232, "native record inventory");
    need(tagged > 30000, "tagged top-word coverage");
    std::printf("PASS: %zu native compact/shared cases, %zu field records, %zu tagged denominators; canaries and partial blocks intact\n",
                cases, records, tagged);
    return 0;
}
