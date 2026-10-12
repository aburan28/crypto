#pragma once

#include <array>
#include <string>

namespace sigma_fused_star {

struct Arm {
    const char *name;
    const char *makeKnob;
    int pairIlp;
    int pairClmul;
    int l2Persist;
    int slotUnroll;
    int fromReduced;
    int invPoly;
    int clmulFlat;
    int aluSquare;
    int lateY;
};

// This table is the executable form of STAR-PROTOCOL.md.  Every non-baseline
// row changes exactly one Make variable from the confirmed fused B16 baseline.
inline constexpr std::array<Arm, 11> kArms{{
    {"baseline", "",                        0, 0, 0, 1, 0, 0, 0, 0, 0},
    {"pair-ilp", "PACKED_PAIR_ILP=1",       1, 0, 0, 1, 0, 0, 0, 0, 0},
    {"pair-clmul", "PACKED_PAIR_CLMUL=1",   0, 1, 0, 1, 0, 0, 0, 0, 0},
    {"l2-persist", "PACKED_L2_PERSIST=1",   0, 0, 1, 1, 0, 0, 0, 0, 0},
    {"unroll2", "UNROLL_SLOTS=2",           0, 0, 0, 2, 0, 0, 0, 0, 0},
    {"from-reduced", "PACKED_FROM_REDUCED=1", 0, 0, 0, 1, 1, 0, 0, 0, 0},
    {"inv-poly1", "PACKED_INV_POLY=1",      0, 0, 0, 1, 0, 1, 0, 0, 0},
    {"inv-poly2", "PACKED_INV_POLY=2",      0, 0, 0, 1, 0, 2, 0, 0, 0},
    {"clmul-flat", "PACKED_CLMUL_FLAT=1",   0, 0, 0, 1, 0, 0, 1, 0, 0},
    {"alu-square", "PACKED_ALU_SQUARE=1",   0, 0, 0, 1, 0, 0, 0, 1, 0},
    {"late-y", "SIGMA_FUSED_LATE_Y=1",      0, 0, 0, 1, 0, 0, 0, 0, 1},
}};

inline const Arm *findArm(const std::string &name) {
    for (const Arm &arm : kArms)
        if (name == arm.name) return &arm;
    return nullptr;
}

inline constexpr long long kUpdatesPerSample = 201863462912LL;
inline constexpr double kMaxAaDrift = 0.01;
inline constexpr double kMinGeometricMean = 1.015;
inline constexpr double kMinPairRatio = 1.005;
inline constexpr double kNoiseMargin = 0.005;
inline constexpr double kGoalMillionPerSecond = 26000.0;

}  // namespace sigma_fused_star
