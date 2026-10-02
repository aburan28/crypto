// Exact host-only screen for a point-dependent halving/table-add rho map.
//
// This enumerates the generated GF(2^23) prime-order subgroup and its
// Frobenius/negation quotient. It does not recover a challenge logarithm or
// make a GPU-throughput claim.
#define ECC_NO_CUDA
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <queue>
#include <string>
#include <vector>

#include "../../include/curveparams.h"

using R = Ref<CfgF23>;
using Elem = R::Elem;
using Point = R::Point;

namespace {

constexpr uint32_t BAD = std::numeric_limits<uint32_t>::max();

uint64_t small(const U192 &x) {
    assert(x.v[1] == 0 && x.v[2] == 0);
    return x.v[0];
}

uint64_t mulMod(uint64_t a, uint64_t b, uint64_t modulus) {
    return uint64_t((unsigned __int128)a * b % modulus);
}

uint32_t mix32(uint32_t x) {
    x ^= x >> 16;
    x *= 0x85ebca6bu;
    x ^= x >> 13;
    x *= 0xc2b2ae35u;
    x ^= x >> 16;
    return x;
}

bool bit(const Elem &x, int i) {
    return (x.v[i >> 6] >> (i & 63)) & 1ull;
}

uint32_t necklaceWord(const Elem &x, const TableWalkConsts<23> &constants) {
    uint32_t phaseBits = 0;
    for (int coordinate = 1; coordinate <= 23; ++coordinate)
        if (bit(x, coordinate - 1)) phaseBits |= 1u << constants.L[coordinate];

    static constexpr int offsets[7][4] = {
        {0, 1, -1, -1}, {0, 2, -1, -1}, {0, 4, -1, -1},
        {0, 1, 3, -1},  {0, 1, 5, -1},  {0, 2, 7, -1},
        {0, 1, 4, 9},
    };
    uint32_t word = 0;
    for (int invariant = 0; invariant < 7; ++invariant) {
        unsigned parity = 0;
        for (int i = 0; i < 23; ++i) {
            unsigned product = 1;
            for (int j = 0; j < 4 && offsets[invariant][j] >= 0; ++j)
                product &= (phaseBits >> ((i + offsets[invariant][j]) % 23)) & 1u;
            parity ^= product;
        }
        word |= parity << invariant;
    }
    return word;
}

U192 coefficient(int tableSeed, int branch, int which, const U192 &ell) {
    U192 r;
    const uint64_t domain = 0x7ab1e0000000ull + uint64_t(tableSeed) * 0x100000ull
                          + uint64_t(branch) * 2 + uint64_t(which);
    for (int i = 0; i < 3; ++i) r.v[i] = R::eccPrfHost(domain, i);
    r.v[2] &= ~(1ull << 63);
    r = mod_reduce(r, ell);
    if (u192_is_zero(r)) r = u192_from(1);
    return r;
}

struct PointData {
    uint32_t scalar;
    uint8_t k;
    uint8_t eps;
    uint32_t necklace;
    uint64_t canonical;
};

struct GraphStats {
    uint64_t states = 0;
    uint64_t exceptional = 0;
    uint64_t cycles = 0;
    uint64_t fixed = 0;
    uint64_t cycles2 = 0;
    uint64_t cycles4 = 0;
    uint64_t cyclesLe8 = 0;
    uint64_t cycleNodesLe8 = 0;
    uint64_t basinLe8 = 0;
    uint64_t maxCycle = 0;
    uint64_t indegreeZero = 0;
    uint64_t indegreeOne = 0;
    uint64_t indegreeGe2 = 0;
    uint64_t indegreeCollisionPairs = 0;
    double meanFirstRepeat = 0;
    double medianFirstRepeat = 0;
    double p90FirstRepeat = 0;
    double meanFirstRepeatOverSqrtN = 0;
};

GraphStats analyse(const std::vector<uint32_t> &next) {
    GraphStats result;
    result.states = next.size();
    std::vector<uint32_t> indegree(next.size(), 0), removed;
    removed.reserve(next.size());
    for (uint32_t n : next) {
        if (n == BAD) ++result.exceptional;
        else ++indegree[n];
    }
    for (uint32_t degree : indegree) {
        result.indegreeZero += degree == 0;
        result.indegreeOne += degree == 1;
        result.indegreeGe2 += degree >= 2;
        if (degree >= 2)
            result.indegreeCollisionPairs += uint64_t(degree) * (degree - 1) / 2;
    }
    std::queue<uint32_t> pending;
    for (uint32_t i = 0; i < next.size(); ++i)
        if (indegree[i] == 0) pending.push(i);
    while (!pending.empty()) {
        const uint32_t here = pending.front();
        pending.pop();
        removed.push_back(here);
        const uint32_t there = next[here];
        if (there != BAD && --indegree[there] == 0) pending.push(there);
    }

    std::vector<uint32_t> cycleLength(next.size(), 0), tailLength(next.size(), 0);
    std::vector<uint8_t> seen(next.size(), 0);
    for (uint32_t i = 0; i < next.size(); ++i) {
        if (indegree[i] == 0 || seen[i]) continue;
        uint32_t length = 0, cur = i;
        do {
            seen[cur] = 1;
            ++length;
            cur = next[cur];
        } while (cur != i);
        ++result.cycles;
        result.maxCycle = std::max<uint64_t>(result.maxCycle, length);
        result.fixed += length == 1;
        result.cycles2 += length == 2;
        result.cycles4 += length == 4;
        result.cyclesLe8 += length <= 8;
        if (length <= 8) result.cycleNodesLe8 += length;
        cur = i;
        do {
            cycleLength[cur] = length;
            cur = next[cur];
        } while (cur != i);
    }
    for (auto it = removed.rbegin(); it != removed.rend(); ++it) {
        const uint32_t there = next[*it];
        if (there != BAD) {
            cycleLength[*it] = cycleLength[there];
            tailLength[*it] = tailLength[there] + 1;
        }
    }
    std::vector<uint32_t> firstRepeat;
    firstRepeat.reserve(next.size());
    long double total = 0;
    for (uint32_t i = 0; i < next.size(); ++i) {
        const uint32_t length = cycleLength[i];
        result.basinLe8 += length && length <= 8;
        if (!length) continue;
        const uint32_t repeat = tailLength[i] + length;
        firstRepeat.push_back(repeat);
        total += repeat;
    }
    if (!firstRepeat.empty()) {
        std::sort(firstRepeat.begin(), firstRepeat.end());
        result.meanFirstRepeat = double(total / firstRepeat.size());
        result.medianFirstRepeat = firstRepeat[firstRepeat.size() / 2];
        result.p90FirstRepeat = firstRepeat[(firstRepeat.size() * 9) / 10];
        result.meanFirstRepeatOverSqrtN = result.meanFirstRepeat / std::sqrt(double(next.size()));
    }
    return result;
}

void writeGraph(std::ostream &out, const char *prefix, const GraphStats &g) {
    out << ",\"" << prefix << "_states\":" << g.states
        << ",\"" << prefix << "_exceptional\":" << g.exceptional
        << ",\"" << prefix << "_cycles\":" << g.cycles
        << ",\"" << prefix << "_fixed\":" << g.fixed
        << ",\"" << prefix << "_cycles2\":" << g.cycles2
        << ",\"" << prefix << "_cycles4\":" << g.cycles4
        << ",\"" << prefix << "_cycles_le8\":" << g.cyclesLe8
        << ",\"" << prefix << "_cycle_nodes_le8\":" << g.cycleNodesLe8
        << ",\"" << prefix << "_basin_le8\":" << g.basinLe8
        << ",\"" << prefix << "_max_cycle\":" << g.maxCycle
        << ",\"" << prefix << "_indegree_zero\":" << g.indegreeZero
        << ",\"" << prefix << "_indegree_one\":" << g.indegreeOne
        << ",\"" << prefix << "_indegree_ge2\":" << g.indegreeGe2
        << ",\"" << prefix << "_indegree_collision_pairs\":" << g.indegreeCollisionPairs
        << ",\"" << prefix << "_mean_first_repeat\":" << g.meanFirstRepeat
        << ",\"" << prefix << "_median_first_repeat\":" << g.medianFirstRepeat
        << ",\"" << prefix << "_p90_first_repeat\":" << g.p90FirstRepeat
        << ",\"" << prefix << "_mean_first_repeat_over_sqrt_n\":"
        << g.meanFirstRepeatOverSqrtN;
}

} // namespace

int main(int argc, char **argv) {
    const char *outPath = argc > 1 ? argv[1] : "stateless-v4-f23.jsonl";
    const int seedCount = argc > 2 ? std::atoi(argv[2]) : 16;
    if (seedCount < 1 || seedCount > 16) {
        std::fprintf(stderr, "seed count must be in 1..16\n");
        return 2;
    }

    const U192 ellBig = u192_from_dec(eccF23::ELL_DEC);
    const uint64_t ell = small(ellBig);
    const uint64_t sigma = std::strtoull(eccF23::S_DEC, nullptr, 10);
    const Point basis = R::make(R::fromLimbs(eccF23::INSTANCE_PX[0]),
                                R::fromLimbs(eccF23::INSTANCE_PY[0]));
    if (!R::onCurve(basis) || !R::scalarMul(basis, ellBig).inf) {
        std::fprintf(stderr, "invalid generated basis\n");
        return 1;
    }

    std::array<uint64_t, 23> sigmaPower{};
    sigmaPower[0] = 1;
    for (int i = 1; i < 23; ++i) sigmaPower[i] = mulMod(sigmaPower[i - 1], sigma, ell);
    if (mulMod(sigmaPower[22], sigma, ell) != 1) {
        std::fprintf(stderr, "Frobenius scalar does not have order 23\n");
        return 1;
    }

    // Scalar orbit metadata. Code j means sigma^j; code 23+j means -sigma^j.
    std::vector<uint32_t> representative(ell, 0), representativeIndex(ell, BAD);
    std::vector<uint8_t> transform(ell, 0);
    std::vector<uint32_t> representatives;
    for (uint32_t a = 1; a < ell; ++a) {
        if (representative[a]) continue;
        const uint32_t index = representatives.size();
        representatives.push_back(a);
        uint64_t value = a;
        for (int j = 0; j < 23; ++j) {
            const uint32_t pos = uint32_t(value), neg = uint32_t(ell - value);
            if (representative[pos] || representative[neg]) {
                std::fprintf(stderr, "short or overlapping automorphism orbit\n");
                return 1;
            }
            representative[pos] = representative[neg] = a;
            representativeIndex[pos] = representativeIndex[neg] = index;
            transform[pos] = uint8_t(j);
            transform[neg] = uint8_t(23 + j);
            value = mulMod(value, sigma, ell);
        }
    }

    TableWalkConsts<23> constants;
    constants.build();
    if (!constants.consistent()) {
        std::fprintf(stderr, "invalid phase constants\n");
        return 1;
    }

    std::vector<PointData> pointData;
    pointData.reserve(representatives.size());
    uint64_t covarianceChecks = 0, covarianceFailures = 0;
    for (uint32_t scalar : representatives) {
        const Point point = R::scalarMul(basis, u192_from(scalar));
        const Elem x = R::nbCoords(point.x), y = R::nbCoords(point.y);
        const int weight = R::weight(point.x);
        const int k = constants.phase(x.v, weight);
        const int eps = constants.negationBit(x.v, y.v, k);
        const uint32_t nb = mix32(necklaceWord(x, constants)
                                  ^ (uint32_t(weight) * 0x9e3779b9u));
        const uint64_t cb = R::hashPoint(R::canonical(x));
        pointData.push_back(PointData{scalar, uint8_t(k), uint8_t(eps), nb, cb});
        for (int j = 0; j < 23; ++j) {
            const Point p = R::frob(point, j);
            const Elem px = R::nbCoords(p.x), py = R::nbCoords(p.y);
            const int pw = R::weight(p.x), pk = constants.phase(px.v, pw);
            const int pe = constants.negationBit(px.v, py.v, pk);
            covarianceFailures += mix32(necklaceWord(px, constants)
                                        ^ (uint32_t(pw) * 0x9e3779b9u)) != nb;
            covarianceFailures += R::hashPoint(R::canonical(px)) != cb;
            covarianceFailures += pk != (k + j) % 23;
            covarianceFailures += pe != eps;
            const Point n = R::neg(p);
            const Elem nx = R::nbCoords(n.x), ny = R::nbCoords(n.y);
            const int nw = R::weight(n.x), nk = constants.phase(nx.v, nw);
            const int ne = constants.negationBit(nx.v, ny.v, nk);
            covarianceFailures += mix32(necklaceWord(nx, constants)
                                        ^ (uint32_t(nw) * 0x9e3779b9u)) != nb;
            covarianceFailures += R::hashPoint(R::canonical(nx)) != cb;
            covarianceFailures += nk != pk;
            covarianceFailures += ne != (pe ^ 1);
            covarianceChecks += 8;
        }
    }
    if (covarianceFailures) {
        std::fprintf(stderr, "covariance failures: %llu / %llu\n",
                     (unsigned long long)covarianceFailures,
                     (unsigned long long)covarianceChecks);
        return 1;
    }

    std::ofstream out(outPath);
    if (!out) {
        std::perror(outPath);
        return 1;
    }
    out << "{\"kind\":\"meta\",\"degree\":23,\"subgroup_order\":" << ell
        << ",\"quotient_states\":" << representatives.size()
        << ",\"covariance_checks\":" << covarianceChecks
        << ",\"covariance_failures\":" << covarianceFailures << "}\n";

    constexpr int branches = 8;
    const uint64_t inverseTwo = (ell + 1) / 2;
    const double halvingRate = 28.953, additionRate = 20.1343255;
    const std::array<int, 6> halfThresholds = {0, 8, 12, 14, 15, 16};
    const std::array<const char *, 2> selectors = {"necklace7", "canonical_hash"};
    std::vector<uint32_t> quotientNext(representatives.size()), exactNext(ell - 1);

    // Same-size random functional graphs calibrate the finite first-repeat and
    // short-cycle statistics.  They are controls, not candidate maps.
    for (int seed = 0; seed < seedCount; ++seed) {
        for (uint32_t i = 0; i < quotientNext.size(); ++i) {
            const uint64_t random = R::eccPrfHost(0x4d4958454448414cull + uint64_t(seed), i);
            quotientNext[i] = uint32_t(random % quotientNext.size());
        }
        const GraphStats q = analyse(quotientNext);
        out << std::setprecision(17)
            << "{\"kind\":\"random_control\",\"degree\":23,\"seed\":" << seed;
        writeGraph(out, "quotient", q);
        out << "}\n";
    }

    for (const char *selector : selectors) {
        for (int threshold : halfThresholds) {
            for (int seed = 0; seed < seedCount; ++seed) {
                const uint64_t target = std::strtoull(eccF23::INSTANCE_K[seed], nullptr, 10);
                std::array<uint64_t, branches> tableScalar{};
                for (int h = 0; h < branches; ++h) {
                    const uint64_t a = small(coefficient(seed, h, 0, ellBig));
                    const uint64_t b = small(coefficient(seed, h, 1, ellBig));
                    tableScalar[h] = (a + mulMod(b, target, ell)) % ell;
                }

                uint64_t halfCount = 0, addCount = 0;
                std::array<uint64_t, branches> branchCount{};
                std::vector<uint32_t> representativeNext(representatives.size());
                for (uint32_t i = 0; i < pointData.size(); ++i) {
                    const PointData &p = pointData[i];
                    const uint64_t word = std::string(selector) == "necklace7"
                                            ? uint64_t(p.necklace) : p.canonical;
                    uint64_t next;
                    if (int(word & 15u) < threshold) {
                        ++halfCount;
                        next = mulMod(inverseTwo, p.scalar, ell);
                    } else {
                        ++addCount;
                        const int h = int((word >> 4) & (branches - 1));
                        ++branchCount[h];
                        uint64_t delta = mulMod(sigmaPower[p.k], tableScalar[h], ell);
                        if (p.eps) delta = delta ? ell - delta : 0;
                        next = (p.scalar + delta) % ell;
                    }
                    representativeNext[i] = uint32_t(next);
                    quotientNext[i] = next ? representativeIndex[next] : BAD;
                }

                for (uint32_t scalar = 1; scalar < ell; ++scalar) {
                    const uint32_t i = representativeIndex[scalar];
                    const uint8_t code = transform[scalar];
                    const int j = code >= 23 ? code - 23 : code;
                    uint64_t next = mulMod(sigmaPower[j], representativeNext[i], ell);
                    if (code >= 23 && next) next = ell - next;
                    exactNext[scalar - 1] = next ? uint32_t(next - 1) : BAD;
                }

                const GraphStats q = analyse(quotientNext), e = analyse(exactNext);
                const double halfFraction = double(halfCount) / pointData.size();
                const double optimisticRawRate = 1.0 /
                    (halfFraction / halvingRate + (1.0 - halfFraction) / additionRate);
                long double addSum2 = 0, addSum4 = 0;
                if (addCount) {
                    for (uint64_t count : branchCount) {
                        const long double p = (long double)count / addCount;
                        addSum2 += p * p;
                        addSum4 += p * p * p * p;
                    }
                }
                out << std::setprecision(17)
                    << "{\"kind\":\"mixed_map\",\"degree\":23,\"selector\":\"" << selector
                    << "\",\"branches\":" << branches << ",\"seed\":" << seed
                    << ",\"half_threshold_16\":" << threshold
                    << ",\"target_scalar\":" << target
                    << ",\"half_count\":" << halfCount << ",\"add_count\":" << addCount
                    << ",\"half_fraction\":" << halfFraction
                    << ",\"ideal_free_dispatch_raw_rate_billion\":" << optimisticRawRate
                    << ",\"add_sum_p2\":" << double(addSum2)
                    << ",\"add_sum_p4\":" << double(addSum4)
                    << ",\"add_sum_p2_uniform_ratio\":" << double(addSum2 * branches)
                    << ",\"add_sum_p4_uniform_ratio\":"
                    << double(addSum4 * branches * branches * branches);
                writeGraph(out, "quotient", q);
                writeGraph(out, "point", e);
                out << "}\n";
                out.flush();
            }
        }
    }
    std::printf("PASS: %zu quotient classes; %llu covariance checks; %d seeds; output %s\n",
                representatives.size(), (unsigned long long)covarianceChecks, seedCount, outPath);
    return 0;
}
