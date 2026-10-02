#define main eccSparseStaticAuditEmbeddedMain
#include "../sparse-basis-compound/audit.cpp"
#undef main

#include <cmath>
#include <map>
#include <set>
#include <unordered_map>

namespace {

constexpr int kDegree = 131;
constexpr int kMapCount = 8;

struct DiagonalStats {
    std::array<std::array<uint32_t, 5>, 261> masks{};
    uint64_t matrixOnes = 0;
    int nonemptyDiagonals = 0;
    int wordTerms = 0;
    int zeroDiagonalTerms = 0;
    int exactShiftLop3Ops = 0;
};

DiagonalStats diagonals(const Columns &columns) {
    DiagonalStats out;
    for (int source = 0; source < kDegree; ++source)
        for (int destination = 0; destination < kDegree; ++destination) {
            if (!bit(columns[source], destination)) continue;
            ++out.matrixOnes;
            const int diagonal = destination - source + (kDegree - 1);
            out.masks[diagonal][destination / 32] |= 1u << (destination % 32);
        }
    for (int diagonal = 0; diagonal < 261; ++diagonal) {
        bool any = false;
        for (uint32_t mask : out.masks[diagonal]) {
            if (!mask) continue;
            any = true;
            ++out.wordTerms;
            if (diagonal == kDegree - 1) ++out.zeroDiagonalTerms;
        }
        if (any) ++out.nonemptyDiagonals;
    }
    out.exactShiftLop3Ops = out.wordTerms + (out.wordTerms - out.zeroDiagonalTerms);
    return out;
}

Columns xorColumns(const Columns &a, const Columns &b) {
    Columns out{};
    for (int i = 0; i < kDegree; ++i) out[i] = a[i] ^ b[i];
    return out;
}

Columns compose(const Columns &outer, const Columns &inner) {
    Columns out{};
    for (int i = 0; i < kDegree; ++i) out[i] = applyColumns(outer, inner[i]);
    return out;
}

bool sameColumns(const Columns &a, const Columns &b) {
    for (int i = 0; i < kDegree; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

std::vector<std::vector<int>> matrixRows(const Columns &columns) {
    std::vector<std::vector<int>> rows(kDegree);
    for (int row = 0; row < kDegree; ++row)
        for (int column = 0; column < kDegree; ++column)
            if (bit(columns[column], row)) rows[row].push_back(column);
    return rows;
}

struct PairHash {
    size_t operator()(uint64_t value) const {
        value ^= value >> 33;
        value *= 0xff51afd7ed558ccdull;
        value ^= value >> 33;
        return size_t(value);
    }
};

struct XorSynthesis {
    int outputs = 0;
    int temporaries = 0;
    int finalTerms = 0;
    int finalXors = 0;
    int totalXors = 0;
    int signals = 0;
};

XorSynthesis synthesizeXors(std::vector<std::vector<int>> rows) {
    int signals = kDegree, temporaries = 0;
    for (;;) {
        std::unordered_map<uint64_t, int, PairHash> counts;
        size_t pairs = 0;
        for (const auto &row : rows) pairs += row.size() * (row.size() - 1) / 2;
        counts.reserve(pairs);
        for (const auto &row : rows)
            for (size_t i = 0; i < row.size(); ++i)
                for (size_t j = i + 1; j < row.size(); ++j) {
                    uint32_t a = unsigned(row[i]), b = unsigned(row[j]);
                    if (a > b) std::swap(a, b);
                    ++counts[(uint64_t(a) << 32) | b];
                }
        int bestCount = 1;
        uint64_t bestPair = ~uint64_t(0);
        for (const auto &[pair, count] : counts)
            if (count > bestCount || (count == bestCount && pair < bestPair)) {
                bestCount = count;
                bestPair = pair;
            }
        if (bestCount <= 1) break;
        const int a = int(bestPair >> 32), b = int(uint32_t(bestPair));
        const int temporary = signals++;
        for (auto &row : rows) {
            auto ia = std::find(row.begin(), row.end(), a);
            auto ib = std::find(row.begin(), row.end(), b);
            if (ia == row.end() || ib == row.end()) continue;
            if (ia > ib) std::swap(ia, ib);
            row.erase(ib);
            row.erase(ia);
            row.push_back(temporary);
        }
        ++temporaries;
    }
    int finalTerms = 0, finalXors = 0;
    for (const auto &row : rows) {
        finalTerms += int(row.size());
        if (!row.empty()) finalXors += int(row.size()) - 1;
    }
    return XorSynthesis{int(rows.size()), temporaries, finalTerms, finalXors,
                        temporaries + finalXors, signals};
}

long double binomialProbability(int weight) {
    long double coefficient = 1;
    for (int i = 1; i <= weight; ++i)
        coefficient = coefficient * (kDegree - weight + i) / i;
    return std::ldexp(coefficient, -kDegree);
}

void printDiagonal(std::ostream &out, const char *name, const DiagonalStats &stats,
                   bool comma) {
    out << "    \"" << name << "\": {\"matrixOnes\": " << stats.matrixOnes
        << ", \"nonemptyDiagonals\": " << stats.nonemptyDiagonals
        << ", \"wordTerms\": " << stats.wordTerms
        << ", \"zeroDiagonalTerms\": " << stats.zeroDiagonalTerms
        << ", \"optimisticFreeShiftLop3Ops\": " << stats.wordTerms
        << ", \"exactShiftPlusLop3Ops\": " << stats.exactShiftLop3Ops << "}"
        << (comma ? "," : "") << "\n";
}

void printSynthesis(std::ostream &out, const char *name, const XorSynthesis &result,
                    bool comma) {
    out << "    \"" << name << "\": {\"outputs\": " << result.outputs
        << ", \"temporaries\": " << result.temporaries
        << ", \"signals\": " << result.signals
        << ", \"finalTerms\": " << result.finalTerms
        << ", \"finalXors\": " << result.finalXors
        << ", \"totalXors\": " << result.totalXors << "}"
        << (comma ? "," : "") << "\n";
}

}  // namespace

int main() {
    const std::vector<int> betaLower = {0, 2, 3, 64, 66, 67, 96, 98, 99,
                                        112, 114, 115, 120, 122, 123, 124, 128, 130};
    Field sparse({0, 2, 3, 8});
    const E betaRoot{{0x962fc4e3ddc388ebull, 0x0a16693fefe59e60ull, 3}};
    if (!sparse.rabinPrimeDegree() || !zero(evalPolynomial(sparse, betaRoot, betaLower)))
        fail("field/root preflight");

    Columns betaToSparse{};
    betaToSparse[0] = basis(0);
    for (int i = 1; i < kDegree; ++i)
        betaToSparse[i] = sparse.mul(betaToSparse[i - 1], betaRoot);
    const Columns sparseToBeta = inverse(betaToSparse);
    if (rank(betaToSparse) != kDegree || rank(sparseToBeta) != kDegree)
        fail("basis conversion rank");

    Columns sparseToNormal{};
    std::array<Columns, 11> l{};
    for (int power = 1; power <= 10; ++power)
        for (int i = 0; i < kDegree; ++i) {
            const E input = basis(i);
            l[power][i] = input ^ sparse.frob(input, power);
        }
    for (int i = 0; i < kDegree; ++i)
        sparseToNormal[i] = repositoryFromPolynomial(applyColumns(sparseToBeta, basis(i)));
    if (rank(sparseToNormal) != kDegree) fail("normal conversion rank");
    for (int power = 3; power <= 10; ++power)
        if (rank(l[power]) != 130) fail("L_j rank");

    int compositionChecks = 0;
    const Columns frobenius = [&] {
        Columns out{};
        for (int i = 0; i < kDegree; ++i) out[i] = sparse.frob(basis(i), 1);
        return out;
    }();
    for (int a = 1; a <= 9; ++a)
        for (int b = 1; a + b <= 10; ++b) {
            const Columns rhs = xorColumns(l[a], compose([&] {
                Columns fa{};
                for (int i = 0; i < kDegree; ++i) fa[i] = sparse.frob(basis(i), a);
                return fa;
            }(), l[b]));
            if (!sameColumns(l[a + b], rhs)) fail("L_(a+b) identity");
            ++compositionChecks;
        }
    for (int a = 1; 2 * a <= 10; ++a) {
        if (!sameColumns(l[2 * a], compose(l[a], l[a]))) fail("L_2a identity");
        ++compositionChecks;
    }
    for (int power = 1; power < 10; ++power) {
        if (!sameColumns(l[power + 1], xorColumns(compose(frobenius, l[power]), l[1])))
            fail("L_(j+1) identity");
        ++compositionChecks;
    }

    uint64_t state = 0x414c554d41505331ull;
    int denseChecks = 0;
    for (int test = 0; test < 4096; ++test) {
        const E input = randomElement(state);
        const E beta = applyColumns(sparseToBeta, input);
        if (applyColumns(betaToSparse, beta) != input ||
            applyColumns(sparseToNormal, input) != repositoryFromPolynomial(beta))
            fail("dense conversion check");
        for (int power = 3; power <= 10; ++power) {
            const E want = input ^ sparse.frob(input, power);
            const E repository =
                applyColumns(betaToSparse, repositoryL(beta, power));
            if (applyColumns(l[power], input) != want || want != repository)
                fail("dense L_j check");
            ++denseChecks;
        }
    }

    const DiagonalStats normalDiag = diagonals(sparseToNormal);
    const DiagonalStats toBetaDiag = diagonals(sparseToBeta);
    const DiagonalStats fromBetaDiag = diagonals(betaToSparse);
    std::array<DiagonalStats, kMapCount> lDiag{};
    for (int map = 0; map < kMapCount; ++map) lDiag[map] = diagonals(l[map + 3]);

    int unionTerms = 0, unionZeroTerms = 0, totalMemberships = 0;
    int uniqueMasks = 0, identicalAcrossAll = 0;
    for (int diagonal = 0; diagonal < 261; ++diagonal)
        for (int word = 0; word < 5; ++word) {
            std::set<uint32_t> masks;
            bool any = false, all = true, same = true;
            uint32_t first = lDiag[0].masks[diagonal][word];
            for (int map = 0; map < kMapCount; ++map) {
                const uint32_t mask = lDiag[map].masks[diagonal][word];
                if (mask) {
                    any = true;
                    masks.insert(mask);
                    ++totalMemberships;
                } else {
                    all = false;
                }
                if (mask != first) same = false;
            }
            if (!any) continue;
            ++unionTerms;
            if (diagonal == kDegree - 1) ++unionZeroTerms;
            uniqueMasks += int(masks.size());
            if (all && same) ++identicalAcrossAll;
        }

    std::map<std::array<uint64_t, 3>, int> outputRows;
    int duplicateRows = 0;
    for (int map = 0; map < kMapCount; ++map)
        for (int row = 0; row < kDegree; ++row) {
            std::array<uint64_t, 3> bits{};
            for (int column = 0; column < kDegree; ++column)
                if (bit(l[map + 3][column], row)) bits[column / 64] |= 1ull << (column % 64);
            if (++outputRows[bits] > 1) ++duplicateRows;
        }

    const XorSynthesis normalXor = synthesizeXors(matrixRows(sparseToNormal));
    const XorSynthesis l3Xor = synthesizeXors(matrixRows(l[3]));
    auto normalPlusL3Rows = matrixRows(sparseToNormal);
    const auto l3Rows = matrixRows(l[3]);
    normalPlusL3Rows.insert(normalPlusL3Rows.end(), l3Rows.begin(), l3Rows.end());
    const XorSynthesis normalPlusL3Xor = synthesizeXors(normalPlusL3Rows);

    std::array<long double, kMapCount> selectorProbability{};
    for (int weightValue = 0; weightValue <= kDegree; ++weightValue)
        selectorProbability[(weightValue >> 1) & 7] += binomialProbability(weightValue);
    std::array<long double, kMapCount> warpPresence{};
    long double expectedDistinct = 0;
    for (int map = 0; map < kMapCount; ++map) {
        warpPresence[map] = 1 - std::pow(1 - selectorProbability[map], 32);
        expectedDistinct += warpPresence[map];
    }

    const int cheapestMap = int(std::min_element(
        lDiag.begin(), lDiag.end(),
        [](const DiagonalStats &a, const DiagonalStats &b) {
            return a.wordTerms < b.wordTerms;
        }) - lDiag.begin());
    const double conversionOptimistic =
        double(toBetaDiag.wordTerms + fromBetaDiag.wordTerms) / 16.0;
    const double oracleMapOnly = normalDiag.wordTerms + 2.0 * lDiag[cheapestMap].wordTerms +
                                 conversionOptimistic;
    const double oracleExactMapOnly = normalDiag.exactShiftLop3Ops +
        2.0 * lDiag[cheapestMap].exactShiftLop3Ops +
        double(toBetaDiag.exactShiftLop3Ops + fromBetaDiag.exactShiftLop3Ops) / 16.0;
    double meanTerms = 0, divergentExact = 0;
    for (int map = 0; map < kMapCount; ++map) {
        meanTerms += double(selectorProbability[map]) * lDiag[map].wordTerms;
        divergentExact += double(warpPresence[map]) * lDiag[map].exactShiftLop3Ops;
    }
    const double compactedMeanMapOnly =
        normalDiag.wordTerms + 2.0 * meanTerms + conversionOptimistic;
    const double divergentMapOnly = normalDiag.exactShiftLop3Ops +
        2.0 * divergentExact +
        double(toBetaDiag.exactShiftLop3Ops + fromBetaDiag.exactShiftLop3Ops) / 16.0;
    const int sharedShiftOps = unionTerms - unionZeroTerms;
    const double allOutputPerCoordinate = sharedShiftOps + totalMemberships + 35;
    const double allOutputMapOnly = normalDiag.exactShiftLop3Ops +
        2.0 * allOutputPerCoordinate +
        double(toBetaDiag.exactShiftLop3Ops + fromBetaDiag.exactShiftLop3Ops) / 16.0;

    constexpr double sparseReducerOps = 476.625;
    constexpr double aluSquareSpreadOps = 75.0;
    constexpr double referenceSubledger = 2276.25;
    const double oracleCandidateSubledger =
        oracleMapOnly + sparseReducerOps + aluSquareSpreadOps;

    constexpr double sms = 188, clockGhz = 2.430, aluLanes = 64;
    constexpr double fusedRate = 15.436677, targetRate = 26.0;
    const double fusedAluBudget = sms * clockGhz * aluLanes / fusedRate;
    const double targetAluBudget = sms * clockGhz * aluLanes / targetRate;
    const double oracleMapAluCeiling = sms * clockGhz * aluLanes / oracleMapOnly;
    const double oracleSubledgerAluCeiling = sms * clockGhz * aluLanes / oracleCandidateSubledger;

    constexpr int squareSourceOps = 5 * 15 + 82;
    const double compactedMeanSquares = 2.0 * 6.5;
    const double compactedSquareMapOps = normalDiag.exactShiftLop3Ops +
        compactedMeanSquares * squareSourceOps + conversionOptimistic;
    const double branchlessSquareMapOps = normalDiag.exactShiftLop3Ops +
        20.0 * squareSourceOps + conversionOptimistic + 70.0;
    const double compactedClmadSquares = compactedMeanSquares * 5.0;
    const double compactedClmadCeiling = sms * clockGhz * 2.0 / compactedClmadSquares;

    const bool survivesTarget =
        oracleMapOnly <= targetAluBudget && oracleCandidateSubledger < referenceSubledger;

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "{\n"
              << "  \"schema\": \"ecc2k130-sparse-alu-synthesis-v1\",\n"
              << "  \"valid\": true,\n"
              << "  \"sourceParent\": \"b3d0095f6fb74c180c71780ead8ab692b4831cc1\",\n"
              << "  \"checks\": {\"basisVectors\": 131, \"denseVectors\": 4096, \"denseLjComparisons\": "
              << denseChecks << ", \"compositionIdentities\": " << compositionChecks << "},\n"
              << "  \"aluBudgetsPerUpdate\": {\"fused15_436677\": " << fusedAluBudget
              << ", \"target26\": " << targetAluBudget << "},\n"
              << "  \"diagonalCircuits\": {\n";
    printDiagonal(std::cout, "sparseToNormal", normalDiag, true);
    printDiagonal(std::cout, "sparseToBeta", toBetaDiag, true);
    printDiagonal(std::cout, "betaToSparse", fromBetaDiag, true);
    for (int map = 0; map < kMapCount; ++map) {
        const std::string name = "L" + std::to_string(map + 3);
        printDiagonal(std::cout, name.c_str(), lDiag[map], map + 1 != kMapCount);
    }
    std::cout << "  },\n"
              << "  \"jointCommonStructure\": {\"unionWordTerms\": " << unionTerms
              << ", \"unionZeroDiagonalTerms\": " << unionZeroTerms
              << ", \"totalMapMemberships\": " << totalMemberships
              << ", \"uniqueNonzeroMasksAcrossTerms\": " << uniqueMasks
              << ", \"identicalMasksAcrossAllEight\": " << identicalAcrossAll
              << ", \"uniqueOutputRows\": " << outputRows.size()
              << ", \"duplicateOutputRows\": " << duplicateRows << "},\n"
              << "  \"xorCse\": {\n";
    printSynthesis(std::cout, "sparseToNormal", normalXor, true);
    printSynthesis(std::cout, "L3", l3Xor, true);
    printSynthesis(std::cout, "sparseToNormalPlusL3", normalPlusL3Xor, false);
    std::cout << "  },\n"
              << "  \"dynamicSelector\": {\"expectedDistinctBranchesPerWarp\": "
              << double(expectedDistinct) << ", \"selectorProbabilities\": [";
    for (int map = 0; map < kMapCount; ++map) {
        if (map) std::cout << ',';
        std::cout << double(selectorProbability[map]);
    }
    std::cout << "], \"warpPresenceProbabilities\": [";
    for (int map = 0; map < kMapCount; ++map) {
        if (map) std::cout << ',';
        std::cout << double(warpPresence[map]);
    }
    std::cout << "]},\n"
              << "  \"mapOnlySourceOpsPerUpdate\": {\n"
              << "    \"oracleCheapestFreeShiftFreeSelection\": " << oracleMapOnly << ",\n"
              << "    \"oracleCheapestExactDiagonal\": " << oracleExactMapOnly << ",\n"
              << "    \"perfectCompactionWeightedMeanFreeShift\": " << compactedMeanMapOnly << ",\n"
              << "    \"warpDivergentExactDiagonal\": " << divergentMapOnly << ",\n"
              << "    \"allOutputsSharedShiftsAndBranchlessSelect\": " << allOutputMapOnly << ",\n"
              << "    \"greedyXorNormalOnly\": " << normalXor.totalXors << ",\n"
              << "    \"greedyXorNormalPlusL3XPlusL3Y\": "
              << normalPlusL3Xor.totalXors + l3Xor.totalXors << "\n"
              << "  },\n"
              << "  \"frobeniusComposition\": {\"sparseAluSquareSourceOps\": "
              << squareSourceOps << ", \"perfectCompactionMeanSquaresPerUpdate\": "
              << compactedMeanSquares << ", \"perfectCompactionMapOps\": "
              << compactedSquareMapOps << ", \"branchlessTenSquareMapOps\": "
              << branchlessSquareMapOps << ", \"perfectCompactionClmadSquaresPerUpdate\": "
              << compactedClmadSquares << ", \"perfectCompactionClmadOnlyCeilingBillionPerSecond\": "
              << compactedClmadCeiling << "},\n"
              << "  \"completeB16ComparableSubledger\": {\n"
              << "    \"referenceSelectorAndReductionSourceOps\": " << referenceSubledger << ",\n"
              << "    \"oracleCandidateMaps\": " << oracleMapOnly << ",\n"
              << "    \"sparseReducers\": " << sparseReducerOps << ",\n"
              << "    \"aluSquareSpread\": " << aluSquareSpreadOps << ",\n"
              << "    \"oracleCandidateTotal\": " << oracleCandidateSubledger << ",\n"
              << "    \"candidateOverReference\": " << oracleCandidateSubledger / referenceSubledger
              << "\n  },\n"
              << "  \"ceilingsBillionUpdatesPerSecond\": {\"oracleMapOnlyAlu\": "
              << oracleMapAluCeiling << ", \"oracleComparableSubledgerAlu\": "
              << oracleSubledgerAluCeiling << "},\n"
              << "  \"decision\": {\"survivesStaticGate\": "
              << (survivesTarget ? "true" : "false")
              << ", \"classification\": \"NEGATIVE_STATIC_EVIDENCE\", \"reason\": \"Even the impossible oracle that assigns L3 to every update, makes every shift and dynamic selection free, and ignores all other walk work spends more ALU on maps than the entire 26 B/s budget. Adding the frozen sparse reducer and ALU-square spread is also worse than the reference comparable source subledger.\"},\n"
              << "  \"scope\": \"Bounded native static circuit feasibility only; no CUDA compile, GPU run, sparse walk, search, solver or key recovery.\"\n"
              << "}\n";
    return survivesTarget ? 2 : 0;
}
