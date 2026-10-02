#include "../../include/packed131.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using eccPacked131::P131;

namespace {

constexpr int kBits = 131;
constexpr int kWords32 = 5;
constexpr int kChunks = 44;
constexpr int kChunkEntries = 8;
constexpr int kMaps = 8;
constexpr int kDenseCases = 4096;
constexpr int kRawDenseCases = 4096;

struct E {
    std::array<uint64_t, 3> w{};
};

bool operator==(const E &a, const E &b) { return a.w == b.w; }
bool operator!=(const E &a, const E &b) { return !(a == b); }

E operator^(E a, const E &b) {
    for (int i = 0; i < 3; ++i) a.w[i] ^= b.w[i];
    return a;
}

E &operator^=(E &a, const E &b) {
    for (int i = 0; i < 3; ++i) a.w[i] ^= b.w[i];
    return a;
}

bool bit(const E &a, int i) { return ((a.w[i / 64] >> (i % 64)) & 1u) != 0; }

void toggle(E &a, int i) { a.w[i / 64] ^= uint64_t(1) << (i % 64); }

E basis(int i) {
    E a;
    toggle(a, i);
    return a;
}

bool zero(const E &a) { return a.w[0] == 0 && a.w[1] == 0 && a.w[2] == 0; }

unsigned weight(const E &a) {
    return unsigned(__builtin_popcountll(a.w[0]) + __builtin_popcountll(a.w[1]) +
                    __builtin_popcountll(a.w[2]));
}

std::string hex(const E &a) {
    std::ostringstream out;
    out << "0x" << std::hex << a.w[2] << std::setfill('0') << std::setw(16) << a.w[1]
        << std::setw(16) << a.w[0];
    return out.str();
}

struct Wide {
    std::array<uint64_t, 5> w{};
};

bool wideBit(const Wide &a, int i) { return ((a.w[i / 64] >> (i % 64)) & 1u) != 0; }

void wideToggle(Wide &a, int i) { a.w[i / 64] ^= uint64_t(1) << (i % 64); }

class Field {
  public:
    explicit Field(std::vector<int> lower) : lower_(std::move(lower)) {}

    E reduce(Wide a) const {
        for (int i = 260; i >= kBits; --i) {
            if (!wideBit(a, i)) continue;
            wideToggle(a, i);
            const int shift = i - kBits;
            for (int t : lower_) wideToggle(a, shift + t);
        }
        E out{{a.w[0], a.w[1], a.w[2] & 7u}};
        return out;
    }

    Wide rawProduct(const E &a, const E &b) const {
        Wide out;
        for (int i = 0; i < kBits; ++i) {
            if (!bit(b, i)) continue;
            for (int j = 0; j < kBits; ++j)
                if (bit(a, j)) wideToggle(out, i + j);
        }
        return out;
    }

    E mul(const E &a, const E &b) const { return reduce(rawProduct(a, b)); }

    E sqr(const E &a) const {
        Wide out;
        for (int i = 0; i < kBits; ++i)
            if (bit(a, i)) wideToggle(out, 2 * i);
        return reduce(out);
    }

    E frob(E a, int j) const {
        for (int i = 0; i < j; ++i) a = sqr(a);
        return a;
    }

    E powSmall(E a, unsigned exponent) const {
        E out = basis(0);
        while (exponent) {
            if (exponent & 1u) out = mul(out, a);
            exponent >>= 1;
            if (exponent) a = sqr(a);
        }
        return out;
    }

    bool rabinPrimeDegree() const {
        E x = basis(1), y = x;
        for (int i = 0; i < kBits; ++i) y = sqr(y);
        if (y != x) return false;
        // For prime degree 131, the only proper Rabin divisor is one.  The
        // remaining gcd is gcd(f, x^2+x), so f(0) and f(1) must both be one.
        const bool atOne = ((1u + unsigned(lower_.size())) & 1u) != 0;
        const bool atZero = std::find(lower_.begin(), lower_.end(), 0) != lower_.end();
        return atZero && atOne;
    }

    const std::vector<int> &lower() const { return lower_; }

  private:
    std::vector<int> lower_;
};

std::array<bool, 131> quotientBits() {
    // Binary long division of 2^131-1 by 263, high bit first.  The type-II
    // construction here has ord_263(2)=131, so the roots are in the base
    // field; no quadratic extension is involved.
    std::array<bool, 131> q{};
    unsigned remainder = 0;
    for (int bit = 130; bit >= 0; --bit) {
        remainder = 2 * remainder + 1;
        if (remainder >= 263) {
            q[bit] = true;
            remainder -= 263;
        }
    }
    if (remainder != 0) throw std::runtime_error("263 does not divide 2^131-1");
    return q;
}

E powBits(const Field &field, const E &x, const std::array<bool, 131> &bits) {
    E out = basis(0);
    for (int bit = 130; bit >= 0; --bit) {
        out = field.sqr(out);
        if (bits[bit]) out = field.mul(out, x);
    }
    return out;
}

uint64_t splitmix(uint64_t &state) {
    uint64_t z = (state += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

E randomElement(uint64_t &state) {
    E out{{splitmix(state), splitmix(state), splitmix(state) & 7u}};
    return out;
}

E evalPolynomial(const Field &field, const E &x, const std::vector<int> &lower) {
    E acc = basis(0);  // monic coefficient of x^131
    for (int i = 130; i >= 0; --i) {
        acc = field.mul(acc, x);
        if (std::find(lower.begin(), lower.end(), i) != lower.end()) acc ^= basis(0);
    }
    return acc;
}

E deterministicBetaRoot(const Field &field, const std::vector<int> &oldLower,
                        unsigned *attempts) {
    const auto exponent = quotientBits();
    uint64_t state = 0x5350415253453133ull;
    const E one = basis(0);
    for (unsigned attempt = 1; attempt <= 1024; ++attempt) {
        E candidate = randomElement(state);
        if (zero(candidate)) candidate = basis(1);
        const E zeta = powBits(field, candidate, exponent);
        if (zeta == one || field.powSmall(zeta, 263) != one) continue;
        const E root = zeta ^ field.powSmall(zeta, 262);
        if (zero(root) || !zero(evalPolynomial(field, root, oldLower))) continue;
        *attempts = attempt;
        return root;
    }
    throw std::runtime_error("failed to construct beta-basis root");
}

using Columns = std::array<E, kBits>;

E applyColumns(const Columns &columns, const E &input) {
    E out;
    for (int i = 0; i < kBits; ++i)
        if (bit(input, i)) out ^= columns[i];
    return out;
}

int rank(const Columns &columns) {
    std::array<E, kBits> rows{};
    for (int c = 0; c < kBits; ++c)
        for (int r = 0; r < kBits; ++r)
            if (bit(columns[c], r)) toggle(rows[r], c);
    int pivot = 0;
    for (int c = 0; c < kBits; ++c) {
        int chosen = pivot;
        while (chosen < kBits && !bit(rows[chosen], c)) ++chosen;
        if (chosen == kBits) continue;
        std::swap(rows[pivot], rows[chosen]);
        for (int r = 0; r < kBits; ++r)
            if (r != pivot && bit(rows[r], c)) rows[r] ^= rows[pivot];
        ++pivot;
    }
    return pivot;
}

Columns inverse(const Columns &columns) {
    struct Row {
        E left, right;
    };
    std::array<Row, kBits> rows{};
    for (int r = 0; r < kBits; ++r) {
        toggle(rows[r].right, r);
        for (int c = 0; c < kBits; ++c)
            if (bit(columns[c], r)) toggle(rows[r].left, c);
    }
    for (int c = 0; c < kBits; ++c) {
        int chosen = c;
        while (chosen < kBits && !bit(rows[chosen].left, c)) ++chosen;
        if (chosen == kBits) throw std::runtime_error("singular conversion matrix");
        std::swap(rows[c], rows[chosen]);
        for (int r = 0; r < kBits; ++r) {
            if (r == c || !bit(rows[r].left, c)) continue;
            rows[r].left ^= rows[c].left;
            rows[r].right ^= rows[c].right;
        }
    }
    Columns out{};
    for (int c = 0; c < kBits; ++c)
        for (int r = 0; r < kBits; ++r)
            if (bit(rows[r].right, c)) toggle(out[c], r);
    return out;
}

uint64_t matrixOnes(const Columns &columns) {
    uint64_t total = 0;
    for (const E &column : columns) total += weight(column);
    return total;
}

struct ChunkTable {
    std::array<E, kChunks * kChunkEntries> entries{};
};

ChunkTable makeTable(const Columns &columns) {
    ChunkTable table;
    for (int chunk = 0; chunk < kChunks; ++chunk) {
        for (int value = 0; value < kChunkEntries; ++value) {
            E &out = table.entries[chunk * kChunkEntries + value];
            for (int b = 0; b < 3; ++b) {
                const int source = 3 * chunk + b;
                if (source < kBits && ((value >> b) & 1)) out ^= columns[source];
            }
        }
    }
    return table;
}

unsigned chunkValue(const E &a, int chunk) {
    const int start = 3 * chunk;
    unsigned value = 0;
    for (int b = 0; b < 3 && start + b < kBits; ++b)
        value |= unsigned(bit(a, start + b)) << b;
    return value;
}

E applyTable(const ChunkTable &table, const E &input) {
    E out;
    for (int chunk = 0; chunk < kChunks; ++chunk)
        out ^= table.entries[chunk * kChunkEntries + chunkValue(input, chunk)];
    return out;
}

P131 packed(const E &a) {
    return P131{{uint32_t(a.w[0]), uint32_t(a.w[0] >> 32), uint32_t(a.w[1]),
                 uint32_t(a.w[1] >> 32), uint32_t(a.w[2])}};
}

E unpacked(const P131 &a) {
    return E{{uint64_t(a.v[0]) | (uint64_t(a.v[1]) << 32),
              uint64_t(a.v[2]) | (uint64_t(a.v[3]) << 32), uint64_t(a.v[4] & 7u)}};
}

E repositoryFromPolynomial(const E &a) {
    return unpacked(eccPacked131::fromPolynomial131(packed(a)));
}

E repositoryL(const E &betaPolynomial, int j) {
    const P131 normal = eccPacked131::fromPolynomial131(packed(betaPolynomial));
    return unpacked(eccPacked131::toPolynomial131(
        eccPacked131::add131(normal, eccPacked131::sigma131(normal, j))));
}

std::array<uint32_t, 9> words32(const Wide &a) {
    std::array<uint32_t, 9> out{};
    for (int i = 0; i < 9; ++i) out[i] = uint32_t(a.w[i / 2] >> (32 * (i & 1)));
    return out;
}

uint32_t shiftedLeft(const std::array<uint32_t, 5> &a, int word, int shift) {
    const uint32_t low = a[word] << shift;
    const uint32_t high = word ? a[word - 1] >> (32 - shift) : 0u;
    return low | high;
}

E sparseCircuit(const std::array<uint32_t, 9> &h, int middle) {
    std::array<uint32_t, 5> d{};
    for (int i = 0; i < 4; ++i) d[i] = (h[i + 4] >> 3) | (h[i + 5] << 29);
    d[4] = (h[8] >> 3) & 3u;
    std::array<uint32_t, 5> folded{};
    for (int i = 0; i < 5; ++i)
        folded[i] = d[i] ^ shiftedLeft(d, i, 2) ^ shiftedLeft(d, i, middle) ^
                    shiftedLeft(d, i, 8);
    const uint32_t overflow = folded[4] >> 3;
    std::array<uint32_t, 5> out{};
    for (int i = 0; i < 5; ++i) out[i] = h[i] ^ folded[i];
    out[0] ^= overflow ^ (overflow << 2) ^ (overflow << middle) ^ (overflow << 8);
    out[4] &= 7u;
    return E{{uint64_t(out[0]) | (uint64_t(out[1]) << 32),
              uint64_t(out[2]) | (uint64_t(out[3]) << 32), out[4]}};
}

void fail(const std::string &message) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
}

struct CandidateResult {
    std::string name;
    std::string modulus;
    std::string betaRoot;
    unsigned rootAttempts = 0;
    uint64_t betaToSparseOnes = 0;
    uint64_t sparseToBetaOnes = 0;
    uint64_t sparseToNormalOnes = 0;
    std::array<uint64_t, kMaps> lOnes{};
    int reductionBasisChecks = 0;
    int reductionDenseChecks = 0;
    int conversionBasisChecks = 0;
    int conversionDenseChecks = 0;
    int mapBasisChecks = 0;
    int mapDenseChecks = 0;
    int productDenseChecks = 0;
};

CandidateResult auditCandidate(const std::string &name, const std::string &modulus,
                               const std::vector<int> &lower, const Field &oldField,
                               const std::vector<int> &oldLower) {
    Field field(lower);
    if (!field.rabinPrimeDegree()) fail(name + ": Rabin irreducibility failed");

    CandidateResult result;
    result.name = name;
    result.modulus = modulus;
    const E betaRoot = deterministicBetaRoot(field, oldLower, &result.rootAttempts);
    result.betaRoot = hex(betaRoot);

    Columns betaToSparse{};
    betaToSparse[0] = basis(0);
    for (int i = 1; i < kBits; ++i) betaToSparse[i] = field.mul(betaToSparse[i - 1], betaRoot);
    if (rank(betaToSparse) != kBits) fail(name + ": beta-to-sparse matrix rank");
    const Columns sparseToBeta = inverse(betaToSparse);
    if (rank(sparseToBeta) != kBits) fail(name + ": sparse-to-beta matrix rank");
    result.betaToSparseOnes = matrixOnes(betaToSparse);
    result.sparseToBetaOnes = matrixOnes(sparseToBeta);

    Columns sparseToNormal{};
    std::array<Columns, kMaps> linear{};
    for (int i = 0; i < kBits; ++i) {
        const E s = basis(i);
        const E beta = applyColumns(sparseToBeta, s);
        sparseToNormal[i] = repositoryFromPolynomial(beta);
        for (int m = 0; m < kMaps; ++m) {
            const int j = m + 3;
            const E direct = s ^ field.frob(s, j);
            const E oracleSparse = applyColumns(betaToSparse, repositoryL(beta, j));
            if (direct != oracleSparse) fail(name + ": sparse Frobenius basis mismatch");
            linear[m][i] = direct;
        }
    }
    if (rank(sparseToNormal) != kBits) fail(name + ": sparse-to-normal matrix rank");
    result.sparseToNormalOnes = matrixOnes(sparseToNormal);
    for (int m = 0; m < kMaps; ++m) {
        if (rank(linear[m]) != 130) fail(name + ": L_j rank is not 130");
        result.lOnes[m] = matrixOnes(linear[m]);
    }

    const ChunkTable normalTable = makeTable(sparseToNormal);
    const ChunkTable toBetaTable = makeTable(sparseToBeta);
    const ChunkTable fromBetaTable = makeTable(betaToSparse);
    std::array<ChunkTable, kMaps> lTables{};
    for (int m = 0; m < kMaps; ++m) lTables[m] = makeTable(linear[m]);

    for (int i = 0; i < kBits; ++i) {
        const E s = basis(i);
        const E beta = applyTable(toBetaTable, s);
        if (applyTable(fromBetaTable, beta) != s || applyColumns(sparseToBeta, s) != beta)
            fail(name + ": conversion table basis mismatch");
        ++result.conversionBasisChecks;
        if (applyTable(normalTable, s) != repositoryFromPolynomial(beta))
            fail(name + ": normal table basis mismatch");
        ++result.mapBasisChecks;
        for (int m = 0; m < kMaps; ++m) {
            if (applyTable(lTables[m], s) != linear[m][i])
                fail(name + ": L_j table basis mismatch");
            ++result.mapBasisChecks;
        }
    }

    // Every raw product basis vector, including degrees 131..260, exercises
    // the two-fold sparse reducer independently of field multiplication.
    for (int degree = 0; degree <= 260; ++degree) {
        Wide raw;
        wideToggle(raw, degree);
        const auto h = words32(raw);
        const int middle = lower[2] == 3 ? 3 : 5;
        if (sparseCircuit(h, middle) != field.reduce(raw))
            fail(name + ": sparse reducer basis mismatch");
        ++result.reductionBasisChecks;
    }

    uint64_t state = 0x434f4d504f554e44ull ^ uint64_t(lower[2]);
    for (int test = 0; test < kRawDenseCases; ++test) {
        Wide raw{{splitmix(state), splitmix(state), splitmix(state), splitmix(state),
                  splitmix(state) & 31u}};
        const auto h = words32(raw);
        const int middle = lower[2] == 3 ? 3 : 5;
        if (sparseCircuit(h, middle) != field.reduce(raw))
            fail(name + ": sparse reducer dense mismatch");
        ++result.reductionDenseChecks;
    }

    for (int test = 0; test < kDenseCases; ++test) {
        const E a = randomElement(state), b = randomElement(state);
        const E sa = applyColumns(betaToSparse, a), sb = applyColumns(betaToSparse, b);
        if (applyColumns(sparseToBeta, sa) != a || applyColumns(sparseToBeta, sb) != b)
            fail(name + ": dense conversion round trip mismatch");
        if (applyTable(toBetaTable, sa) != a || applyTable(fromBetaTable, a) != sa)
            fail(name + ": dense conversion table mismatch");
        result.conversionDenseChecks += 2;

        const E sparseProduct = field.mul(sa, sb);
        const E betaProduct = oldField.mul(a, b);
        if (sparseProduct != applyColumns(betaToSparse, betaProduct))
            fail(name + ": dense product isomorphism mismatch");
        const E repositoryProduct = unpacked(eccPacked131::mulPolynomial131(packed(a), packed(b)));
        if (repositoryProduct != betaProduct)
            fail(name + ": repository beta product mismatch");
        ++result.productDenseChecks;

        if (applyTable(normalTable, sa) != repositoryFromPolynomial(a))
            fail(name + ": dense normal table mismatch");
        ++result.mapDenseChecks;
        for (int m = 0; m < kMaps; ++m) {
            const int j = m + 3;
            const E got = applyTable(lTables[m], sa);
            const E want = applyColumns(betaToSparse, repositoryL(a, j));
            if (got != want || got != (sa ^ field.frob(sa, j)))
                fail(name + ": dense L_j mismatch");
            ++result.mapDenseChecks;
        }
    }

    return result;
}

void printCandidate(const CandidateResult &r, bool comma) {
    std::cout << "    {\n"
              << "      \"name\": \"" << r.name << "\",\n"
              << "      \"modulus\": \"" << r.modulus << "\",\n"
              << "      \"rabinIrreducible\": true,\n"
              << "      \"betaRoot\": \"" << r.betaRoot << "\",\n"
              << "      \"order263ConstructionAttempts\": " << r.rootAttempts << ",\n"
              << "      \"matrixRanks\": {\"betaToSparse\": 131, \"sparseToBeta\": 131, \"sparseToNormal\": 131, \"Lj\": [130,130,130,130,130,130,130,130]},\n"
              << "      \"matrixOnes\": {\"betaToSparse\": " << r.betaToSparseOnes
              << ", \"sparseToBeta\": " << r.sparseToBetaOnes
              << ", \"sparseToNormal\": " << r.sparseToNormalOnes << ", \"Lj\": [";
    for (int i = 0; i < kMaps; ++i) {
        if (i) std::cout << ',';
        std::cout << r.lOnes[i];
    }
    std::cout << "]},\n"
              << "      \"checks\": {\"reductionBasis\": " << r.reductionBasisChecks
              << ", \"reductionDense\": " << r.reductionDenseChecks
              << ", \"conversionBasis\": " << r.conversionBasisChecks
              << ", \"conversionDense\": " << r.conversionDenseChecks
              << ", \"mapBasis\": " << r.mapBasisChecks
              << ", \"mapDense\": " << r.mapDenseChecks
              << ", \"denseProducts\": " << r.productDenseChecks << "}\n"
              << "    }" << (comma ? "," : "") << "\n";
}

}  // namespace

int main() {
    try {
        const std::vector<int> oldLower = {0, 2, 3, 64, 66, 67, 96, 98, 99, 112,
                                           114, 115, 120, 122, 123, 124, 128, 130};
        const Field oldField(oldLower);
        if (!oldField.rabinPrimeDegree()) fail("repository beta modulus failed Rabin check");

        // Bind the generic old-field reducer to the actual repository reducer
        // on every raw-product direction and a deterministic dense panel.
        int oldBasisChecks = 0, oldDenseChecks = 0;
        for (int degree = 0; degree <= 260; ++degree) {
            Wide raw;
            wideToggle(raw, degree);
            const auto h = words32(raw);
            if (unpacked(eccPacked131::reducePolynomial131(h.data())) != oldField.reduce(raw))
                fail("repository reducer basis mismatch");
            ++oldBasisChecks;
        }
        uint64_t rawState = 0x4245544152454455ull;
        for (int test = 0; test < kRawDenseCases; ++test) {
            Wide raw{{splitmix(rawState), splitmix(rawState), splitmix(rawState),
                      splitmix(rawState), splitmix(rawState) & 31u}};
            const auto h = words32(raw);
            if (unpacked(eccPacked131::reducePolynomial131(h.data())) != oldField.reduce(raw))
                fail("repository reducer dense mismatch");
            ++oldDenseChecks;
        }

        // lower vectors are kept in ascending order; element two is the
        // candidate's distinguishing middle exponent used by the circuit.
        const CandidateResult first =
            auditCandidate("pentanomial-8-3-2", "0x80000000000000000000000000000010d",
                           {0, 2, 3, 8}, oldField, oldLower);
        const CandidateResult second =
            auditCandidate("pentanomial-8-5-2", "0x800000000000000000000000000000125",
                           {0, 2, 5, 8}, oldField, oldLower);

        constexpr int mapBytes = kChunks * kChunkEntries * kWords32 * 4;
        constexpr int runtimeMaps = 11;
        constexpr int tableBytes = runtimeMaps * mapBytes;
        constexpr int driverSharedBytes = 1024;
        constexpr int sharedBytes = tableBytes + driverSharedBytes;
        constexpr int ldsPerMap = 88;
        constexpr int scalarWordsPerMap = 220;
        constexpr int aluPerMap = 317;
        constexpr double mapsPerUpdate = 3.0 + 2.0 / 16.0;
        constexpr double ldsPerUpdate = ldsPerMap * mapsPerUpdate;
        constexpr double scalarWordsPerUpdate = scalarWordsPerMap * mapsPerUpdate;
        constexpr double mapAluPerUpdate = aluPerMap * mapsPerUpdate;

        constexpr int currentReducerShifts = 86, currentReducerOrs = 39;
        constexpr int currentReducerXors = 53, currentReducerAnds = 2;
        constexpr int sparseReducerShifts = 40, sparseReducerOrs = 16;
        constexpr int sparseReducerXors = 24, sparseReducerAnds = 2;
        constexpr int currentReducerOps = currentReducerShifts + currentReducerOrs +
                                          currentReducerXors + currentReducerAnds;
        constexpr int sparseReducerOps = sparseReducerShifts + sparseReducerOrs +
                                         sparseReducerXors + sparseReducerAnds;

        constexpr int batch = 16;
        constexpr int forwardProducts = 30;
        constexpr int inverseProducts = 8;
        constexpr int reversePrefixProducts = 31;
        constexpr int lambdaProducts = 16;
        constexpr int finalProducts = 16;
        constexpr int totalProducts = forwardProducts + inverseProducts + reversePrefixProducts +
                                      lambdaProducts + finalProducts;
        constexpr int sparseProductReductions = totalProducts - inverseProducts;
        constexpr int sparseSquares = batch;
        constexpr double sparseReductionsPerUpdate =
            double(sparseProductReductions + sparseSquares) / batch;
        constexpr double betaReductionsPerUpdate = double(inverseProducts) / batch;
        constexpr double clmadPerUpdate = double(totalProducts * 6 + 5 * 4) / batch;

        constexpr double sms = 188.0, maxClockGhz = 2.430;
        constexpr double randomLdsLanes = 9.2, aluLanes = 64.0, clmadLanes = 2.0;
        constexpr double referenceBps = 15.436677;
        const double ldsCeiling = sms * maxClockGhz * randomLdsLanes / ldsPerUpdate;
        const double clmadCeiling = sms * maxClockGhz * clmadLanes / clmadPerUpdate;
        const double tableAluCeiling = sms * maxClockGhz * aluLanes / mapAluPerUpdate;
        const double admission = 1.05 * referenceBps;
        const bool exactChecks = true;
        const bool sharedFits = sharedBytes <= 101376;
        constexpr int fusedRegisters = 126, tableAddedLiveWords = 11;
        const bool registerEstimateFits = fusedRegisters + tableAddedLiveWords <= 255;
        const bool rooflinePass = ldsCeiling >= admission && clmadCeiling >= admission &&
                                  tableAluCeiling >= admission;
        const bool gpuPrototype = exactChecks && sharedFits && registerEstimateFits && rooflinePass;

        const double currentReductionSourceOps =
            double((totalProducts + sparseSquares) * currentReducerOps) / batch;
        const double candidateReductionSourceOps =
            sparseReductionsPerUpdate * sparseReducerOps + betaReductionsPerUpdate * currentReducerOps;
        constexpr double currentSelectionSourceOps = 1230.0;
        const double candidateSelectionSourceOps = mapAluPerUpdate;

        std::cout << std::fixed << std::setprecision(9);
        std::cout << "{\n"
                  << "  \"schema\": \"ecc2k130-sparse-basis-compound-static-v1\",\n"
                  << "  \"valid\": true,\n"
                  << "  \"sourceParent\": \"5434e208953527d2f37e14bde9f6c470542922cf\",\n"
                  << "  \"referenceBillionUpdatesPerSecond\": " << referenceBps << ",\n"
                  << "  \"repositoryBetaModulus\": {\"hex\": \"0xd1d0d000d0000000d000000000000000d\", \"weight\": 19, \"rabinIrreducible\": true, \"reductionBasisChecks\": "
                  << oldBasisChecks << ", \"reductionDenseChecks\": " << oldDenseChecks << "},\n"
                  << "  \"candidates\": [\n";
        printCandidate(first, true);
        printCandidate(second, false);
        std::cout << "  ],\n"
                  << "  \"tableRoute\": {\n"
                  << "    \"chunkBits\": 3, \"chunks\": 44, \"entriesPerChunk\": 8, \"wordsPerEntry\": 5,\n"
                  << "    \"bytesPerMap\": " << mapBytes << ", \"hotMaps\": 9, \"batchBoundaryMaps\": 2, \"runtimeMaps\": " << runtimeMaps << ",\n"
                  << "    \"tableBytes\": " << tableBytes << ", \"driverSharedBytes\": " << driverSharedBytes << ", \"totalSharedBytes\": " << sharedBytes << ", \"optinLimitBytes\": 101376, \"fitsOneBlock\": " << (sharedFits ? "true" : "false") << ",\n"
                  << "    \"perMap\": {\"lds128\": 44, \"lds32\": 44, \"ldsInstructions\": " << ldsPerMap << ", \"scalarWordEquivalents\": " << scalarWordsPerMap << ", \"dataAlu\": " << aluPerMap << ", \"minimumLiveWords\": " << tableAddedLiveWords << "},\n"
                  << "    \"mapsPerUpdateAtBatch16\": " << mapsPerUpdate << ", \"ldsInstructionsPerUpdate\": " << ldsPerUpdate << ", \"scalarWordEquivalentsPerUpdate\": " << scalarWordsPerUpdate << ", \"dataAluPerUpdate\": " << mapAluPerUpdate << "\n"
                  << "  },\n"
                  << "  \"reducerCircuits\": {\n"
                  << "    \"currentDirect\": {\"shifts\": " << currentReducerShifts << ", \"ors\": " << currentReducerOrs << ", \"xors\": " << currentReducerXors << ", \"ands\": " << currentReducerAnds << ", \"total32BitSourceOps\": " << currentReducerOps << "},\n"
                  << "    \"sparseTwoFold\": {\"shifts\": " << sparseReducerShifts << ", \"ors\": " << sparseReducerOrs << ", \"xors\": " << sparseReducerXors << ", \"ands\": " << sparseReducerAnds << ", \"total32BitSourceOps\": " << sparseReducerOps << ", \"ratioToCurrent\": " << double(sparseReducerOps) / currentReducerOps << "}\n"
                  << "  },\n"
                  << "  \"batch16Arithmetic\": {\n"
                  << "    \"productsPerBatch\": {\"forward\": " << forwardProducts << ", \"inverse\": " << inverseProducts << ", \"reversePrefix\": " << reversePrefixProducts << ", \"lambda\": " << lambdaProducts << ", \"finalXY\": " << finalProducts << ", \"total\": " << totalProducts << "},\n"
                  << "    \"sparseReductionsPerUpdateIncludingLambdaSquare\": " << sparseReductionsPerUpdate << ", \"betaReductionsPerUpdate\": " << betaReductionsPerUpdate << ", \"clmadPerUpdate\": " << clmadPerUpdate << ",\n"
                  << "    \"currentReductionSourceOpsPerUpdate\": " << currentReductionSourceOps << ", \"candidateReductionSourceOpsPerUpdate\": " << candidateReductionSourceOps << ",\n"
                  << "    \"currentSelectionAndLjSourceOpsPerUpdate\": " << currentSelectionSourceOps << ", \"candidateSelectionAndLjPlusBatchMapsSourceOpsPerUpdate\": " << candidateSelectionSourceOps << ",\n"
                  << "    \"sourceOpsSavedAcrossThoseStagesPerUpdate\": " << (currentReductionSourceOps + currentSelectionSourceOps - candidateReductionSourceOps - candidateSelectionSourceOps) << "\n"
                  << "  },\n"
                  << "  \"storage\": {\"fusedBaselineRegisters\": " << fusedRegisters << ", \"tableEvaluatorAddedLiveWords\": " << tableAddedLiveWords << ", \"naiveRegisterUpperEstimate\": " << fusedRegisters + tableAddedLiveWords << ", \"registerLimit\": 255, \"reducerTemporaryWordsCurrent\": 20, \"reducerTemporaryWordsSparse\": 10},\n"
                  << "  \"optimisticRooflineBillionPerSecond\": {\n"
                  << "    \"assumptions\": {\"sms\": 188, \"maxClockGhz\": " << maxClockGhz << ", \"randomSharedLanesPerSmClock\": " << randomLdsLanes << ", \"aluLanesPerSmClock\": " << aluLanes << ", \"clmadLanesPerSmClock\": " << clmadLanes << "},\n"
                  << "    \"lookupOnly\": " << ldsCeiling << ", \"tableAluOnly\": " << tableAluCeiling << ", \"clmadOnly\": " << clmadCeiling << ", \"admissionThreshold\": " << admission << ", \"lookupOnlyRatioToReference\": " << ldsCeiling / referenceBps << "\n"
                  << "  },\n"
                  << "  \"decision\": {\"boundedGpuPrototypeJustified\": " << (gpuPrototype ? "true" : "false") << ", \"classification\": \"NEGATIVE_STATIC_EVIDENCE\", \"reason\": \"Even with all multiplication, reduction, state traffic and control made free, the admitted three-bit shared-table route has a 15.283375 B/s lookup-only ceiling at the most favourable recorded clock, below the confirmed 15.436677 B/s reference and the frozen 5% admission margin.\"},\n"
                  << "  \"scope\": \"Exact native field/isomorphism/linear-map evidence and a conservative static pipe ledger; no CUDA compile, GPU run, search, solver, collision recovery or key recovery.\"\n"
                  << "}\n";
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
