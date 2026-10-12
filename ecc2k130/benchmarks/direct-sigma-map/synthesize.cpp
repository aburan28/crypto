#include "../../include/packed131.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using eccPacked131::P131;

namespace {

constexpr int kBits = 131;
constexpr int kWords = 5;
constexpr int kMaps = 8;

using Columns = std::array<P131, kBits>;
using Maps = std::array<Columns, kMaps>;

bool same(const P131 &a, const P131 &b) {
    for (int i = 0; i < kWords; ++i)
        if (a.v[i] != b.v[i]) return false;
    return true;
}

bool canonical(const P131 &a) { return (a.v[4] & ~7u) == 0; }

bool bit(const P131 &a, int i) {
    return ((a.v[i / 32] >> (i % 32)) & 1u) != 0;
}

void setBit(P131 &a, int i) { a.v[i / 32] |= 1u << (i % 32); }

P131 basis(int i) {
    P131 a{{0, 0, 0, 0, 0}};
    setBit(a, i);
    return a;
}

P131 composed(P131 p, int j) {
    const P131 x = eccPacked131::fromPolynomial131(p);
    return eccPacked131::toPolynomial131(
        eccPacked131::add131(x, eccPacked131::sigma131(x, j)));
}

P131 applyColumns(const Columns &columns, const P131 &input) {
    P131 out{{0, 0, 0, 0, 0}};
    for (int i = 0; i < kBits; ++i) {
        if (!bit(input, i)) continue;
        for (int w = 0; w < kWords; ++w) out.v[w] ^= columns[i].v[w];
    }
    return out;
}

int rank(const Columns &columns) {
    std::array<std::array<uint64_t, 3>, kBits> rows{};
    for (int column = 0; column < kBits; ++column)
        for (int row = 0; row < kBits; ++row)
            if (bit(columns[column], row)) rows[row][column / 64] ^= 1ull << (column % 64);
    int pivot = 0;
    for (int column = 0; column < kBits && pivot < kBits; ++column) {
        int chosen = pivot;
        while (chosen < kBits && ((rows[chosen][column / 64] >> (column % 64)) & 1u) == 0) ++chosen;
        if (chosen == kBits) continue;
        std::swap(rows[pivot], rows[chosen]);
        for (int r = 0; r < kBits; ++r) {
            if (r == pivot || ((rows[r][column / 64] >> (column % 64)) & 1u) == 0) continue;
            for (int w = 0; w < 3; ++w) rows[r][w] ^= rows[pivot][w];
        }
        ++pivot;
    }
    return pivot;
}

P131 shift131(const P131 &a, int displacement) {
    P131 out{{0, 0, 0, 0, 0}};
    for (int source = 0; source < kBits; ++source) {
        const int destination = source + displacement;
        if (destination >= 0 && destination < kBits && bit(a, source)) setBit(out, destination);
    }
    return out;
}

struct DiagonalFamily {
    std::array<P131, 2 * kBits - 1> masks{};
    int diagonals = 0;
    int wordTerms = 0;
    int matrixOnes = 0;
};

DiagonalFamily makeDiagonals(const Columns &columns) {
    DiagonalFamily family;
    for (int source = 0; source < kBits; ++source) {
        for (int destination = 0; destination < kBits; ++destination) {
            if (!bit(columns[source], destination)) continue;
            ++family.matrixOnes;
            setBit(family.masks[destination - source + kBits - 1], destination);
        }
    }
    for (const P131 &mask : family.masks) {
        bool any = false;
        for (uint32_t word : mask.v) {
            if (word) {
                any = true;
                ++family.wordTerms;
            }
        }
        family.diagonals += any;
    }
    return family;
}

P131 applyDiagonals(const DiagonalFamily &family, const P131 &input) {
    P131 out{{0, 0, 0, 0, 0}};
    for (int i = 0; i < 2 * kBits - 1; ++i) {
        const P131 shifted = shift131(input, i - (kBits - 1));
        for (int w = 0; w < kWords; ++w) out.v[w] ^= shifted.v[w] & family.masks[i].v[w];
    }
    return out;
}

struct ChunkFamily {
    int bits = 0;
    int chunks = 0;
    int entries = 0;
    std::vector<P131> table;
};

struct Half5Family {
    std::array<size_t, 27> offsets{};
    std::array<int, 27> entries{};
    std::vector<P131> lowTable;
    std::array<P131, 27> highColumns{};
};

ChunkFamily makeChunks(const Columns &columns, int chunkBits) {
    ChunkFamily family;
    family.bits = chunkBits;
    family.chunks = (kBits + chunkBits - 1) / chunkBits;
    family.entries = 1 << chunkBits;
    family.table.resize(size_t(family.chunks) * family.entries);
    for (int chunk = 0; chunk < family.chunks; ++chunk) {
        for (int value = 0; value < family.entries; ++value) {
            P131 &out = family.table[size_t(chunk) * family.entries + value];
            for (int b = 0; b < chunkBits; ++b) {
                const int source = chunk * chunkBits + b;
                if (source >= kBits || ((value >> b) & 1) == 0) continue;
                for (int w = 0; w < kWords; ++w) out.v[w] ^= columns[source].v[w];
            }
        }
    }
    return family;
}

unsigned extractChunk(const P131 &input, int start, int width) {
    unsigned out = 0;
    for (int i = 0; i < width && start + i < kBits; ++i)
        out |= unsigned(bit(input, start + i)) << i;
    return out;
}

P131 applyChunks(const ChunkFamily &family, const P131 &input) {
    P131 out{{0, 0, 0, 0, 0}};
    for (int chunk = 0; chunk < family.chunks; ++chunk) {
        const unsigned value = extractChunk(input, chunk * family.bits, family.bits);
        const P131 &part = family.table[size_t(chunk) * family.entries + value];
        for (int w = 0; w < kWords; ++w) out.v[w] ^= part.v[w];
    }
    return out;
}

Half5Family makeHalf5(const Columns &columns) {
    Half5Family family;
    for (int chunk = 0; chunk < 27; ++chunk) {
        const int start = chunk * 5;
        const int entries = start + 4 < kBits ? 16 : 2;
        family.offsets[chunk] = family.lowTable.size();
        family.entries[chunk] = entries;
        for (int value = 0; value < entries; ++value) {
            P131 out{{0, 0, 0, 0, 0}};
            for (int b = 0; b < 4; ++b) {
                const int source = start + b;
                if (source >= kBits || ((value >> b) & 1) == 0) continue;
                for (int w = 0; w < kWords; ++w) out.v[w] ^= columns[source].v[w];
            }
            family.lowTable.push_back(out);
        }
        if (start + 4 < kBits) family.highColumns[chunk] = columns[start + 4];
    }
    return family;
}

P131 applyHalf5(const Half5Family &family, const P131 &input) {
    P131 out{{0, 0, 0, 0, 0}};
    for (int chunk = 0; chunk < 27; ++chunk) {
        const unsigned value = extractChunk(input, chunk * 5, 5);
        const unsigned low = value & unsigned(family.entries[chunk] - 1);
        const P131 &part = family.lowTable[family.offsets[chunk] + low];
        for (int w = 0; w < kWords; ++w) out.v[w] ^= part.v[w];
        if ((value & 16u) != 0)
            for (int w = 0; w < kWords; ++w) out.v[w] ^= family.highColumns[chunk].v[w];
    }
    return out;
}

uint32_t nextWord(uint32_t &state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

void fail(const std::string &message) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
}

std::string hexWord(uint32_t word) {
    std::ostringstream out;
    out << "0x" << std::hex << std::setw(8) << std::setfill('0') << word << "u";
    return out.str();
}

void emitTable3(const std::string &path, const std::array<ChunkFamily, kMaps> &families) {
    std::ofstream out(path);
    if (!out) fail("cannot open generated header " + path);
    const int chunks = families[0].chunks;
    const int entries = families[0].entries;
    out << "// Generated by benchmarks/direct-sigma-map/synthesize.cpp.\n";
    out << "#pragma once\n";
    out << "namespace eccDirectSigma131 {\n";
    out << "alignas(32) static const uint32_t table3[" << kMaps << "][" << chunks
        << "][" << entries << "][" << kWords << "] = {\n";
    for (int map = 0; map < kMaps; ++map) {
        out << " {\n";
        for (int chunk = 0; chunk < chunks; ++chunk) {
            out << "  {\n";
            for (int value = 0; value < entries; ++value) {
                const P131 &p = families[map].table[size_t(chunk) * entries + value];
                out << "   {";
                for (int w = 0; w < kWords; ++w) {
                    if (w) out << ',';
                    out << hexWord(p.v[w]);
                }
                out << "},\n";
            }
            out << "  },\n";
        }
        out << " },\n";
    }
    out << "};\n";
    out << "inline eccPacked131::P131 apply3(eccPacked131::P131 a, int index) {\n";
    out << " eccPacked131::P131 out{{0,0,0,0,0}};\n";
    out << " for (int chunk=0; chunk<" << chunks << "; ++chunk) {\n";
    out << "  const int start=chunk*3, word=start>>5, offset=start&31;\n";
    out << "  uint32_t value=a.v[word]>>offset;\n";
    out << "  if(offset>29 && word<4) value|=a.v[word+1]<<(32-offset);\n";
    out << "  value&=7u;\n";
    out << "  for(int w=0;w<5;++w) out.v[w]^=table3[index][chunk][value][w];\n";
    out << " }\n";
    out << " out.v[4]&=7u; return out;\n";
    out << "}\n";
    for (int map = 0; map < kMaps; ++map) {
        out << "inline eccPacked131::P131 apply3_j" << (map + 3)
            << "(eccPacked131::P131 a) {\n";
        out << " eccPacked131::P131 out{{0,0,0,0,0}};\n";
        for (int chunk = 0; chunk < chunks; ++chunk) {
            const int start = chunk * 3;
            const int word = start / 32;
            const int offset = start % 32;
            out << " { uint32_t value=a.v[" << word << "]>>" << offset << ";";
            if (offset > 29 && word < 4)
                out << " value|=a.v[" << (word + 1) << "]<<" << (32 - offset) << ";";
            out << " value&=7u;\n";
            for (int w = 0; w < kWords; ++w)
                out << "  out.v[" << w << "]^=table3[" << map << "][" << chunk
                    << "][value][" << w << "];\n";
            out << " }\n";
        }
        out << " out.v[4]&=7u; return out;\n";
        out << "}\n";
    }
    out << "}\n";
}

void emitHalf5(const std::string &path, const std::array<Half5Family, kMaps> &families) {
    std::ofstream out(path);
    if (!out) fail("cannot open generated header " + path);
    constexpr int entries = 26 * 16 + 2;
    out << "// Generated by benchmarks/direct-sigma-map/synthesize.cpp.\n";
    out << "#pragma once\n";
    out << "namespace eccDirectSigmaHalf5 {\n";
    out << "alignas(32) static const uint32_t low128[" << kMaps << "][" << entries << "][4] = {\n";
    for (int map = 0; map < kMaps; ++map) {
        out << " {\n";
        for (const P131 &p : families[map].lowTable) {
            out << "  {";
            for (int w = 0; w < 4; ++w) {
                if (w) out << ',';
                out << hexWord(p.v[w]);
            }
            out << "},\n";
        }
        out << " },\n";
    }
    out << "};\n";
    out << "alignas(32) static const uint8_t top3[" << kMaps << "][" << entries << "] = {\n";
    for (int map = 0; map < kMaps; ++map) {
        out << " {";
        for (size_t i = 0; i < families[map].lowTable.size(); ++i) {
            if (i) out << ',';
            out << unsigned(families[map].lowTable[i].v[4]);
        }
        out << "},\n";
    }
    out << "};\n";
    for (int map = 0; map < kMaps; ++map) {
        out << "inline eccPacked131::P131 apply5_j" << (map + 3)
            << "(eccPacked131::P131 a) {\n";
        out << " eccPacked131::P131 out{{0,0,0,0,0}};\n";
        for (int chunk = 0; chunk < 27; ++chunk) {
            const int start = chunk * 5;
            const int word = start / 32;
            const int offset = start % 32;
            const size_t tableOffset = families[map].offsets[chunk];
            const int entryMask = families[map].entries[chunk] - 1;
            out << " { uint32_t value=a.v[" << word << "]>>" << offset << ";";
            if (offset > 27 && word < 4)
                out << " value|=a.v[" << (word + 1) << "]<<" << (32 - offset) << ";";
            out << " const uint32_t low=value&" << entryMask << "u;\n";
            for (int w = 0; w < 4; ++w)
                out << "  out.v[" << w << "]^=low128[" << map << "]["
                    << tableOffset << "+low][" << w << "];\n";
            out << "  out.v[4]^=top3[" << map << "][" << tableOffset << "+low];\n";
            if (start + 4 < kBits) {
                out << "  if(value&16u) {\n";
                for (int w = 0; w < kWords; ++w)
                    if (families[map].highColumns[chunk].v[w] != 0)
                        out << "   out.v[" << w << "]^="
                            << hexWord(families[map].highColumns[chunk].v[w]) << ";\n";
                out << "  }\n";
            }
            out << " }\n";
        }
        out << " out.v[4]&=7u; return out;\n";
        out << "}\n";
    }
    out << "}\n";
}

}  // namespace

int main(int argc, char **argv) {
    const std::string output = argc > 1 ? argv[1] : "direct_sigma_table3.generated.h";
    const std::string output5 = argc > 2 ? argv[2] : "direct_sigma_half5.generated.h";
    Maps maps{};
    std::array<DiagonalFamily, kMaps> diagonals;
    std::array<ChunkFamily, kMaps> table3;
    std::array<Half5Family, kMaps> half5;
    uint32_t rng = 0x1315a17u;

    std::cout << "source_operation_model=v1\n";
    std::cout << "from_reduced=" << ECC_PACKED_FROM_REDUCED << "\n";
    for (int map = 0; map < kMaps; ++map) {
        const int j = map + 3;
        for (int source = 0; source < kBits; ++source) {
            maps[map][source] = composed(basis(source), j);
            if (!canonical(maps[map][source])) fail("non-canonical oracle basis output");
        }
        if (rank(maps[map]) != 130) fail("unexpected rank for j=" + std::to_string(j));
        diagonals[map] = makeDiagonals(maps[map]);
        table3[map] = makeChunks(maps[map], 3);
        half5[map] = makeHalf5(maps[map]);

        for (int source = 0; source < kBits; ++source) {
            const P131 input = basis(source);
            const P131 want = maps[map][source];
            if (!same(applyColumns(maps[map], input), want)) fail("column oracle basis mismatch");
            if (!same(applyDiagonals(diagonals[map], input), want)) fail("diagonal basis mismatch");
            if (!same(applyChunks(table3[map], input), want)) fail("table3 basis mismatch");
            if (!same(applyHalf5(half5[map], input), want)) fail("half5 basis mismatch");
        }
        for (int test = 0; test < 256; ++test) {
            P131 input{{nextWord(rng), nextWord(rng), nextWord(rng), nextWord(rng), nextWord(rng) & 7u}};
            if (test == 0) input = P131{{0, 0, 0, 0, 0}};
            if (test == 1) input = P131{{~0u, ~0u, ~0u, ~0u, 7u}};
            const P131 want = composed(input, j);
            if (!same(applyColumns(maps[map], input), want)) fail("dense column mismatch");
            if (!same(applyDiagonals(diagonals[map], input), want)) fail("dense diagonal mismatch");
            if (!same(applyChunks(table3[map], input), want)) fail("dense table3 mismatch");
            if (!same(applyHalf5(half5[map], input), want)) fail("dense half5 mismatch");
        }

        std::cout << "j=" << j << " rank=130 matrix_ones=" << diagonals[map].matrixOnes
                  << " nonempty_diagonals=" << diagonals[map].diagonals
                  << " diagonal_word_terms=" << diagonals[map].wordTerms << '\n';
    }

    for (int bits = 1; bits <= 8; ++bits) {
        const ChunkFamily family = makeChunks(maps[0], bits);
        const size_t bytes = size_t(kMaps) * family.chunks * family.entries * kWords * sizeof(uint32_t);
        const int loads = family.chunks * kWords;
        const int xors = loads - kWords;
        std::cout << "chunk_bits=" << bits << " chunks=" << family.chunks
                  << " all_maps_bytes=" << bytes << " loads_per_coordinate=" << loads
                  << " accumulator_xors=" << xors << '\n';
    }

    emitTable3(output, table3);
    emitHalf5(output5, half5);
    std::cout << "half5_all_maps_bytes=56848 loads_per_coordinate=135"
              << " conditional_high_xors_per_coordinate_at_warp_divergence=130\n";
    std::cout << "PASS: 1048 basis-vector map cases and 2048 deterministic dense/edge cases"
              << " for oracle, diagonal, table3, and half5 families\n";
    std::cout << "generated=" << output << '\n';
    std::cout << "generated=" << output5 << '\n';
    return 0;
}
