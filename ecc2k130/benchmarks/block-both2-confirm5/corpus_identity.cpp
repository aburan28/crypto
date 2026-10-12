#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Record = std::array<unsigned char, 32>;

struct Corpus {
    std::array<unsigned char, 16> header{};
    std::vector<Record> records;
    std::size_t duplicates = 0;
};

std::uint32_t le32(const unsigned char *p) {
    return std::uint32_t(p[0]) | (std::uint32_t(p[1]) << 8) |
           (std::uint32_t(p[2]) << 16) | (std::uint32_t(p[3]) << 24);
}

Corpus readCorpus(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    in.seekg(0, std::ios::end);
    const auto end = in.tellg();
    if (end < 16) throw std::runtime_error(path + ": shorter than v3 header");
    const std::size_t bytes = static_cast<std::size_t>(end);
    if ((bytes - 16) % 32) throw std::runtime_error(path + ": partial record");
    in.seekg(0);
    Corpus out;
    in.read(reinterpret_cast<char *>(out.header.data()), out.header.size());
    const unsigned char magic[8] = {'E','C','C','2','K','D','T','3'};
    if (std::memcmp(out.header.data(), magic, 8) != 0 ||
        le32(out.header.data() + 8) != 3 || le32(out.header.data() + 12) != 32)
        throw std::runtime_error(path + ": expected ECC2KDT3/v3/32 header");
    out.records.resize((bytes - 16) / 32);
    for (auto &record : out.records)
        in.read(reinterpret_cast<char *>(record.data()), record.size());
    if (!in) throw std::runtime_error(path + ": truncated payload");
    std::sort(out.records.begin(), out.records.end());
    out.duplicates = std::adjacent_find(out.records.begin(), out.records.end()) == out.records.end()
        ? 0 : std::size_t(-1);
    if (out.duplicates == std::size_t(-1)) {
        out.duplicates = 0;
        for (std::size_t i = 1; i < out.records.size(); ++i)
            if (out.records[i] == out.records[i - 1]) ++out.duplicates;
    }
    return out;
}

void writeSorted(const std::string &path, const Corpus &corpus) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);
    for (const auto &record : corpus.records)
        out.write(reinterpret_cast<const char *>(record.data()), record.size());
    if (!out) throw std::runtime_error("failed writing " + path);
}
}

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        const unsigned char header[16] = {'E','C','C','2','K','D','T','3',3,0,0,0,32,0,0,0};
        if (le32(header + 8) != 3 || le32(header + 12) != 32) return 1;
        std::cout << "PASS corpus_identity self-test\n";
        return 0;
    }
    if (argc != 5) {
        std::cerr << "usage: corpus_identity control.bin candidate.bin sorted-control sorted-candidate\n";
        return 2;
    }
    try {
        Corpus control = readCorpus(argv[1]);
        Corpus candidate = readCorpus(argv[2]);
        writeSorted(argv[3], control);
        writeSorted(argv[4], candidate);
        const bool sameHeader = control.header == candidate.header;
        const bool sameRecords = control.records == candidate.records;
        std::cout << "control\t" << control.records.size() << "\t"
                  << (control.records.size() - control.duplicates) << "\t"
                  << control.duplicates << "\n";
        std::cout << "candidate\t" << candidate.records.size() << "\t"
                  << (candidate.records.size() - candidate.duplicates) << "\t"
                  << candidate.duplicates << "\n";
        std::cout << "sameHeader\t" << (sameHeader ? 1 : 0) << "\n";
        std::cout << "sameSortedRecords\t" << (sameRecords ? 1 : 0) << "\n";
        return sameHeader && sameRecords ? 0 : 1;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
