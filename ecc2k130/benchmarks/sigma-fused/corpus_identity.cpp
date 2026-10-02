#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace {
struct Header { char magic[8]; uint32_t version, recordBytes; };

bool read(const std::string &path, std::vector<std::string> *records) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { std::fprintf(stderr, "%s: missing\n", path.c_str()); return false; }
    std::vector<char> bytes((std::istreambuf_iterator<char>(in)), {});
    if (bytes.size() < sizeof(Header)) {
        std::fprintf(stderr, "%s: truncated v2 corpus\n", path.c_str()); return false;
    }
    Header h{};
    std::memcpy(&h, bytes.data(), sizeof(h));
    if (std::memcmp(h.magic, "ECC2KDP2", 8) || h.version != 2 || h.recordBytes != 72) {
        std::fprintf(stderr, "%s: expected ECC2KDP2 version 2, 72-byte records\n", path.c_str());
        return false;
    }
    const size_t payload = bytes.size() - sizeof(Header);
    if (payload % h.recordBytes) {
        std::fprintf(stderr, "%s: partial record\n", path.c_str()); return false;
    }
    for (size_t off = sizeof(Header); off < bytes.size(); off += h.recordBytes)
        records->emplace_back(bytes.data() + off, h.recordBytes);
    std::sort(records->begin(), records->end());
    return true;
}
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: corpus_identity CORPUS CORPUS [CORPUS...]\n");
        return 2;
    }
    std::vector<std::string> reference;
    if (!read(argv[1], &reference)) return 1;
    std::printf("%s: %zu sorted v2 records (reference)\n", argv[1], reference.size());
    for (int i = 2; i < argc; ++i) {
        std::vector<std::string> candidate;
        if (!read(argv[i], &candidate)) return 1;
        const bool equal = candidate == reference;
        std::printf("%s: %zu sorted v2 records (%s)\n", argv[i], candidate.size(),
                    equal ? "IDENTICAL" : "DIFFERS");
        if (!equal) return 1;
    }
    return 0;
}
