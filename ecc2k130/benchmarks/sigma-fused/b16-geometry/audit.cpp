#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

std::string read(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("missing " + path);
    return std::string(std::istreambuf_iterator<char>(input), {});
}

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

void contains(const std::string &text, const std::string &marker) {
    require(text.find(marker) != std::string::npos, "missing marker: " + marker);
}

struct Arm {
    const char *name;
    int blockThreads;
    int minBlocks;
    int residentBlocks;
};

constexpr Arm kArms[] = {{"t256", 256, 2, 2}, {"t512", 512, 1, 1}};

void commonMarkers(const std::string &log, const Arm &a, bool verification) {
    contains(log, "packed launch bounds: " + std::to_string(a.blockThreads) +
                  " threads, " + std::to_string(a.minBlocks) + " min blocks\n");
    contains(log, "packed sigma fused: 1\n");
    contains(log, "packed sigma fused late y: 0\n");
    contains(log, "packed witness: 0\n");
    contains(log, "packed shared sigma: 1\n");
    contains(log, "packed table walk: 0 (0 branches, 0 shared bytes)\n");
    contains(log, "packed kernel: 126 registers/thread, 0 local bytes/thread, 1792 shared bytes/block");
    if (verification)
        contains(log, "device: NVIDIA RTX PRO 6000 Blackwell Server Edition, 188 SMs, " +
                      std::to_string(a.residentBlocks) + " block(s) of " +
                      std::to_string(a.blockThreads) + " packed threads resident per SM\n");
    require(log.find("MISMATCH") == std::string::npos && log.find("OVERFLOW") == std::string::npos &&
            log.find("collision found") == std::string::npos && log.find("solved") == std::string::npos,
            "forbidden runtime marker");
}

uint64_t finalWork(const std::string &log) {
    const std::regex progress("([0-9]+) iterations[ ]+([0-9]+) dp[ ]+([0-9]+) stored[ ]+([0-9]+) dropped");
    uint64_t iterations = 0;
    for (auto it = std::sregex_iterator(log.begin(), log.end(), progress);
         it != std::sregex_iterator(); ++it) {
        iterations = std::stoull((*it)[1]);
        require(std::stoull((*it)[4]) == 0, "nonzero dropped count");
    }
    return iterations;
}

uint64_t verificationRecords(const std::string &log) {
    const std::regex terminal("finished: [0-9.]+ M it/s, ([0-9]+) distinguished points \\(300 verified against the reference, 0 dropped\\)");
    std::smatch match;
    require(std::regex_search(log, match, terminal), "missing verification terminal");
    return std::stoull(match[1]);
}

}  // namespace

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "usage: audit RESULTS\n";
        return 2;
    }
    try {
        const std::string dir = argv[1];
        const std::string host = read(dir + "/host.txt");
        contains(host, "NVIDIA RTX PRO 6000 Blackwell Server Edition");
        contains(host, "V13.3.73");

        uint64_t records = 0;
        for (const Arm &a : kArms) {
            const std::string build = read(dir + "/build-" + a.name + ".log");
            contains(build, "-DECC_THREADS=" + std::to_string(a.blockThreads));
            contains(build, "-DECC_MINBLOCKS=" + std::to_string(a.minBlocks));
            contains(build, "-DECC_PACKED_INLINE_POLY=3");
            contains(build, "-DECC_PACKED_SHARED_SIGMA=1");
            contains(build, "-DECC_SIGMA_FUSED=1");
            contains(build, "-DECC_SIGMA_FUSED_LATE_Y=0");
            contains(build, "Function properties for _ZN12eccPacked1314walkE10WalkParamsIjEPj\n"
                            "    0 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads");
            const std::string verify = read(dir + "/verify-" + a.name + ".log");
            commonMarkers(verify, a, true);
            contains(verify, "backend cuda-packed131: 96256 threads x 16 slots x 1 lanes = "
                             "1540096 walks, dp weight 48, 95 steps per launch\n");
            require(finalWork(verify) == 1024163840ull, "incorrect verification work");
            const uint64_t thisRecords = verificationRecords(verify);
            if (!records) records = thisRecords;
            require(thisRecords == records, "verification corpus count differs");
        }

        const std::string whole = read(dir + "/t256-whole.ck");
        require(read(dir + "/t512-whole.ck") == whole, "whole checkpoint differs");
        require(read(dir + "/t256-to-t512.ck") == whole, "t256-to-t512 checkpoint differs");
        require(read(dir + "/t512-to-t256.ck") == whole, "t512-to-t256 checkpoint differs");

        std::ifstream corpus(dir + "/corpus-sorted.bin", std::ios::binary);
        require(bool(corpus), "missing sorted corpus");
        std::array<unsigned char, 32> record{}, previous{};
        uint64_t corpusRecords = 0;
        while (corpus.read(reinterpret_cast<char *>(record.data()), record.size())) {
            if (corpusRecords) require(previous <= record, "corpus is not sorted");
            require(record[6] == 33 && record[7] == 0, "foreign run id");
            require(record[24] <= 7, "noncanonical 131-bit key");
            for (int i = 25; i < 32; ++i) require(record[i] == 0, "nonzero key padding");
            previous = record;
            ++corpusRecords;
        }
        require(corpus.eof() && corpus.gcount() == 0, "partial corpus record");
        require(corpusRecords == records, "sorted corpus count differs");

        std::ifstream samples(dir + "/samples.tsv");
        require(bool(samples), "missing samples");
        std::string line;
        std::getline(samples, line);
        int warmups = 0, aa = 0, ab = 0;
        std::map<std::string, int> variants;
        while (std::getline(samples, line)) {
            std::istringstream fields(line);
            std::array<std::string, 9> column;
            for (int i = 0; i < 8; ++i)
                require(bool(std::getline(fields, column[i], '\t')), "sample framing");
            std::getline(fields, column[8]);
            const std::string &phase = column[0];
            const std::string &variant = column[3];
            const Arm &a = (variant == "candidate") ? kArms[1] : kArms[0];
            require(column[5] == "403726925824" && column[6] == "0", "unequal timing work");
            const std::string logPath = dir + "/" + phase + "-" + column[1] + "-" +
                                        column[2] + "-" + variant + ".log";
            const std::string log = read(logPath);
            commonMarkers(log, a, false);
            contains(log, "backend cuda-packed131: 385024 threads x 16 slots x 1 lanes = "
                          "6160384 walks, dp weight 0, 1024 steps per launch\n");
            require(finalWork(log) == 403726925824ull, "incorrect timing work");
            const std::regex terminal("finished: ([0-9.]+) M it/s, 0 distinguished points \\(0 verified against the reference, 0 dropped\\)");
            std::smatch match;
            require(std::regex_search(log, match, terminal), "missing timing terminal");
            require(std::stod(match[1]) == std::stod(column[4]), "TSV rate differs from log");
            if (phase == "warmup") ++warmups;
            else if (phase == "aa") ++aa;
            else if (phase == "ab") ++ab;
            else throw std::runtime_error("unknown phase");
            ++variants[variant];
        }
        require(warmups == 4 && aa == 10 && ab == 10, "incomplete timing panel");
        require(variants["control"] == 7 && variants["candidate"] == 7 &&
                variants["control_a"] == 5 && variants["control_b"] == 5,
                "unexpected variant counts");

        std::cout << "{\"valid\":true,\"verification_arms\":2,\"replayed_per_arm\":300,"
                     "\"warmup_rows\":4,\"aa_rows\":10,\"ab_rows\":10,"
                     "\"updates_per_timing_row\":403726925824,\"sorted_records\":"
                  << corpusRecords
                  << ",\"run_id\":33,\"registers_per_thread\":126,"
                     "\"local_bytes_per_thread\":0,\"shared_bytes_per_block\":1792}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
