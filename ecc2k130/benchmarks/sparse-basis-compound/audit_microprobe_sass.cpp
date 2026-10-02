#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

namespace {

struct Counts {
    int lds128 = 0;
    int ldsScalar = 0;
    int local = 0;
    int calls = 0;
};

[[noreturn]] void fail(const std::string &message) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
}

void finish(const std::string &name, const Counts &counts, int *kernels,
            std::ostream &out) {
    if (name.find("benchmarkKernel") == std::string::npos) return;
    ++*kernels;
    out << name << '\t' << counts.lds128 << '\t' << counts.ldsScalar << '\t'
        << counts.local << '\t' << counts.calls << '\n';
    if (counts.lds128 < 44 || counts.ldsScalar < 44 || counts.local != 0 ||
        counts.calls != 0)
        fail("instruction gate for " + name);
}

}  // namespace

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "usage: audit_microprobe_sass PROBE_SASS REPORT_TSV\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    if (!input || !output) fail("open input/output");
    output << "function\tlds128\tldsScalar\tlocalInstructions\tcalls\n";
    std::string current, line;
    Counts counts;
    int kernels = 0;
    while (std::getline(input, line)) {
        const std::string marker = "Function : ";
        const size_t at = line.find(marker);
        if (at != std::string::npos) {
            finish(current, counts, &kernels, output);
            current = line.substr(at + marker.size());
            counts = Counts{};
            continue;
        }
        if (line.find("benchmarkKernel") != std::string::npos &&
            line.find("/*") == std::string::npos) {
            finish(current, counts, &kernels, output);
            current = line;
            counts = Counts{};
            continue;
        }
        if (current.find("benchmarkKernel") == std::string::npos) continue;
        if (line.find("LDS.128") != std::string::npos) ++counts.lds128;
        else if (line.find("LDS") != std::string::npos &&
                 line.find("ULDS") == std::string::npos)
            ++counts.ldsScalar;
        if (line.find("LDL") != std::string::npos || line.find("STL") != std::string::npos)
            ++counts.local;
        if (line.find("CALL") != std::string::npos) ++counts.calls;
    }
    finish(current, counts, &kernels, output);
    if (kernels != 7) fail("expected seven benchmark kernels, found " + std::to_string(kernels));
    output << "PASS\tseven kernels\t44+ vector/scalar loads each\tzero local\tzero calls\n";
    std::cout << "PASS: seven table kernels contain LDS.128 and scalar LDS, with no local instructions or calls\n";
    return 0;
}
