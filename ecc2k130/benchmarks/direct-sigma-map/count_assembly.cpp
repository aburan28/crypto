#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

namespace {

bool startsWith(const std::string &text, const std::string &prefix) {
    return text.size() >= prefix.size() && text.compare(0, prefix.size(), prefix) == 0;
}

std::string trim(std::string text) {
    while (!text.empty() && std::isspace(static_cast<unsigned char>(text.front()))) text.erase(text.begin());
    while (!text.empty() && std::isspace(static_cast<unsigned char>(text.back()))) text.pop_back();
    return text;
}

}  // namespace

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "usage: count_assembly assembly.s\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "cannot open " << argv[1] << '\n';
        return 2;
    }
    std::map<std::string, long> counts;
    std::string current;
    std::string line;
    while (std::getline(in, line)) {
        const std::string stripped = trim(line);
        const size_t colon = stripped.find(':');
        if (colon != std::string::npos && stripped.substr(0, colon).find_first_of(" \t") == std::string::npos) {
            std::string label = stripped.substr(0, colon);
            if (!label.empty() && label.front() == '_') label.erase(label.begin());
            if (startsWith(label, "composed_j") || startsWith(label, "table3_j") ||
                startsWith(label, "half5_j")) current = label;
            else if (label.find("apply3_j") != std::string::npos) current = "generated_" + label;
            else if (label.find("toPolynomial131") != std::string::npos) current = "helper_toPolynomial131";
            else if (label.find("sigma131") != std::string::npos) current = "helper_sigma131";
            else if (label.find("fromPolynomialProduct131") != std::string::npos) current = "helper_fromPolynomialProduct131";
            else if (label.find("fromPolynomialReduced131") != std::string::npos) current = "helper_fromPolynomialReduced131";
            else if (!label.empty() && label.front() != 'L') current.clear();
            continue;
        }
        if (current.empty() || stripped.empty() || stripped.front() == '.' || stripped.front() == '#') continue;
        if (stripped.find(':') != std::string::npos && stripped.front() == 'L') continue;
        ++counts[current];
    }
    for (const auto &[name, count] : counts) std::cout << name << ' ' << count << '\n';
    return counts.empty() ? 1 : 0;
}
