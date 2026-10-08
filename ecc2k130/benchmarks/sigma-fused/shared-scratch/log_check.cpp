// Fail-closed correctness-log checker for the frozen shared-scratch producer.
#include <charconv>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

void need(bool condition, const std::string &message) {
    if (!condition) fail(message);
}

template<class T>
T integer(const std::string &value, const std::string &label) {
    T result{};
    const char *first = value.data();
    const char *last = first + value.size();
    const auto parsed = std::from_chars(first, last, result);
    need(first != last && parsed.ec == std::errc{} && parsed.ptr == last,
         "invalid " + label);
    return result;
}

std::string read(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    need(bool(input), "cannot read " + path);
    return std::string(std::istreambuf_iterator<char>(input), {});
}

void validate(const std::string &text, std::uint64_t expectedIterations,
              unsigned expectedVerified, long long expectedResume) {
    static const std::regex finishedPattern(
        R"(^\s*finished: [0-9]+(?:\.[0-9]+)? M it/s.*$)");
    static const std::regex iterationPattern(R"(([0-9]+) iterations)");
    static const std::regex countsPattern(
        R"(\(([0-9]+) verified against the reference, ([0-9]+) dropped\))");
    static const std::regex resumePattern(R"( at iteration ([0-9]+)$)");
    std::istringstream input(text);
    std::string line;
    unsigned finishedCount = 0, countsCount = 0, resumeCount = 0;
    std::uint64_t finalIterations = 0;
    unsigned verified = ~0u, dropped = ~0u;
    long long resume = -1;
    bool sawIterations = false;
    while (std::getline(input, line)) {
        std::smatch match;
        if (std::regex_match(line, match, finishedPattern)) ++finishedCount;
        for (std::sregex_iterator iterator(line.begin(), line.end(), iterationPattern), end;
             iterator != end; ++iterator) {
            finalIterations = integer<std::uint64_t>(
                (*iterator)[1].str(), "iteration count");
            sawIterations = true;
        }
        for (std::sregex_iterator iterator(line.begin(), line.end(), countsPattern), end;
             iterator != end; ++iterator) {
            verified = integer<unsigned>((*iterator)[1].str(), "verified count");
            dropped = integer<unsigned>((*iterator)[2].str(), "dropped count");
            ++countsCount;
        }
        if (std::regex_search(line, match, resumePattern)) {
            resume = integer<long long>(match[1].str(), "resume iteration");
            ++resumeCount;
        }
    }
    need(finishedCount == 1, "expected exactly one finished row");
    need(sawIterations && finalIterations == expectedIterations,
         "final iteration count mismatch");
    need(countsCount == 1 && verified == expectedVerified && dropped == 0,
         "verification/drop count mismatch");
    if (expectedResume < 0)
        need(resumeCount == 0, "unexpected resume marker");
    else
        need(resumeCount == 1 && resume == expectedResume,
             "resume boundary mismatch");
    for (const char *bad : {"MISMATCH", "OVERFLOW", "collision found", "solved",
                            "stopping:", "unusable"})
        need(text.find(bad) == std::string::npos,
             "forbidden runtime marker " + std::string(bad));
}

std::string fixture(std::uint64_t iterations, unsigned verified, unsigned dropped,
                    long long resume = -1) {
    std::ostringstream output;
    if (resume >= 0) output << "resumed from checkpoint at iteration " << resume << '\n';
    output << "  1.0 s  1000.000 M it/s  " << iterations
           << " iterations  0 dp  0 stored  " << dropped << " dropped\n"
           << "  finished: 1000.000 M it/s, 0 distinguished points ("
           << verified << " verified against the reference, " << dropped
           << " dropped)\n";
    return output.str();
}

void rejected(const std::string &text, std::uint64_t iterations,
              unsigned verified, long long resume, const std::string &label) {
    bool failed = false;
    try {
        validate(text, iterations, verified, resume);
    } catch (...) {
        failed = true;
    }
    need(failed, label + " was accepted");
}

void selfTest() {
    constexpr std::uint64_t work = 1024163840ull;
    validate(fixture(work, 300, 0), work, 300, -1);
    validate(fixture(2339280, 0, 0, 380), 2339280, 0, 380);
    rejected(fixture(work - 1, 300, 0), work, 300, -1, "truncated work");
    rejected(fixture(work, 300, 0) + fixture(work, 300, 0),
             work, 300, -1, "duplicate finish");
    rejected(fixture(work, 300, 1), work, 300, -1, "drop");
    rejected(fixture(work, 299, 0), work, 300, -1, "verify drift");
    rejected(fixture(2339280, 0, 0, 379), 2339280, 0, 380, "resume drift");
    rejected(fixture(work, 300, 0, 380), work, 300, -1, "unexpected resume");
    std::cout << "PASS: exact work, one finish, verification/drop and resume-boundary negative controls\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            selfTest();
            return 0;
        }
        if (argc != 5) {
            std::cerr << "usage: log_check LOG EXPECTED_ITERATIONS EXPECTED_VERIFIED EXPECTED_RESUME_OR_MINUS1\n";
            return 2;
        }
        validate(read(argv[1]), integer<std::uint64_t>(argv[2], "expected iterations"),
                 integer<unsigned>(argv[3], "expected verified"),
                 integer<long long>(argv[4], "expected resume"));
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "LOG CHECK FAIL: " << error.what() << '\n';
        return 1;
    }
}
