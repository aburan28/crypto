#include <cstdint>
#include <iostream>

// Independent direct root enumeration for the separately corrected field fixture.
int main() {
    constexpr std::int64_t p = 37;
    std::uint64_t tested = 0;
    std::uint64_t roots = 0;
    for (std::int64_t x = 0; x < p; ++x) {
        ++tested;
        if ((x * x + 2) % p == 0) {
            ++roots;
            std::cout << "root=" << x << '\n';
        }
    }
    std::cout << "polynomial=x^2+2\np=" << p << '\n'
              << "field_elements_tested=" << tested << '\n'
              << "roots=" << roots << '\n';
    // A quadratic over a field is irreducible precisely when it has no root.
    return tested == 37 && roots == 0 ? 0 : 1;
}
