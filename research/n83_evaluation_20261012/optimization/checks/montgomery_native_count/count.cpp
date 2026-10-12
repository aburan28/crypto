#include <cstdint>
#include <iostream>

// Standalone enumeration of the declared Montgomery model.  This program
// imports no repository code and uses no elliptic-curve library or conversion.
int main() {
    constexpr std::int64_t p = 101;
    constexpr std::int64_t A = 3;
    constexpr std::int64_t B = 1;
    std::uint64_t tested = 0;
    std::uint64_t affine = 0;
    std::uint64_t zero_y_points = 0;

    for (std::int64_t x = 0; x < p; ++x) {
        for (std::int64_t y = 0; y < p; ++y) {
            ++tested;
            const std::int64_t lhs = (B * y * y) % p;
            const std::int64_t rhs = (x * x * x + A * x * x + x) % p;
            if (lhs == rhs) {
                ++affine;
                if (y == 0) {
                    ++zero_y_points;
                }
            }
        }
    }

    // The nonsingular projective Montgomery cubic has one point at infinity.
    constexpr std::uint64_t infinity = 1;
    const std::uint64_t total = affine + infinity;
    std::cout << "equation=B*y^2=x^3+A*x^2+x (mod p)\n"
              << "p=" << p << "\nA=" << A << "\nB=" << B << '\n'
              << "affine_pairs_tested=" << tested << '\n'
              << "affine_points=" << affine << '\n'
              << "zero_y_points=" << zero_y_points << '\n'
              << "points_at_infinity=" << infinity << '\n'
              << "total_points=" << total << '\n';
    return tested == 10201 && total == 96 ? 0 : 1;
}
