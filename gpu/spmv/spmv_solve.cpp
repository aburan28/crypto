#include "spmv_io.hpp"

#include <iostream>

static std::uint64_t add_mod(
    std::uint64_t a, std::uint64_t b, std::uint64_t modulus) {
    const auto sum = a + b; // both inputs are below 2^63
    return sum >= modulus ? sum - modulus : sum;
}

static std::uint64_t multiply_mod(
    std::uint64_t a, std::uint64_t b, std::uint64_t modulus) {
    std::uint64_t result = 0;
    while (b != 0) {
        if (b & 1) result = add_mod(result, a, modulus);
        b >>= 1;
        if (b != 0) a = add_mod(a, a, modulus);
    }
    return result;
}

int main(int argc, char **argv) {
    try {
        const auto p = read_spmv1(input_path(argc, argv));
        std::vector<std::uint64_t> y(p.rows * p.lanes, 0);
        for (std::uint64_t row = 0; row < p.rows; ++row) {
            for (std::uint64_t at = p.row_ptr[row]; at < p.row_ptr[row + 1]; ++at) {
                const auto column = p.column_index[at];
                const auto coefficient = p.coefficient[at] % p.modulus;
                for (std::uint64_t lane = 0; lane < p.lanes; ++lane) {
                    auto &out = y[row * p.lanes + lane];
                    const auto product = multiply_mod(
                        coefficient, p.x[column * p.lanes + lane], p.modulus);
                    out = add_mod(out, product, p.modulus);
                }
            }
        }
        print_spmv1(p, y);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
