#pragma once

#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct SpmvProblem {
    std::uint64_t rows{};
    std::uint64_t columns{};
    std::uint64_t lanes{};
    std::uint64_t modulus{};
    std::uint64_t nonzeros{};
    std::vector<std::uint64_t> row_ptr;
    std::vector<std::uint64_t> column_index;
    std::vector<std::uint64_t> coefficient;
    std::vector<std::uint64_t> x;
};

inline SpmvProblem read_spmv1(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open SPMV1 input");
    std::string magic;
    SpmvProblem p;
    in >> magic >> p.rows >> p.columns >> p.lanes >> p.modulus >> p.nonzeros;
    if (magic != "SPMV1") throw std::runtime_error("bad SPMV1 header");
    if (p.modulus < 2 || p.modulus >= (std::uint64_t{1} << 63))
        throw std::runtime_error("SPMV1 modulus is outside [2, 2^63)");
    p.row_ptr.resize(p.rows + 1);
    p.column_index.resize(p.nonzeros);
    p.coefficient.resize(p.nonzeros);
    p.x.resize(p.columns * p.lanes);
    for (auto &v : p.row_ptr) in >> v;
    for (auto &v : p.column_index) in >> v;
    for (auto &v : p.coefficient) in >> v;
    for (auto &v : p.x) in >> v;
    if (!in || p.row_ptr.front() != 0 || p.row_ptr.back() != p.nonzeros)
        throw std::runtime_error("truncated or inconsistent SPMV1 input");
    for (std::uint64_t row = 0; row < p.rows; ++row)
        if (p.row_ptr[row] > p.row_ptr[row + 1])
            throw std::runtime_error("SPMV1 row pointers are not monotone");
    for (auto column : p.column_index)
        if (column >= p.columns) throw std::runtime_error("SPMV1 column out of range");
    return p;
}

inline std::string input_path(int argc, char **argv) {
    if (argc != 3 || std::string(argv[1]) != "--in")
        throw std::runtime_error("usage: worker --in FILE");
    return argv[2];
}

inline void print_spmv1(const SpmvProblem &p, const std::vector<std::uint64_t> &y) {
    if (y.size() != p.rows * p.lanes) throw std::runtime_error("wrong output size");
    std::cout << "SPMV1 " << p.rows << ' ' << p.lanes << '\n';
    for (std::size_t i = 0; i < y.size(); ++i) {
        if (i) std::cout << ' ';
        std::cout << y[i];
    }
    std::cout << '\n';
}
