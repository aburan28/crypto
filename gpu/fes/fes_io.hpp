// gpu/fes/fes_io.hpp -- the stable text I/O contract every FES worker speaks,
// so the Rust binary can drive a CUDA, Metal, or host-emulation worker through
// the exact same wire. Host-only (plain C++), included by fes_solve.cpp and by
// fes.cu's host side.
//
// Input file (whitespace-separated decimal integers):
//     n m
//     <const>
//     <lin_0> ... <lin_{n-1}>
//     <quad_tri_0> ... <quad_tri_{T-1}>        T = n*(n+1)/2
// where each coefficient is a u64 whose bit e is that term's coefficient in
// equation e (see fes.cuh / fref.py). Output on stdout:
//     SOLUTIONS <k>
//     <x_0>
//     ...
//     <x_{k-1}>
// each x a decimal u64 assignment that zeroes the system. A worker may cap how
// many it lists; the Rust side re-verifies every returned x, so a worker can
// only ever propose candidates.
#pragma once

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct FesSystem {
  int n = 0;
  int m = 0;
  uint64_t cst = 0;
  std::vector<uint64_t> lin;      // length n
  std::vector<uint64_t> quad_tri; // length n*(n+1)/2
};

// Parse a system file. Returns false on any malformed input.
inline bool fes_read_system(const std::string &path, FesSystem &sys) {
  std::ifstream f(path);
  if (!f) return false;
  if (!(f >> sys.n >> sys.m)) return false;
  if (sys.n < 0 || sys.n > 62 || sys.m < 0 || sys.m > 64) return false;
  if (!(f >> sys.cst)) return false;
  sys.lin.assign(sys.n, 0);
  for (int i = 0; i < sys.n; ++i)
    if (!(f >> sys.lin[i])) return false;
  const int t = sys.n * (sys.n + 1) / 2;
  sys.quad_tri.assign(t, 0);
  for (int i = 0; i < t; ++i)
    if (!(f >> sys.quad_tri[i])) return false;
  return true;
}

// Print solutions in the contract's format.
inline void fes_write_solutions(const std::vector<uint64_t> &sols) {
  printf("SOLUTIONS %zu\n", sols.size());
  for (uint64_t x : sols) printf("%llu\n", (unsigned long long)x);
}
