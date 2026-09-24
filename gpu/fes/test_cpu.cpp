// gpu/fes/test_cpu.cpp -- verify the FES Gray-code search on the host, with no
// GPU.  It compiles the same fes.cuh the CUDA/Metal kernels use, runs the
// per-sub-cube search for every prefix (exactly what one GPU thread does), and
// checks the collected solutions against (a) a from-scratch brute force and
// (b) the independently brute-forced set the generator emitted.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "fes.cuh"
#include "system.h"

static std::vector<uint64_t> brute_force() {
  std::vector<uint64_t> sols;
  for (uint64_t x = 0; x < (1ULL << FES_N); ++x) {
    if (fes_eval(FES_CONST, FES_LIN, FES_QUAD_TRI, FES_N, x) == 0) {
      sols.push_back(x);
    }
  }
  return sols;
}

static std::vector<uint64_t> fes_all(int prefix_bits) {
  const int low = FES_N - prefix_bits;
  std::vector<uint64_t> out(1u << 16);
  std::vector<uint64_t> df(FES_N);
  std::vector<uint64_t> sols;
  for (uint64_t p = 0; p < (1ULL << prefix_bits); ++p) {
    int found = 0;
    fes_search_subcube(FES_CONST, FES_LIN, FES_QUAD_TRI, FES_N, low, p,
                       out.data(), &found, (int)out.size(), df.data());
    for (int i = 0; i < found && i < (int)out.size(); ++i) sols.push_back(out[i]);
  }
  std::sort(sols.begin(), sols.end());
  return sols;
}

int main() {
  int failures = 0;

  auto brute = brute_force();
  std::sort(brute.begin(), brute.end());

  // The generator's independent solution set must match our brute force.
  std::vector<uint64_t> gen(FES_SOLUTIONS, FES_SOLUTIONS + FES_NUM_SOLUTIONS);
  std::sort(gen.begin(), gen.end());
  if (gen != brute) {
    printf("FAIL: generator solution set disagrees with brute force (%zu vs %zu)\n",
           gen.size(), brute.size());
    ++failures;
  }

  // The planted solution must be a genuine zero.
  if (fes_eval(FES_CONST, FES_LIN, FES_QUAD_TRI, FES_N, FES_PLANTED) != 0) {
    printf("FAIL: planted solution is not a zero\n");
    ++failures;
  }

  // The sharded Gray-code search must find exactly the brute-force set, at
  // several shard granularities (different thread counts).
  for (int pb : {0, 2, 4, 6}) {
    if (pb >= FES_N) continue;
    auto got = fes_all(pb);
    if (got != brute) {
      printf("FAIL: FES search with %d prefix bits found %zu, expected %zu\n", pb,
             got.size(), brute.size());
      ++failures;
    } else {
      printf("ok: %d prefix bits (%d threads) -> %zu solutions\n", pb, 1 << pb,
             got.size());
    }
  }

  if (failures == 0) {
    printf("PASS: FES host verification (n=%d, m=%d, %zu solutions)\n", FES_N,
           FES_M, brute.size());
  }
  return failures == 0 ? 0 : 1;
}
