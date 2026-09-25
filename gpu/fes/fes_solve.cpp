// gpu/fes/fes_solve.cpp -- host-emulation FES worker. Reads a system in the
// shared contract (fes_io.hpp), runs the exact fes.cuh Gray-code search on the
// CPU, and prints the solutions. It needs no GPU, so it is the worker the CI
// exercises the Rust binary's subprocess path against, and the CPU fallback a
// GPU-less machine uses. The CUDA (fes.cu) and Metal (fes_metal.mm) workers
// speak the same contract, so the Rust side drives all three identically.
//
//   fes_solve --in <system-file>      (also accepts a bare positional path)
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "fes.cuh"
#include "fes_io.hpp"

int main(int argc, char **argv) {
  // Accept `--in <file>` (the shared worker convention) or a bare positional
  // path, so every worker — fes_solve, fes_cuda, fes_metal — takes the same
  // arguments from the Rust driver.
  const char *in_path = nullptr;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--in") && i + 1 < argc) {
      in_path = argv[++i];
    } else if (argv[i][0] != '-') {
      in_path = argv[i];
    }
  }
  if (!in_path) {
    fprintf(stderr, "usage: %s --in <system-file>\n", argv[0]);
    return 2;
  }
  FesSystem sys;
  if (!fes_read_system(in_path, sys)) {
    fprintf(stderr, "fes_solve: could not parse %s\n", in_path);
    return 2;
  }

  // One sub-cube covering the whole space (low_bits = n), which is what a
  // single GPU thread with no prefix would run; the algorithm is identical.
  const int max_out = 1 << 20;
  std::vector<uint64_t> out(max_out);
  std::vector<uint64_t> df(sys.n > 0 ? sys.n : 1);
  int found = 0;
  fes_search_subcube(sys.cst, sys.lin.data(), sys.quad_tri.data(), sys.n, sys.n,
                     0, out.data(), &found, max_out, df.data());
  int listed = found < max_out ? found : max_out;
  std::vector<uint64_t> sols(out.begin(), out.begin() + listed);
  fes_write_solutions(sols);
  return 0;
}
