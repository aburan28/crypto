# Module-LWE CUDA reference

Computes `t = A*s + e` in `R_q^k`, where `R_q = Z_3329[x]/(x^256+1)` and
`k` is 2, 3, or 4. It accepts arrays of integer coefficients and supports
independent batches. Pass `A` as canonical residues in `[0, 3328]` and `s,e`
with absolute coefficient values at most 3328. All outputs are canonical
residues in `[0, 3328]`.

This is the mathematical key-generation equation used by Kyber/ML-KEM, not a
complete or production-ready ML-KEM implementation. It does not implement
seed expansion, centered-binomial sampling, NTT representation, packing,
encapsulation, decapsulation, or constant-time handling of secret material.
The direct convolution costs O(batch*k^2*256^2); use an ML-KEM NTT pipeline
before making throughput claims.

```
cmake -S . -B build
cmake --build build -j
./build/mlwe_cpu_check
./build/mlwe_cuda_check  # created if a CUDA compiler is available
```

If CMake is unavailable, run the CPU check with
`g++ -O2 -std=c++17 cpu_check.cpp -o cpu_check && ./cpu_check`.

Direct CUDA build: `nvcc -O2 -std=c++17 -arch=sm_80 mlwe_cuda.cu -o mlwe_cuda_check`.
Choose an architecture supported by the target GPU and CUDA toolkit.

The CUDA checker compares all output coefficients to the independent CPU
convolution for batch=8 and k=2,3,4. The CPU checker includes a negacyclic
wrap fixture, large-residue tests against a closed-form answer, and selected
coefficient checks on deterministic random inputs. The CUDA launcher supports
up to 65535 instances per batch.

The repository CI runs the CPU check and compiles the CUDA code. Compilation
does not verify GPU execution; run `mlwe_cuda_check` on a CUDA GPU before
accepting device correctness. No throughput or complete ML-KEM claims are made.
