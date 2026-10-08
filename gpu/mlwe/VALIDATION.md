# Validation status

This addition implements the forward Module-LWE equation. It makes no key
recovery, cryptanalytic advantage, or performance claim.

Local command:

```
g++ -O2 -std=c++17 -Wall -Wextra -pedantic cpu_check.cpp -o cpu_check
./cpu_check
```

Observed output, 2026-10-01 UTC:

```
PASS k=2 batch=3
PASS k=3 batch=3
PASS k=4 batch=3
```

These checks cover a known negacyclic wrap, large residues with a closed-form
answer, and selected coefficients from deterministic seeded fixtures. No
wall-clock measurement was taken.

CUDA compilation and GPU execution were unavailable locally: neither nvcc
nor an NVIDIA GPU was present. The PR adds a CUDA compile CI job; its result
must be evaluated separately from GPU execution. All-coefficient CUDA/CPU
comparison on a GPU remains pending and is the acceptance gate for device
correctness. The PR should remain a draft until that gate is satisfied.
