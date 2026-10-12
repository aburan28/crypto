#!/usr/bin/env bash
set -euo pipefail
cd /workspace/n83_optimization/checks/montgomery_native_count
g++ -std=c++17 -O2 -Wall -Wextra -Werror -pedantic modulus37.cpp -o modulus37 2> compile_modulus37.log
./modulus37 > run_modulus37.log
sha256sum modulus37.cpp modulus37 > modulus37.SHA256SUMS
