#!/usr/bin/env bash
set -euo pipefail
cd /workspace/n83_optimization/checks/montgomery_native_count
g++ -std=c++17 -O2 -Wall -Wextra -Werror -pedantic count.cpp -o count 2> compile.log
./count > run.log
sha256sum count.cpp count > SHA256SUMS
