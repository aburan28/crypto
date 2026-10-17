#!/bin/bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
cc -O3 -march=native -std=c11 -D_POSIX_C_SOURCE=200809L -fPIC -shared -I"$root/upstream/include" -I"$root/upstream/src" "$root/upstream"/src/{modarith,util,group_common,group_zp,group_ec,ec_tau,rho,bsgs}.c -pthread -lm -o "$root/native/libcryptanalysis.so"
