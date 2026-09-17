# Round 3 (ITERATION-FUNCTION.md section 7): build, verify and bench the
# inversion / denominator / batch variants against the section 6 kernels.
# Run inside nvidia/cuda:13.3.1-devel-ubuntu24.04 with the ecc2k130 tree
# mounted at /src:
#   docker run --rm --gpus all -v $PWD:/src -w /src nvidia/cuda:13.3.1-devel-ubuntu24.04 bash /src/gpujob-round3.sh
cd /src
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq make g++ python3 >/dev/null 2>&1
nvidia-smi --query-gpu=name,clocks.max.sm,power.limit,memory.total --format=csv,noheader
F='ARCH=-gencode arch=compute_120,code=sm_120'
K="THREADS=256 MINBLOCKS=2 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1"
echo "=== unit test: device arithmetic incl. the polynomial-basis inversion (test-packed-cuda, POLY_INV=1)"
make test-packed-cuda "$F" $K BATCH=16 PACKED_POLY_INV=1 2>&1 | grep -E "PASS|FAIL|mismatch|error" | head -8
echo "=== unit test: device arithmetic, normal-basis inversion (control)"
make test-packed-cuda "$F" $K BATCH=16 PACKED_POLY_INV=0 2>&1 | grep -E "inversion|FAIL|mismatch|error" | head -4
echo "=== unit test: device table-walk primitives"
make test-table-walk-cuda "$F" $K BATCH=16 2>&1 | grep -E "PASS|FAIL|mismatch|error" | head -4
build() { make -B ecc2k130 "$F" $K $2 2>&1 | grep -E "error|Used [0-9]+ registers|spill" | grep -A1 -B1 "eccPacked131.*walk" | grep -E "registers|spill" | head -2; mv ecc2k130 ecc2k130-$1 || exit 1; }
echo "=== build legacy16      (shipping walk, control)";                 build legacy16    "BATCH=16 WALK_TABLE=0"
echo "=== build legacy16pi    (H: shipping walk + POLY_INV=1)";          build legacy16pi  "BATCH=16 WALK_TABLE=0 PACKED_POLY_INV=1"
echo "=== build table16       (A: section 6 table kernel, control)";     build table16     "BATCH=16 WALK_TABLE=1"
echo "=== build table16pi     (B: + POLY_INV=1)";                        build table16pi   "BATCH=16 WALK_TABLE=1 PACKED_POLY_INV=1"
echo "=== build table16pi2    (C: + POLY_INV=2)";                        build table16pi2  "BATCH=16 WALK_TABLE=1 PACKED_POLY_INV=2"
echo "=== build table16pind   (D: B + DENOM_STORE=0)";                   build table16pind "BATCH=16 WALK_TABLE=1 PACKED_POLY_INV=1 TABLE_DENOM_STORE=0"
echo "=== build table24pind   (E: D at batch 24)";                       build table24pind "BATCH=24 WALK_TABLE=1 PACKED_POLY_INV=1 TABLE_DENOM_STORE=0"
echo "=== build table24pi2nd  (E2: E with POLY_INV=2)";                  build table24pi2nd "BATCH=24 WALK_TABLE=1 PACKED_POLY_INV=2 TABLE_DENOM_STORE=0"
echo "=== build table32pind   (F: D at batch 32)";                       build table32pind "BATCH=32 WALK_TABLE=1 PACKED_POLY_INV=1 TABLE_DENOM_STORE=0"
echo "=== build table24pi     (G: B at batch 24, denominators stored)";  build table24pi   "BATCH=24 WALK_TABLE=1 PACKED_POLY_INV=1"
ALL="legacy16 legacy16pi table16 table16pi table16pi2 table16pind table24pind table24pi2nd table32pind table24pi"
for b in $ALL; do
  echo "=== verify $b: 300 device reports re-walked by the host reference"
  timeout 1500 ./ecc2k130-$b --curve 131 --packed --dp-weight 48 --dp-cap 1048576 --steps 16 --launches 6 --verify 300 2>&1 | grep -E "MISMATCH|OVERFLOW|finished|table walk|inversion|backend"
  echo "exit: ${PIPESTATUS[0]}"
done
for rep in 1 2 3; do for b in $ALL; do
  echo "=== bench $b rep $rep"
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 32 --verify 0 2>&1 | grep -E "finished"
  nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader
done; done
echo "=== audited geometry: --threads 385024, 1024 steps x 32 launches"
for rep in 1 2 3; do for b in legacy16 table16 table16pi table16pind table24pind table24pi2nd table32pind table24pi; do
  echo "=== bench385k $b rep $rep"
  ./ecc2k130-$b --curve 131 --packed --bench --threads 385024 --steps 1024 --launches 32 --verify 0 2>&1 | grep -E "finished|resident"
  nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader
done; done
