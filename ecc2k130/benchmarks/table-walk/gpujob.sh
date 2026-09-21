# Three-way build, verify and bench job behind benchmarks/table-walk/comparison.json.
cd /src
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq make g++ python3 >/dev/null 2>&1
nvidia-smi --query-gpu=name,clocks.max.sm,power.limit --format=csv,noheader
F='ARCH=-gencode arch=compute_120,code=sm_120'
K="BATCH=16 THREADS=256 MINBLOCKS=2 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1"
echo "=== unit test: device table-walk primitives (LUT selection)"
make test-table-walk-cuda "$F" $K 2>&1 | grep -v "^ptxas\|^$\|^nvcc\|^mkdir\|^    -X\|^\./build" | tail -6
build() { make -B ecc2k130 "$F" $K $2 2>&1 | grep -E "error|Used [0-9]+ registers|spill" | grep -A1 -B1 "eccPacked131.*walk" | grep -E "registers|spill" | head -2; mv ecc2k130 ecc2k130-$1 || exit 1; }
echo "=== build legacy";   build legacy "WALK_TABLE=0"
echo "=== build table";    build table  "WALK_TABLE=1"
echo "=== build table-alusq"; build alusq "WALK_TABLE=1 PACKED_ALU_SQUARE=1"
for b in table alusq legacy; do
  echo "=== verify $b: 300 device reports re-walked by the host reference"
  ./ecc2k130-$b --curve 131 --packed --dp-weight 50 --steps 16 --launches 6 --verify 300 2>&1 | grep -E "MISMATCH|finished|table walk|resident"
done
for rep in 1 2 3; do for b in legacy table alusq; do
  echo "=== bench $b rep $rep"
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 32 --verify 0 2>&1 | grep -E "finished"
  nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader
done; done
# Run inside nvidia/cuda:13.3.1-devel-ubuntu24.04 with the ecc2k130 tree mounted at /src:
#   docker run --rm --gpus all -v $PWD:/src -w /src nvidia/cuda:13.3.1-devel-ubuntu24.04 bash /src/gpujob.sh
