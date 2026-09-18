# 29 B/s on this RTX PRO 6000

Walk-rate engineering, not an ECDLP exponent claim. Boundaries and the
falsification target: [THROUGHPUT-29B.md](../../THROUGHPUT-29B.md). Frozen
numbers: [summary.json](summary.json).

The 20 B/s attempt on this SKU stopped at 17.298 B/s. This round asked
whether leftover product ALU, the ONB inverse chain, a slot pipeline, or
(after Nsight) moving the polynomial square off `clmad` can clear **29 B
complete scalar updates/s**.

## Result

**No.** Best verified median **17.414 B/s** (`ALU_SQUARE=1` on the 17.298
control). 0.600 of 29, 0.792 of the 22 B one-add floor. 300/300 reports, 0
dropped, 5.3125 products/update. The first-round hoist / ONB-inv / pipeline
arms were slower than the control.

Nsight Compute on the 17.298 kernel: DRAM 11%, ALU pipe 36%, FP64/`clmad`
71%. That is why `ALU_SQUARE` paid (+0.67%) and the ALU-only knobs did not.

## Acceptance

Success is a verified median **> 29.0 B/s**. Inadmissible: two GPUs as one,
a different SKU, dropping identity, changing the product count.

```sh
bash benchmarks/throughput-29b-gpu/run.sh
```

Best row, same recipe knobs plus the square:

```sh
make -B ecc2k130 ARCH='-gencode arch=compute_120,code=sm_120' \
  BATCH=16 THREADS=256 MINBLOCKS=2 \
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 \
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 \
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 \
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 \
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 \
  PACKED_PAIR_ILP=1 PACKED_L2_PERSIST=1 UNROLL_SLOTS=2 \
  WALK_TABLE=1 TABLE_PIVOT_BYTES=1 PACKED_ALU_SQUARE=1
./ecc2k130 --curve 131 --packed --bench --steps 1024 --launches 32 --verify 0
```
