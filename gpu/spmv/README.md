# ECDLP sparse matrix-times-block workers

These workers implement the `SPMV1` contract used by
`cryptanalysis::koblitz_sparse_la::SpmvBackend::Worker`.

```sh
make -C gpu/spmv host-worker
IC_SPMV_WORKER=$PWD/gpu/spmv/spmv_solve cargo test --lib koblitz_sparse_la

# With CUDA installed:
make -C gpu/spmv cuda-worker
IC_SPMV_WORKER=$PWD/gpu/spmv/spmv_cuda cargo test --lib koblitz_sparse_la
```

The Rust caller re-computes and compares every worker product with its portable
CPU reference. A worker failure or mismatch therefore cannot create a false
relation-matrix solution; it only causes a CPU fallback.

`SPMV1` is whitespace-delimited:

1. `SPMV1 rows columns lanes modulus nonzeros`
2. `row_ptr[rows + 1]`
3. `column_index[nonzeros]`
4. `coefficient[nonzeros]`
5. `x[columns * lanes]`

The worker prints `SPMV1 rows lanes`, followed by `rows * lanes` residues.
Rows are independent, so a cluster launcher can split the CSR row range across
workers and concatenate results in row order without any modular reduction
between shards. The in-process `sharded` backend uses the same partitioning.

The CUDA kernel uses exact add-and-double multiplication for moduli below
`2^63`. It is intentionally a correctness-first backend; hardware-specific
Barrett/Montgomery tuning needs a separately preregistered performance round.
