# 50 B/s halving boundary

The selected 1024-thread lambda-halving primitive measures **29.053 B/s** and
passes all 512 subgroup differential cases. It does not meet 50 B/s and is not
a rho map by itself.

`result.log` contains the build, correctness check and three samples.
`result.json` records the CLMAD ceiling, ideal mixed-walk optimization and
fruitless-cycle status. The canonical report is
[`../../THROUGHPUT-50B.md`](../../THROUGHPUT-50B.md).
