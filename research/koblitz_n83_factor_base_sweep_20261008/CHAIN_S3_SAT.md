# Higher-arity chained-S3 SAT source gate

`binary_semaev_chain_sat.rs` adds a full-width native-XOR circuit for the
binary-curve relation `S3(x,y,z)=xyz+(xy+xz+yz)^2+b=0`. Four such nodes encode
five summands and five nodes encode six summands. Every field multiplication
is a set of two-input AND gates with exact polynomial-basis reduction; the
linear output and S3 equations are native XOR rows. Rewriting
`xy+xz+yz = xy+(x+y)z` needs only three field products per affine S3 node;
the square is a linear XOR transform. The module accepts field
degrees through 127, including the study's degree-83 polynomial, but this
change constructs **no retained N83 instance**.

For `m` summands, there are `m-2` internal point sums. The identity has no
x-coordinate, so the producer searches all `2^(m-2)` intermediate-identity
patterns: eight for m=5 and sixteen for m=6. An affine node gets an S3
equation. An identity output requires equal input x-coordinates; adding an
affine summand to an identity input copies its x-coordinate to the output.
The impossible identity-plus-affine-to-identity case is refuted immediately.
Every SAT model is re-evaluated over the original binary field, checked
against the exact finite factor-base coordinate set, and lifted over all
available point signs in the original curve group. Internal S3 roots or
trace-inadmissible algebraic solutions cannot become accepted relations.
If a model fails group lifting, all of its intermediate x assignments share
the same summand tuple, so one blocking clause excludes that tuple. A
refutation requires every pattern to finish; model, conflict, variable or
domain-clause caps return an inconclusive result.

The pre-domain all-affine source count at n=83 is 84,743 SAT variables for
m=5 and 105,908 for m=6. The respective circuits have 82,668 and 103,335
two-input AND gates, each represented by three CNF clauses, before exact
finite-domain constraints. These are construction counts, not memory or
solving-time observations. The finite-coordinate trie can dominate memory:
the producer computes its exact missing-prefix clause count on a given base
before installing any domain clauses and applies a caller-supplied cap. An
external hard memory and wall guard remains necessary for any retained-base
construction or search.

`sat_decompose_chained_s3` is a separate experimental library entry point
for ambient-basis explicit-orbit bases and m=5 or m=6. It is not silently
substituted for the frozen sweep's SAT disposition, and no retained-base
relation, natural rank, full cold time or winner follows from its small-curve
checks. The focused correctness suite also fixes both operands above bit 64
in the pinned degree-83 polynomial and compares every product-output bit
against independent field multiplication; it does not build a retained-base
model or time a relation search. Its next empirical gate is a source-pinned,
cgroup-guarded retained
model-construction receipt with the exact stored base domain. A subsequent
search must record unsuccessful patterns, solver work, original field/group
replay and rank novelty before the backend can enter a versioned sweep.

## Retained-base construction gate

`chain_capacity_supervisor.py` now accepts a replayed primary panel object,
K=64/256/600, a frozen policy and seed, arity five or six, and one
intermediate-identity mask. It freezes a clean source commit and the five
source files used by the worker, hashes a Linux ELF binary, and launches the
worker without network under a hard Docker memory ceiling, zero swap and an
external process-wall deadline. The worker independently checks the cgroup
and compiled source snapshot, replays the selected base and public target,
preflights the exact finite-coordinate trie, and either records a variable or
domain-clause cap or constructs that one SAT pattern. It never calls the SAT
solver. The outer receipt retains timeouts, memory exits, changed frozen
inputs and malformed worker receipts as separate outcomes. Synthetic guard
tests exercise those paths, but no retained N83 model has been launched under
this gate. A source estimate or capped outcome cannot rank a factor base.

After a separately authorized resource tranche, use an isolated clean
checkout and a source-matched static Linux worker. The invocation is:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/chain_capacity_supervisor.py \
  --panel /absolute/retained/pilot-01 --columns 64 \
  --policy public_x_hash --seed 2026100801 --summands 5 --identity-mask 0 \
  --max-variables 100000 --max-domain-clauses 1000000 \
  --wall-seconds 60 --memory-mib 4096 --binary /absolute/linux/worker \
  --output-dir /absolute/new/receipt-directory
```

The initial gate should use one pattern and one base. Expansion to every
identity pattern, solving, rank, and a versioned cross-base sweep depends on
its measured memory and wall receipt. The capacity wall time under Docker
emulation is not a matched cold-runtime measurement.

## Generic index-calculus driver path

`DecompositionStrategy::ChainedS3` now passes a checked five- or
six-summand witness to the existing relation-row rewrite, full-width
incremental rank tracker and target verification. Its four hard limits are
explicit in `KoblitzIcOptions::chain_s3_limits`; the default limits are zero,
so selecting this strategy without a resource configuration returns an
inconclusive disposition before sampling trials or launching SAT search. The
driver records solver calls, models, conflicts, inconclusive
targets, original-field/domain model failures and valid x-models rejected by
original-curve point-sign lifting separately. A model replay failure stops
the run before a rank or target claim. The factor-base-log and individual-log
entry points use the same strategy dispatch. A public degree-7
explicit-orbit fixture checks a complete relation-to-rank-to-verified-target
path; this validates plumbing and does not measure the retained degree-83
capacity or cold runtime.

## Guarded retained-base search path

The separate `primary-chain-cold` worker command now selects exactly one
replayed primary public-x object by policy and seed. It accepts the retained
v1 K=64/256/600 panel or a separately replayed v2 one-object panel at
K=1,182/2,048/4,096/8,192/16,627. The v2 manifest must bind the frozen size
design and its matching replay receipt. The public target corpus and
known-answer sidecar are supplied separately from the original public panel,
so a v2 object does not have to copy or alter those fixtures. It rejects
unsupported arity and zero or excessive SAT caps before making a run directory,
checks the declared Linux cgroup memory limit and zero swap, and checks a
clean source snapshot against the compiled exporter, adapter, S3 circuit,
generic index-calculus driver and v2 size design. It passes explicit variable, finite-domain
clause, model and conflict limits to `DecompositionStrategy::ChainedS3`.
The public fixture is the only target solver input; the known-answer sidecar
is opened only after a candidate scalar is produced. Progress events and
separate rank, solver-model, group-rejected-model and failure counts are
retained. A solver cap is `UNKNOWN_solver_cap`, not a refutation.

`chain_search_supervisor.py` is the outer entry point. It freezes the source
files, panel and public-fixture hashes, and binary hash, pins a locally present Docker image ID, disables
network, limits CPU, memory, swap and wall time, and retains an outer receipt
for success, timeout, resource exit, changed inputs or malformed worker
output. Synthetic tests cover the guard and reject a false total-runtime or
rank claim. `PASS_verified_target_only` is a checked target result within this
stage; column-log replay, fully charged factor-base precompute, artifact I/O,
matched rho and independent-host replay remain outside it. The local pilot
budget is exhausted, so no retained N83 search has been launched through this
new path.

After a separate compute grant, a first capacity receipt and a source-matched
Linux worker, the guarded preflight command is:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/chain_search_supervisor.py \
  --panel /absolute/retained/pilot-01 \
  --fixtures /absolute/retained/pilot-01 --columns 64 \
  --policy public_x_hash --seed 2026100801 --summands 5 \
  --max-trials 0 --max-variables 100000 --max-domain-clauses 1000000 \
  --max-models 1 --conflict-budget 10000 \
  --wall-seconds 60 --memory-mib 4096 --binary /absolute/linux/worker \
  --output-dir /absolute/new/search-receipt-directory
```

Only after that preflight and the independent capacity gate pass should a
new, separately bounded run raise `--max-trials` above zero. Its negative
and capped outcomes remain evidence about this source and base only; they do
not identify the minimum-runtime factor base.

For a larger v2 base, point `--panel` at its own freshly replayed one-object
directory and leave `--fixtures` at the frozen public corpus directory. The
guard accepts no run until both panel and fixture receipts are present. The
larger base must first pass construction, generic replay and upload/download
verification under a new allocation; none has been built yet.
