# IC curve records and ICV1 crosswalk

[`curves.yaml`](curves.yaml) and its [schema](curves.schema.json) mirror the exact small curve registry in
`aburan28/cryptanalysis` at
`experiments/ic-candidate-catalog/curves.yaml`. The mirrors have the same
bytes and retain full `EC1` curve UIDs, exact field/curve records,
known/unknown trait statuses and isogeny route references. Evidence paths
inside it are relative to **cryptanalysis**. Change both copies in one
paired update; compare their SHA-256 before merging either update.

Run `python3 docs/curves/ic/validate_semantics.py` to check the mirrored YAML,
JSON Schemas, exact EC1/UID hashes, ICV1 crosswalk, trait statuses and typed
link references. `mirror-lock.json` pins the three mirrored source files.
The `ic-semantic-metadata` CI job runs these checks and compares their bytes
and the validator implementation against cryptanalysis `main`; its daily run
also detects later mirror drift. Keep the negative tests mirrored as well.
Merge a cryptanalysis source change before updating this mirror. The check
validates metadata and evidence references, not mathematical certificates or
measured performance.

`ICV1` is crypto's model/display name. `EC1` plus the full
`urn:ec-record:1:sha256:...` UID identifies the field representation,
model, subgroup and generator used by a DLP experiment. The toy N13 and
ECC2K-130 source entries match existing ICV1 **models**, but their EC1
representations are not currently entries in
[`registry.json`](../registry.json). The degree-263 codomain has an exact
EC1 record and verified route in cryptanalysis; its general binary
Weierstrass model is not yet registered under ICV1. Its ICV1 fields are
therefore `null`. Do not synthesize an ICV1 slug from a translated or
isomorphic model without a recorded isomorphism and registry support.

The [factor-base browser](https://aburan28.github.io/crypto/browser/#factor-bases)
is generated from committed ecbench sessions and lists `FB1` summaries.
It is not a point store. Large bases remain in cryptanalysis's
[`fb-archive`](https://github.com/aburan28/cryptanalysis/tree/main/experiments/fb-archive):
compressed point sets, orbit representatives, or content-addressed shards
with counts and hashes. Only attach an `FB1` to an archive entry after
the curve UID, point-set encoding/digest, quotient rule and actual usable
point count have been reconciled. `fb<B>` cannot be formed from a recipe
with unknown `B`.

[Typed curve-link rules](curve-links/README.md) distinguish twists, same-field
model changes, and field extensions from isogenies. The mirrored curve YAML
keeps explicit `not_enumerated` inventories for links that have not been
constructed or checked. A shared j-invariant is not a base-field point map.

The directed isogeny graph and ordered `IW1` routes live in
cryptanalysis's
[`isogeny_routes.json`](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-candidate-catalog/isogeny_routes.json).
An isogenous neighbor gets its own curve UID. Its conductor, prime-specific
volcano level, map, kernel and target-log transport proof are separate
evidence. An unknown value is YAML `null` with a status, never `0`
or a question mark. A complete IC claim still needs an exact candidate,
frozen workload, verified one-target run, and paired rho run under the
existing measurement contract.

## Walking an isogeny class

**How to run it, including S3 setup:** [`docs/isogeny-walk/README.md`](../../isogeny-walk/README.md).

`src/bin/isogeny_walk.rs` walks the `F_p`-isogeny class of a prime-field
curve (P-256, P-224, or a custom short Weierstrass curve with its order and
generator) over the `ℓ`-isogeny graphs for a set of odd primes `ℓ`, breadth
first, and writes a directory in these formats:

| file | what |
|:--|:--|
| `curves.yaml` | every curve reached, as records of [`curves.yaml`](curves.yaml)'s schema: ICV1 slug and full identity, EC1 alias and UID, the EC1 preimage as `field` and `curve`, traits with statuses, proved volcano levels, a `V<ℓ>L<level>` position alias, and `IW1` route ids |
| `isogeny_routes.json` | cryptanalysis's `isogeny_routes.json` layout: `curve_nodes`, `edges` with kernel polynomials, `routes` with `IW1` ids (one per edge and one root path per curve) |
| `walk.json` | configuration, class invariants (trace, `t² − 4p` and its trial factorisation, twist order, embedding degree, splitting of each `ℓ`), counts, trait distributions |

```bash
cargo build --release --bin isogeny_walk
./target/release/isogeny_walk class  --curve p256 --max-ell 61
./target/release/isogeny_walk walk   --curve p256 --max-ell 61 --max-curves 2000 --out DIR
./target/release/isogeny_walk verify --curve p256 --dir DIR
```

How an edge is found and why it is trusted:

- **Neighbours.** The roots of `Φ_ℓ(j, Y)` in `F_p`.  `Φ_ℓ mod p` is built
  natively from the `q`-expansion of `j`; the system is triangular, and
  every exponent that does not determine a coefficient is checked to vanish.
- **Kernel.** Elkies' construction turns each root into an explicit kernel
  polynomial.
- **Certificate.** `kernel::verify_kernel` accepts it only if it is
  squarefree of degree `(ℓ−1)/2`, divides the `ℓ`-division polynomial, and
  its roots are closed under `[g]` for a generator `g` of `(Z/ℓ)^*/±1`.
  Vélu's codomain must then be `F_p`-isomorphic to the recorded target.
  `verify` replays that certificate from `isogeny_routes.json` alone,
  together with every curve's ICV1 and EC1 identity and an order audit.
- **Order.** For a prime group order the audit proves `#E = n`.
- **Directions.** An `IW1` edge is `h`, `d` or `u` only when both
  endpoint levels are proved:
  - depth 0 (`v_ℓ(t² − 4p) ≤ 1`): every curve is at level 0;
  - depth 1: `ℓ + 1` rational `ℓ`-isogenies means the surface, one means
    the floor.

  Any other edge is `x`.
- **Models and generators.** Each walked curve's model and generator
  follow fixed, recorded rules (`icwalk-canon/v1`, `icwalk-gen/v1`).  The
  generator is **not** the image of the root's generator, so no
  discrete-log transport is established.

Walked curves get computed ICV1 slugs and EC1 identities but are not added
to [`registry.json`](../registry.json); register a slug before citing it in
prose (`AGENTS.md` §11).  The first runs are in
[`research/isogeny_walk_20261004/`](../../../research/isogeny_walk_20261004/README.md).

### Traits and queued detection

Trait detection (`src/cryptanalysis/isogeny_walk/traits.rs`) has two
scopes, because isogenous curves over `F_p` share `p`, `#E` and the trace.

- **Class audits** depend only on `p` and `#E` and run once per walk, on
  the root.  They are `ecc_safety`'s parameter audit, the structural report
  (`p256_structural::curve_structural_report`) and the
  Petit–Kosters–Messeng signals (`pkm_criterion`).  `walk.json#class_audits`
  records them.  The two audits that read a model are re-run on a sample of
  walked curves (`--class-audit-sample`, default 8), and
  `invariance_check.identical_to_root` records whether every verdict
  matched the root's.
- **Curve detectors** (the `Detector` trait) read one curve's model and
  generator.  Each writes one `trait_status` entry in `curves.yaml`.
  - The built-in detectors are `non_singular`, `generator_valid`,
    `a_minus_3_model`, `qr_prefix_64` and `coefficient_bits`.
  - Add a detector to `default_detectors()` and every walk, shard and
    queue job records it.

Detection also runs as queued jobs, using `taskq` (`taskq/README.md`).

```bash
./target/release/isogeny_walk plan --curve p256 --max-ell 61 --max-curves 20000 \
    --commit "$(git rev-parse HEAD)" --shards 16 --out specs/
for f in specs/*.json; do taskq submit --spec "$f"; done
# Fetch each task's artifacts into its own directory, then:
./target/release/isogeny_walk collect --out traits/ shard-dirs/*
```

- **What each job does.** Each spec checks out the pinned commit, builds
  the walker, and rebuilds the walk in its setup step.  The walk is
  deterministic, so no shared storage is needed.  The job then runs
  `isogeny_walk traits --shard i --of N` into `$TASKQ_OUTPUT_DIR`.
- **Ownership.** A shard owns the curves with `index ≡ i (mod N)`.  It
  writes `traits.jsonl`, `metrics.json` (taskq parses it into the run's
  metrics) and, for shard 0, `class_audits.json`.
- **Collecting.** `collect` merges the shards and refuses:
  - a missing or repeated shard;
  - a curve seen twice or not at all;
  - shards built from different walks (the routes file's SHA-256 differs);
  - a `traits.jsonl` that does not match its recorded hash.

  It reports each trait's distribution in `collect.json`.
- **Without a queue.** The same shards run locally in a loop.

### Running on a laptop, offline

Everything except `--store`, `fetch` and `collect --from` runs without a
network.  The walker is plain Rust: no GPU, no Python, no services.

```bash
# Once, while online: fetch the crates so later builds need no network.
git clone https://github.com/aburan28/crypto && cd crypto
cargo fetch
# Offline from here.
cargo build --release --offline --bin isogeny_walk
B=./target/release/isogeny_walk

$B class  --curve p192 --max-ell 61                       # seconds: class invariants
$B walk   --curve p192 --max-ell 61 --max-curves 2000 --out walk-p192
$B verify --curve p192 --dir walk-p192                     # independent replay
for i in 0 1 2 3; do $B traits --curve p192 --dir walk-p192 --shard $i --of 4 --out traits-$i; done
$B collect --out traits-p192 traits-0 traits-1 traits-2 traits-3

# Back online: upload what you made, exactly as `--store` would have.
$B publish --dir walk-p192 --store s3://crypto-autoresearcher/isogeny-walk
$B publish --dir traits-0 --store s3://crypto-autoresearcher/isogeny-walk --run <run id from walk.json>
```

- **Sizing.** Measured on a 14-core M4 Pro, unisolated, so this is a
  practicality note:
  - a 20,000-curve P-256 walk with every odd `ℓ ≤ 61` takes about
    2.5 minutes;
  - P-224 takes about 6.5 minutes (more of its curves are expanded);
  - raw output is about 17 KB per curve (curves.yaml plus
    isogeny_routes.json), so plan disk for `--max-curves`.
- **Fewer cores.** `--threads N` caps them.  Time scales roughly with
  cores.
- **`Φ_ℓ mod p`.** It is rebuilt at each start: about 10 s for `ℓ = 59` at
  256 bits.  `--max-ell 31` keeps startup under a second.
- **Determinism.** The walk is deterministic.  A laptop run and an S3 run
  with the same arguments produce byte-identical files and the same run id,
  so `publish` of a duplicate finds the existing run and uploads nothing.

### Storage in S3

Walk and trait outputs are stored in S3, not on local disk.  The store is
`s3://crypto-autoresearcher/isogeny-walk` (`src/cryptanalysis/isogeny_walk/store.rs`),
laid out on the pattern of PR #1330's campaign contract.

```text
runs/<curve>-<walk key>/attempts/<attempt id>/{curves.yaml,isogeny_routes.json,walk.json}.gz
runs/<curve>-<walk key>/complete.json
runs/<curve>-<walk key>/traits/<N>/shard-<i>/…/complete.json
runs/<curve>-<walk key>/traits/<N>/collected/…/complete.json
```

- **Write-once.** Every object is created with `If-None-Match: *` and
  never overwritten.
- **Hash-checked.** Each upload's SHA-256, as S3 reports it, must equal the
  local hash before the run's `complete.json` is created, also
  write-once.
- **One winner.** That marker lists every object by key, stored SHA-256,
  bytes and uncompressed SHA-256.  If two attempts race, one marker wins and
  the loser's objects are never authoritative.  A run that is already
  complete is not uploaded again.
- **Checked downloads.** `fetch` recomputes both hashes of every object.
- **Transport.** The AWS CLI, with the caller's credentials.

```bash
isogeny_walk walk  --curve p256 --max-ell 61 --max-curves 20000 \
    --store s3://crypto-autoresearcher/isogeny-walk --prune-local --out W   # STORE.json names the run
isogeny_walk fetch --from s3://crypto-autoresearcher/isogeny-walk --run <run id> --out W
isogeny_walk traits --curve p256 --dir W --shard 0 --of 8 \
    --store s3://crypto-autoresearcher/isogeny-walk --run <run id> --out S0
isogeny_walk collect --from s3://crypto-autoresearcher/isogeny-walk --run <run id> --of 8 \
    --store s3://crypto-autoresearcher/isogeny-walk --out T
isogeny_walk plan ... --store s3://crypto-autoresearcher/isogeny-walk --walk-from-store
```

With `plan --store --walk-from-store`, each taskq job fetches the stored
walk instead of rebuilding it and publishes its shard.  The workers then
need AWS credentials that can read and write the prefix.

The bucket is encrypted (SSE-S3) and blocks public access.  Versioning is
off; write-once keys make it unnecessary for these objects.
