# isogeny

`isogeny` is one executable for catalogue-wide prime-degree screening, bounded
isogeny construction, and independent exact verification. Its algorithms and standard
curve data are compiled into the binary. Running it requires no checkout, Cargo,
companion executable, or external registry file.

## Install the released executable

Download the `isogeny` archive matching your platform from
[the crypto releases](https://github.com/aburan28/crypto/releases), extract it, and put
the `isogeny` executable on your PATH. The release includes SHA256SUMS for the archives
and BUILD.json with the source revision and embedded-standard-data hash.

The release workflow builds Linux glibc, Linux musl, Apple silicon, and Windows
executables. Pull requests build and test the archives; publication runs from main
after all required release jobs succeed. ReleaseMe manages publication.

```sh
isogeny --help
isogeny --version
isogeny curves

# Construct and independently verify the recorded P224 case.
isogeny search --curve p224 --ell 1471 --out ./p224-1471

# The recorded P192 case; its independent replay took about 11 minutes on an M4 Pro.
isogeny search --curve p192 --ell 10453 --out ./p192-10453

# Other registered prime curves, selected by name, alias or ICV1 slug.
isogeny search --curve P-384 --ell 19 --out ./p384-19
isogeny search --curve P-521 --ell 7 --out ./p521-7

# Recheck a saved run; print a certificate without changing its input files.
isogeny verify --input ./p224-1471

# Screen for small-order Frobenius eigenvalues above degree 1009.
isogeny screen --curve p192 --from 1010 --to 65537 --max-order 8

# Find distinct maps by advancing only through independently verified codomains.
isogeny walk --curve p192 --ell 2341 --count 50 --out ./p192-large-walk
isogeny walk --curve p224 --ell 1117 --count 50 --out ./p224-large-walk
isogeny verify-walk --input ./p224-large-walk
```

Search requires a fresh output directory and uses seed 1 and the published curve
order automatically. Construction has a default 1800-second limit and an 8 GiB
resident-memory cap. Independent replay follows construction; its cost is separate
from that construction timeout and memory cap.

`walk` requires a fresh directory and a prime degree above 1009. The native kernel
route accepts P192/P224 eigenvalue orders through 16 and requests of 1 through
10,000 maps. Each step retains its exact
map, resource receipt, independently transported subgroup generator, and linked
hashes. `walk.json` records verified progress after every step. The independent
`verify-walk` command replays the entire chain from the embedded public source,
recomputing every destination and generator without modifying its inputs. Distinct
maps may share the same prime degree; report map count and degree count separately.
One eigenline per step retains coverage `PARTIAL`. Interrupted output is evidence;
the current walk command does not resume into an existing directory.

`curves` exposes the complete pinned repository catalogue, including NIST/SEC,
Brainpool, SM2, GOST/ANSSI, Montgomery/Edwards, binary, Koblitz, subfield and extension
models. Unresolved standards-source records remain visible. This is the inventory
of those pinned inputs rather than a worldwide census.

Each model carries separate screening, construction and independent-replay capabilities.
Known-order curves across the listed field families can be screened from their exact
Frobenius trace. Construction supports registered short-Weierstrass prime models through
640 bits, including P521. Replay requires a complete registered public subgroup generator.
Montgomery/Edwards construction still requires checked coordinate-conversion adapters;
binary/extension maps require their own construction and verification adapters. Their
catalogue and screening entries remain available.

Names resolve to exact catalogue models. Ambiguous names require an ICV1 slug.
Subgroup orders and cofactors are preserved, and the degree must be coprime to the
registered subgroup order. Independent replay uses exact-size field specializations,
binary extended-GCD inversion, and the retained exact polynomial checks.
Execution lookup uses a build-generated compact index while `curves` returns the
complete embedded inventory. Receipt/executable SHA-256 uses a pinned native backend
with runtime CPU detection and a portable fallback checked against independent hashes.

`--method auto` uses extension-field kernel construction for P192/P224 at prime
degrees above 1009 when the smaller Frobenius eigenvalue order is at most 16.
The recorded cases include P224/1117, P224/1471, P192/2341 and P192/10453.
The constructor checks the extension field, torsion order, Frobenius eigenline
and descent of the kernel; independent replay then checks the resulting map.
These runs construct one Frobenius eigenline of two and report
`degree_coverage: PARTIAL`. Other supported split prime degrees use full
modular-polynomial enumeration. `--method kernel` explicitly selects the checked
kernel route where supported; `--method modpoly` selects full enumeration.
Search and screening currently accept odd prime degrees up to 1000000. These are
implementation bounds; a bounded construction timeout does not establish nonexistence.

```sh
# Run full enumeration with a chosen construction budget.
isogeny search --curve p192 --ell 1021 --method modpoly \
  --timeout 1800 --out ./p192-1021

# Omit independent replay explicitly; the result is labeled CONSTRUCTED/NOT_RUN.
isogeny search --curve p224 --ell 1471 --construct-only --out ./p224-construction
```

Each command prints JSON on stdout. Stage messages go to stderr. Exit 0 denotes a
passing result or an explicitly requested construction-only result; 1 denotes a
failed construction or replay, 2 a usage error, and 3 a timeout/resource limit.
Construction-only output has status `CONSTRUCTED`, with independent replay `NOT_RUN`.
It carries no independent certificate.

Search retains `search.json`, raw map output and stderr, a resource/exit receipt,
and, after passing independent replay, `replay.json` and `curves.json`. The exact
verifier checks the kernel, codomain, rational-map identity, standard subgroup,
and 20 fresh public scalar-transport cases. Failed runs retain their construction
records. The tool never promotes its results into the canonical catalogue automatically.

## Maintainer build and checks

The repository builds this standalone package in the release workflow. These are
development commands; released-tool users run the executable directly.

```sh
cargo build --release --locked --manifest-path tools/isogeny-cli/Cargo.toml
cargo test --release --locked --manifest-path tools/isogeny-cli/Cargo.toml
```

The package links the existing standalone algorithms and the independently validated
field, polynomial, and kernel verification implementations. The dated research sources
and evidence remain frozen. Dependencies are pinned to the independent replay study's
lockfile versions. Runtime curve data and version/source provenance are embedded.

## Largest prime-degree challenge

`isogeny challenge policy` prints the pinned **Prime Degree Frontier** rules
and all registered prime-field models. Each exact source curve ranks separately;
unsupported models stay visible. `challenge plan --curve all` creates a
candidate-only ladder. `challenge run`, `challenge judge` and
`challenge leaderboard` construct or independently replay explicit maps in
fresh directories. Saved scores are not accepted as submissions. See
[challenge rules and commands](../../docs/challenges/isogeny-prime-degree/README.md).
The new command is included by the existing ReleaseMe build; actual new-asset
publication remains a separate release gate.
