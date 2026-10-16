# Complete ordinary isogeny classes

The native `isogeny-class` command builds a resumable, disk-backed multigraph
for an explicit finite set of prime isogeny degrees. It separately reports
exhaustion of that component, completeness of the full vertex class, conductor
strata, map replay and DLP coverage. It never equates a bounded walk with an
entire class.

Vertices are F_p-isomorphism classes with one fixed trace. Each edge represents
one rational cyclic kernel, with an explicit target-model isomorphism. Expanding
the finite degree set adds kernels; changing models or composing maps does not
require treating every resulting formula as a new vertex or generator edge.

Build with `cargo build --release --locked --manifest-path tools/volcano-cli/Cargo.toml`.
The same operations are available as `crypto isogeny class ...` and
`isogeny class ...` after building those binaries from this revision.

```sh
# JSON uses exact integer values, not floating point or hexadecimal strings.
# {"p":...,"a":...,"b":...,"order":...}; order is optional (native SEA).
isogeny-class plan --curve curve.json
isogeny-class build --curve curve.json --out graph.sqlite --steps all
isogeny-class resume --input graph.sqlite --steps all
isogeny-class status --input graph.sqlite
isogeny-class verify --input graph.sqlite
isogeny-class vertex --input graph.sqlite --id 0
```

`plan` computes the ordinary Frobenius discriminant, fundamental discriminant,
conductor of Z[pi], all possible endomorphism orders, their class numbers and
the full vertex population. It certifies the supplied order using point-order
divisors on the curve/twist (or exact enumeration for small p). The source's
endomorphism conductor remains unknown until graph labeling; it is not set to
the Frobenius conductor. The plan emits the exact source model's ICV1 identity.

Automatic degrees include every conductor prime and, when the maximal order
has class number above one, non-inert primes through floor(sqrt(|D_K|)). This
conservative Minkowski bound generates the surface class group. An explicit
`--degrees 2,3,5` explores exactly those degrees, including rational kernels
whose generators are only defined over extensions. Repeated target vertices,
horizontal loops and parallel kernels remain separate edges. The full vertex
population can be reached without enumerating every prime degree.

The supported backend is ordinary prime-field short Weierstrass curves with
`3 < p < 2^61`, `-10000000 <= D_K < -4`, and prime degrees at most 31, different
from p. Larger fundamental discriminants, special j, supersingular curves,
extension fields and larger degrees explicitly require another backend. The
general division-polynomial kernel search can be expensive for larger degrees;
an expansion budget counts completed vertices, not internal factorization work.
This is not a promise of tractability for arbitrary curves or degree sets.

`build` refuses an existing output. One expansion commits all outgoing kernels
and newly discovered vertices atomically. SQLite uses full synchronous writes;
an interrupted expansion rolls back, and `resume` continues the frontier. Keep
the SQLite file and any journal together while a writer is active. The database
freezes source inputs, selected degrees and the implementation-source hash.
Resume refuses a changed source hash; verification may replay an earlier source.
Only one writer should operate on a database at a time. Plan/status read-outs
are not certificates against simultaneous or malicious edits.

The default **10,000 expansions is a partial checkpoint**. Exit 3 indicates
partial graph/class coverage. An exhausted frontier with fewer vertices than
the exact class number is `COMPLETE_COMPONENT_ONLY`, also exit 3. Use
`--steps all` for an unbounded frontier run. When the exact class is complete,
build/resume automatically performs replay and stratum labeling before exit 0.
The separate `verify` checks SQLite integrity, root reachability, all stored
models, regenerated outgoing kernel sets, exact rational-map polynomial
identities, target model scalings and each conductor population. Replay shares
native arithmetic with construction: independent implementation replay is
explicitly **NOT_RUN**. The small controls also compare every vertex against an
independent exhaustive census of all nonsingular models.

SQLite is the lossless graph export: `nodes(id,a,b,j,expanded,conductor)` and
`edges(source,ell,kernel,target,u)`, with the full plan and replay receipt in
`metadata`. `kernel` lists ascending coefficients. Reconstruct the normalized
quotient with the existing Kohel backend; the stored isomorphism is
`(x,y) -> (u^2*x,u^3*y)`. `vertex` emits the resulting rational maps, exact target
models and identities, including exceptional-input handling. All integer
coordinates are canonical residues. No canonical model is silently substituted.

```sh
isogeny-class solve --input graph.sqlite --vertex 0 \
  --px X --py Y --qx XQ --qy YQ --order R --max-baby 1000000
```

`solve` verifies P's exact order, runs BSGS in that cyclic subgroup and checks
`[k]P=Q`. It solves one explicit target. It does not run a graph-wide experiment,
Pollard rho, PH or index calculus; graph construction reports that workload as
not run. The relation-yield and min/max/mean campaign remains a separate requested
deliverable. Lookup/entry counts from this solver are diagnostics; total cost
and timing are null, and no method ranking is inferred.

Use `plan` before provisioning storage. Its 24-byte-per-vertex number describes
only a fixed three-u64 a,b,j layout, not a SQLite capacity estimate or an
information-theoretic lower bound. Compressed/generative encodings can be smaller. Edges,
indexes, journals, verification and cryptanalysis records require additional
space. Large outputs belong in the repository artifact store, with hashes,
byte counts and provenance, not Git. Five fixed user cases and their coverage
are tracked in [the study](../../research/full_volcano_classes_20261016/PROTOCOL.md).
