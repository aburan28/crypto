# Seven-constructor kernel baseline

This extends the requested inventory-wide baseline campaign beyond the eight
BMSS constructors. It changes the native measurement adapter only; the seven
constructor algorithms remain unchanged. Execution and per-CPU receipts are
still required before this series is measured.

```sh
isogeny-perf kernels --out FRESH.jsonl --revision FULL_SOURCE_SHA \
  --bits 192,224,256,511 --ells 7,31,101
isogeny-perf compare baseline.jsonl candidate.jsonl FRESH_COMPARISON.json
```

| Methods | Supplied input | Work counted in construction |
|---|---|---|
| `velu_from_point`, `velu_cyclic`, `sqrt_velu` | Affine kernel generator | Multiples, constructor and owned result allocation |
| `velu_xonly`, `velu_xonly_fast`, `sqrt_velu_fast` | Kernel generator abscissa | Multiples, constructor and owned result allocation |
| `kohel_given_kernel` | Complete kernel polynomial | Rational-map construction and owned result allocation; deriving the polynomial excluded |

Field and curve generation, finding the torsion generator, reference
construction, input/output fingerprints and verification are outside every
constructor interval. A record retains the full public generator, reference
kernel, complete probe list, input boundary and reference-map digest. The
presence of extra provenance data does not make that data an argument used by
the algorithm. Supplied-kernel and supplied-generator rows have different
boundaries; a ratio between their constructor counts is not an end-to-end
algorithm comparison.

Each field-width/degree fixture has independent Kohel and point-enumerating
Vélu references. Their degree, codomain and full transported points must agree
at `2*ell+1` distinct non-kernel abscissae. Every counted, native and timed
constructor invocation must reproduce that fingerprint. The adapter serializes
the resulting degree, codomain and all those transported points canonically.
These are deterministic fixture fingerprints, distinct from the BMSS suite's
full rational-map coefficient comparison. Do not label the kernel series a
new symbolic map proof or a named-NIST search.

Counts use the same generic Field model as the BMSS series, with inversion
internals opaque. Counts are saved before verification touches the counted
field. Timing uses the original native backend, includes the owned evaluation
representation and excludes its later destruction. Native wall samples remain
diagnostic unless a separate admitted isolation receipt qualifies the entire
run. The local shared Apple host therefore records counts without wall samples.

The complete first matrix requests 84 rows: seven methods, four field widths,
three prime degrees. Missing methods, cells, output fingerprints, source or
binary identities, count models and terminal PASS receipts are rejected by
the native reader. The old BMSS `campaign` command still requests its own
64-case degree-7/31 matrix; it must not silently substitute kernel cases.

Before optimizing any of these methods, retain this implementation's source
and executable plus complete per-CPU baselines. Subsequent changes need the
same exact fixtures and fingerprints, five A/A pairs, five interleaved A/B
pairs, and the native isolation receipts for every admitted wall comparison.
