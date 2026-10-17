# Complete graph explorer

[Download the self-contained HTML](https://drive.google.com/file/d/13AQHFuLnePIfbbS9XbaXNueJY-JMyKCE/view).
It contains all 24,823 vertices, 51 conductor strata and 69,300 directed
degree-2/3 kernel maps of the three additional graphs. Download and open it in
a browser; Drive's preview does not execute the application. No external data
or JavaScript dependency is fetched. Graph enumeration is separate from the
99-vertex dev-build measurement checkpoint, whose vertices have orange rings.

Choose a class, filter a conductor, inspect a vertex ID, pan or zoom. The side
panel retains every outgoing kernel map, including parallel maps and self-loops;
double-click a map button to visit its target. An overview can overlay many
edges or vertices at the same pixel. This is a layout limitation, not data
sampling. Conductor columns are a presentation order, not a claim that a mixed
degree graph is one simple descending volcano. Missing benchmark coverage is
explicit. Curve coefficients, j-invariants, scaling parameters and kernel
strings remain exact beyond JavaScript's safe-integer range.

Completeness is for the declared finite prime-degree sets and Fp-isomorphism
vertices. This does not complete the original A–E request, cover other prime
degrees or enumerate infinitely many isogeny compositions. Source graph hashes
and recorded native verification receipts are embedded. Export checks full
expansion, exact planned vertex/stratum populations, endpoints and degree
membership; it does not replace mathematical map replay or independent replay.

The native read-only exporter is [`export.rs`](export.rs). Link it using matching
rlibs from a built checkout, supplying the explicit paths as with the execution
helpers (do not guess among incompatible cached dependencies):

```sh
rustc --edition=2021 -C strip=symbols export.rs \
  --extern "crypto_lib=$crypto_rlib" --extern "serde_json=$serde_rlib" \
  --extern "rusqlite=$sqlite_rlib" -L "dependency=$deps" -o graph-export
./graph-export /complete/graphs /checkpoint/inventory.json template.html /fresh/graphs.html
node test-ui.mjs /fresh/graphs.html
```

The Node presentation test checks all three populations and graph hashes,
all 51 conductor filters, exact large-number strings, root/floor inspectors,
outgoing map lists, deep links, class switches, ID bounds and checkpoint labels.
The actual Chromium checks also pass with no page errors, and the desktop,
floor and mobile captures in `evidence/` were visually inspected. A class-switch
bug was caught and fixed before delivery; its failed log is retained. Dense
conductor labels are thinned at small zoom levels while every stratum remains
available in the selector. The original blocked file-navigation attempts did
not produce screenshots; managed Chromium required loopback HTTP instead.

To repeat visual capture where `file://` is blocked, serve only a directory
holding the generated artifact, then use the presentation-only capture tool:

```sh
python3 -m http.server 8768 --bind 127.0.0.1 --directory /fresh
python3 capture.py http://127.0.0.1:8768/graphs.html /fresh/screenshots
```

Build/export output and the delivery hash are in `evidence/export.json` and
[`delivery.json`](delivery.json). Remote name, size and parent were checked;
remote whole-file hashing was not performed. The frozen measurement archives
were not rewritten to include this supplemental viewer.
