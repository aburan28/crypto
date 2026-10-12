# WeakCurves criterion filters and isogeny graph

The atlas now joins **1,090 exact curve nodes**, **674 finite-field isogeny
classes**, and **232 independently verified directed isogeny edges**. The graph
contains all 920 leaderboard models and 170 additional models appearing in
the P-192/P-224 walk evidence. Seventy-three indexed classes contain at least
two models. The class and map relationships are deliberately different data.

## Relationship rules

For elliptic curves over a finite field of size `q=p^n`, equal Frobenius trace
means equal point count and, by [Tate's finite-field isogeny theorem](https://pazuki.perso.math.cnrs.fr/index_fichiers/Tate66.pdf),
membership in the same `F_q`-isogeny class. The exporter assigns a class ID
`FQI1:` followed by 16 hex digits of SHA-256 over the UTF-8 string
`p:n:trace`. This is an **abstract existence** relationship. The unique finite
field of each size may have several polynomial representations in the corpus;
an explicit route between their recorded coordinates needs a field
isomorphism as well as a curve map.

Bright graph edges require a named, verified map record. There are 132 edges
from the [P-192 radius-two walk](../isogeny_walk_p192_20261004/runs/p192-ell61-radius2/isogeny_routes.json),
bound to its [independent verification receipt](../isogeny_walk_p192_20261004/runs/p192-ell61-radius2/VERIFY.json),
and 100 from the [P-192/P-224 large-prime map registry](../../docs/curves/large-prime-walks-20261010.json).
The former records degree, kernel polynomial, Vélu codomain and endpoint UIDs;
the latter records degree, certificate digest, exact rational-map check,
subgroup check and public scalar-transport checks. The exporter checks the
receipt hash, every endpoint and class invariant, map status and edge-ID
uniqueness before writing the graph. It also opens all 100 large-prime map
artifacts and checks their SHA-256 digests against the registry certificates;
the graph exposes a direct link for each map. These records give map evidence, not an
automatic DLP advantage. The map network is a bounded walk, so an absent edge
does not imply no isogeny exists.

## Criterion filter semantics

The filter exposes five stages. `condition` is an exact selected-group
mathematical condition. `geometry` is an exact covering-model label. `map`
means an explicit isogeny map is indexed. `candidate` needs additional transfer
or cost evidence. `screen` is a necessary or sufficient structural lead that
does not alone certify a weaker DLP. A curve matches only when the export has a
recorded positive entry at the requested stage. Lack of a tag is unassessed,
not a negative result. The 16-criterion registry describes recognition and
witness requirements even where no atlas model currently meets a check.

| Check | Current automated match | Further evidence for a DLP finding |
| --- | --- | --- |
| MOV and Frey–Rück | Exact selected-order embedding degree, where recorded | Nondegenerate pairing image, target transport, destination-field cost |
| Isogeny transfer | Endpoint of one of 232 verified directed maps | Selected subgroup transport and an independently assessed cheaper destination |
| GLV | Prime-field `j=0, p≡1 mod 3` or `j=1728, p≡1 mod 4` automorphism screen | Exact endomorphism action on `r`, orbit/decomposition method and matched cost |
| GLS | Published quadratic-extension twist construction, currently no exact atlas match | Exact twist/Frobenius map, subgroup eigenvalue and matched cost |
| Cover and decomposition | Exact local covering-model label | Subgroup map, natural relation yield, rank, descent and total cost |
| Cubic norm-one | Exact covering label or positive necessary conductor-depth screen | For screen-only rows, construct and verify the cover; then selected-group route and cost |

GLV and GLS were designed to accelerate scalar multiplication. An
endomorphism can also change an attack's constant factors, but the structural
screen is not a security ranking. The site therefore exposes them as method
leads with exact evidence stages rather than calling every matching curve
weak. The filter can show zero matches for a documented method; GLS currently
has no exact atlas-model screen.

## Reproduction

From the repository root:

```sh
python3 research/weak_curve_leaderboard_20261012/build_isogeny_graph.py --check
cmp research/weak_curve_leaderboard_20261012/isogeny_graph.json website/weakcurves/dist/isogeny_graph.json
node --check website/weakcurves/dist/app.js
node --check website/weakcurves/dist/graph.js
```

Run the builder without `--check` after a frozen source changes, then copy its
JSON to the site. `source_sha256` in the graph records the exact snapshot,
route and verification files read by the builder. The browser graph uses the
same JSON as the filter and details panel; it has no separate inference
engine. The graph's class spokes are not counted as certified map edges.
