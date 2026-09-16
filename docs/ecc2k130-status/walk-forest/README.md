# The walk-forest figure

`../walk-forest.svg`, the drawing under **What the walk builds** on the
status page, is generated from the files in this directory. Nothing in it is
drawn by hand: every node is an orbit one of the campaign client's walks
passed through, every edge one iteration, and every filled node a
distinguished point.

## Why it is drawn on GF(2^23)

On the challenge curve a single trail is about `2^27.9` iterations long and
the campaign records only its endpoint, so the trails themselves cannot be
drawn at that scale. The figure is drawn instead on the smallest curve the
client walks, `GF(2^23)` in the same permuted type-II optimal normal basis,
with the client's own seed schedule (`eccSeedFor(run_id, lane)`), iteration
function `R -> R + sigma^j(R)` and distinguishing test (normal-basis weight
of `x` at most `w`). What changes with the field size is how long a trail is
and how often trails meet; what the picture shows is the mechanism, which is
the same.

## Files

| File | What it is |
|:--|:--|
| `client-run1.bin` | What `ecc2k130-cpu` itself wrote to `--dp-file` for run-id 1 on this instance, 32-byte records `(seed, canonical orbit)`, before its search solved the instance and stopped. |
| `trails.txt` | Every orbit each drawn walk passed through, from `trailforest --generate`: lanes 0 to 63 of run-id 1, first walk each, held to `client-run1.bin` wherever the two overlap (`checked` in the header). |
| `forest.bin` | The endpoints of the drawn walks in the client's record format, written by the same run; `walk_forest.py` refuses to draw a trail that does not end on its record. |
| `README.md` | This file. |

A search on this curve solves inside its first launch, so a client corpus
holds only the handful of walks reported before that, and the shortest ones
first. `trailforest --generate` walks the same seeds with the reference
implementation the client verifies its reports against (`--verify`), but
keeps every lane's first walk instead of stopping at the first collision.
`--check-corpus` then holds the generated trails to the client's real
records: every record for a drawn lane must end on the orbit the replay
reached, or the tool fails instead of printing.

## Regenerating

From the repository root:

```bash
# 1. the client's own records, 512 walks (one slot per thread) on GF(2^23)
make -C ecc2k130 ecc2k130-cpu BATCH=1
ecc2k130/ecc2k130-cpu --curve 23 --instance 0 --threads 1 --steps 64 --launches 1 \
    --verify 4 --dp-weight 7 --run-id 1 \
    --dp-file docs/ecc2k130-status/walk-forest/client-run1.bin

# 2. the trails behind lanes 0..63 of that run, checked against those records
make -C ecc2k130 trailforest
ecc2k130/build/trailforest --curve 23 --instance 0 --dp-weight 7 \
    --generate --run-id 1 --walks 64 \
    --check-corpus docs/ecc2k130-status/walk-forest/client-run1.bin \
    --corpus-out docs/ecc2k130-status/walk-forest/forest.bin \
    > docs/ecc2k130-status/walk-forest/trails.txt

# 3. the drawing
python3 scripts/site/walk_forest.py \
    --trails docs/ecc2k130-status/walk-forest/trails.txt \
    --corpus docs/ecc2k130-status/walk-forest/forest.bin \
    --out docs/ecc2k130-status/walk-forest.svg
```

Step 1 is only needed to refresh the client's records; the default client
build runs 16384 walks, which on a group this small collide at their start
points before a step is taken, hence `BATCH=1`. Steps 2 and 3 are
deterministic: the same trails always render to the same bytes, and
`scripts/site/test_build.py` fails if the committed SVG is not what the
committed trails render to, if a trail does not end on its corpus record, if
a client record disagrees with a drawn trail, or if the caption's counts
are not the ones in the trails file. Changing the run-id, the number of
walks or the cutoff changes the figure; update the caption in
`../index.html` with it, and the test will say which number is stale.
