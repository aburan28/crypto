# The walk-forest figure

`../walk-forest.svg`, the drawing under **What the walk builds** on the
status page, is generated from the files in this directory: real ECC2K-130
walks, made by the campaign client on the challenge curve, replayed with the
client's own kernel and drawn. Nothing in it is drawn by hand: every node is
an orbit one of those walks passed through, every edge a fixed number of
iterations along the trail, and every filled node the distinguished point the
walk reported.

## What the data is

| File | What it is |
|:--|:--|
| `trails.txt` | The drawn walks, from `trailforest --sample`: for each, its step count and the orbits it stood on every `every` steps plus the one it ended on, each orbit replaced by the first 16 hex digits of the SHA-256 of its canonical representative. The header carries the curve, cutoff, stride, cap, how many records were tried, how many finished inside the cap and were drawn, and how many of those were checked against their record (all of them). |
| `forest.hashes` | The name of the orbit each drawn walk's record says it ended on, one line per walk; `walk_forest.py` refuses to draw a trail that does not end there. |
| `reference-check.txt` | The scalar reference implementation's replay of the first records, hashed and sampled the same way with `walk_forest.py --hash-trails`. The kernel's trails for those walks must match it orbit for orbit, and the test suite checks that they do. |
| `gf2-23/` | The tool's regression set on the `GF(2^23)` test curve, where orbits can be committed in the clear: the client's own records for a run (`client-run1.bin`), the reference walk of that run's seed schedule (`trails.txt`) and its endpoints (`forest.bin`). Not drawn. |

The walks are the client's: `ecc2k130-cpu --curve 131` on the Certicom
ECC2K-130 parameters, distinguishing at normal-basis weight `w <= 34`, the
cutoff in `ecc2k130/aws/campaign.json`, writing the same 32-byte
`(seed, canonical orbit)` records the fleet uploads. A walk on this curve is
expected to take `2^25.27` iterations to reach a distinguished point, which is
far too long a trail to draw, so the figure draws the walks that finished
early: every record's seed is replayed for at most `cap` steps and the ones
that reach their distinguished point inside that are kept. That is a real
sample of real walks, biased towards short ones by construction, and the
caption says so.

## Why the orbits are hashed

The campaign publishes counts and nothing else: not point keys, not walk
coefficients, not seeds (`scripts/rho_status/snapshot.py` refuses to). A
distinguished point's orbit *is* its key, so the trails here name orbits by a
64-bit prefix of a hash of a 131-bit value, which identifies an orbit for
drawing and testing without revealing it. The records themselves stay off the
repository; `forest.hashes` carries only the names. The walks drawn are also
under a run-id far outside the fleet's slot registry and were never uploaded,
so they are not part of the campaign's corpus. The same pipeline, pointed at
the fleet's corpus, would draw the fleet's walks (below).

## Regenerating

From the repository root:

```bash
# 1. real walks: the client on the challenge curve, ten minutes on four cores,
#    records to a file that stays out of the repository
make -C ecc2k130 cpu trailforest
ecc2k130/ecc2k130-cpu --curve 131 --dp-weight 34 --run-id 65000 --threads 4 \
    --steps 1024 --launches 400 --verify 0 --dp-file /tmp/real131.bin

# 2. the trails behind the first 256 records, replayed with the client's
#    kernel, kept when they finish within the cap, sampled every 512 steps,
#    named by hash, each checked against its record
ecc2k130/build/trailforest --curve 131 --dp-weight 34 --sample \
    --corpus /tmp/real131.bin --max 256 --cap 65536 --every 512 \
    --hashes-out docs/ecc2k130-status/walk-forest/forest.hashes \
    > docs/ecc2k130-status/walk-forest/trails.txt

# 3. the scalar reference over the first records, hashed the same way
ecc2k130/build/trailforest --curve 131 --dp-weight 34 \
    --corpus /tmp/real131.bin --max 8 > /tmp/ref131.txt
python3 scripts/site/walk_forest.py --trails /tmp/ref131.txt --hash-trails 512 \
    --out docs/ecc2k130-status/walk-forest/reference-check.txt

# 4. the drawing, and the same forest as a graph for the page's explorer
python3 scripts/site/walk_forest.py \
    --trails docs/ecc2k130-status/walk-forest/trails.txt \
    --corpus docs/ecc2k130-status/walk-forest/forest.hashes \
    --out docs/ecc2k130-status/walk-forest.svg
python3 scripts/site/walk_forest.py \
    --trails docs/ecc2k130-status/walk-forest/trails.txt \
    --corpus docs/ecc2k130-status/walk-forest/forest.hashes \
    --json --title "ECC2K-130, real walks" --out docs/ecc2k130-status/walk-forest.json
python3 scripts/site/walk_forest.py \
    --trails docs/ecc2k130-status/walk-forest/gf2-23/trails.txt \
    --corpus docs/ecc2k130-status/walk-forest/gf2-23/forest.bin \
    --json --title "GF(2^23) test curve" --out docs/ecc2k130-status/walk-forest-gf2-23.json
```

## The explorer

`../walk-forest.js` loads `../walk-forest.json` (the real curve) and
`../walk-forest-gf2-23.json` (the test curve, where trails meet) and draws
them on a canvas in place of the figure: drag to pan, wheel or pinch to
zoom, hover a node for the walk it is on and how far along it sits, click a
node to light every path through it to its distinguished point, and play the
walks to watch each set off along its trail. The graphs carry the static
figure's own layout, so the two agree and the page computes no layout; a
graph is byte-identical to what its trails export to, and the tests check
that. Without JavaScript, or if a graph fails to load, the figure stays.

Steps 2 to 4 are deterministic given the records. The reference in step 3
walks a few thousand steps a second, so it is only run over the first
records, which are the shortest; `--max` bounds it. `scripts/site/test_build.py`
fails if the committed SVG is not what the committed trails render to, if a
trail does not end on its recorded endpoint, if the reference's trails
disagree with the kernel's, or if the caption's counts are not the ones in
the trails file. Changing the cap, stride or record count changes the figure;
update the caption in `../index.html` with it, and the test will say which
number is stale.

## Drawing the fleet's own walks

The fleet's records are the same format, under `dp/slot-N/*.bin` in the
campaign bucket (`ecc2k130/aws/README.md`). From a host with access to it,
concatenate any set of delta files into one corpus and run steps 2 to 4 on
it; nothing else changes, and `trailforest` checks every drawn trail against
its record the same way. At the fleet's live cutoff of weight 32 a walk is
expected to take `2^27.9` iterations, so `--cap` selects the drawable ones
just as it does here.

## The `GF(2^23)` set

`gf2-23/` was the figure's first source and stays as the tool's test: a curve
small enough that the client's search solves inside its first launch, so
`trailforest --generate` walks the run's seed schedule without stopping and
`--check-corpus` holds the result to the client's records
(`gf2-23/client-run1.bin`). To regenerate it:

```bash
make -C ecc2k130 ecc2k130-cpu BATCH=1
ecc2k130/ecc2k130-cpu --curve 23 --instance 0 --threads 1 --steps 64 --launches 1 \
    --verify 4 --dp-weight 7 --run-id 1 --dp-file docs/ecc2k130-status/walk-forest/gf2-23/client-run1.bin
ecc2k130/build/trailforest --curve 23 --instance 0 --dp-weight 7 --generate --run-id 1 --walks 64 \
    --check-corpus docs/ecc2k130-status/walk-forest/gf2-23/client-run1.bin \
    --corpus-out docs/ecc2k130-status/walk-forest/gf2-23/forest.bin \
    > docs/ecc2k130-status/walk-forest/gf2-23/trails.txt
```
