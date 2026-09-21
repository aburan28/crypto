# ECC2K-130 status page

Static dashboard for the `ecc2k-130` distinguished-point campaign,
published at <https://aburan28.github.io/crypto/status/>. The Action
overwrites `status.json` and `history.json`; this HTML only renders those
files, and `scripts/site/build.py` lays them out for publishing.

The operations panel above the counts is the walk total the walkers
checkpointed, as `2^n`, taken from the snapshot's `work.iterations`. Under
it is the **collision-odds panel**: rho has no intermediate progress, so
the page draws the one thing that accumulates, the probability that a
collision has already happened. After `W` iterations that is the birthday
bound `1 - exp(-pi W^2 / 4 E^2)` with its mean set to the expected cost
`E = 2^60.9`, so the curve's mean is the same expectation the ETA counts
to; the page says it is conditional on the walk behaving as a random map.
The curve runs from now to the date the odds reach 99% at the measured
rate, with the fleet doubled drawn dashed, and four tiles give the odds
over the next 30 days (conditional on no collision so far), the dates the
odds pass one half and nine tenths, and the share of the expected work
walked, as a number rather than a bar: that share was `2^-21` when the
page was first published, and a bar filled from it is empty while one
filled from the exponents lies. With no measured rate the curve is drawn
against work instead of dates and the date tiles say why they are missing.

Beneath the bars are the two numbers that depend on time. **Iterations per
second** is measured, not projected: that same total's increase between two
published snapshots, divided by the time between them, with the span it
used printed beside it. The **ETA** is the work still expected at exactly
that rate. A snapshot carrying no iteration total falls back to the point
count at one point per `2^25.27` iterations at `HW(x) <= 34`, which reads
about six times low for this campaign and is labelled on the page as the
fallback; a campaign with one total so far says the rate arrives with the
next snapshot rather than showing a number it cannot measure yet. See
`scripts/rho_status/README.md` for where the total comes from and why the
derivation is only a fallback.

**What the walk builds**, between the cumulative chart and the worker
table, is `walk-forest.svg`: real ECC2K-130 walks the campaign client made
on the challenge curve, replayed with the client's own kernel, sampled along
their length and drawn, generated from the hashed trails in `walk-forest/`
by `scripts/site/walk_forest.py`. Only walks short enough to draw are shown,
the caption says so and states the counts, and the orbits are named by hash
because a distinguished point's orbit is its key and the page publishes
counts only. `walk-forest/README.md` has the exact commands, including how
to point the same pipeline at the fleet's corpus, and
`scripts/site/test_build.py` fails if the SVG, the trails, their recorded
endpoints, the scalar reference's replay and the caption disagree.
`scripts/site/build.py` copies the SVG next to the page so the relative
`src` resolves both in the working tree and once published.

The **Contribute compute** section under the worker table is static: it
carries the [cairn](https://github.com/aburan28/cairn) download link
(`releases/latest`, plus the one-line installer the cairn README
documents) for the paid piecework path, and the `ecc2k130/` client
commands for the unpaid one. It states plainly that no ECC2K-130
objective is posted on cairn yet — the binary-field checker
`GF(2^131)` needs is not in cairn's `examples/certicom-ecdlp/` — so the
download is an invitation to be ready and to work the live rungs, not a
claim that points on this campaign are payable today. If that changes,
this section is what has to change with it.

Open `index.html` from the working tree next to the two JSON files and it
renders exactly as published; only the site navigation links resolve solely
on the published site.

`how.html` is the readable walk: the same iteration as the GPU client, on
`GF(2^23)`, with `rho-toy.js` recovering a planted `k` in the browser. The
Python original is `ecc2k130/examples/rho_toy.py`.

The dashboard fetches the live S3 `status.json` and the Pages copy together
and fills missing fields (walk rate, walking slots, per-worker rows) from
whichever document has them, then measures the rate from `history.json`
when neither snapshot carries one. `published_at` is kept from whichever
source wrote more recently: a healthy Pages job on a frozen ingest feed
must not inherit the feed's old stamp, or the banner blames the publisher
for a dead feed. Past the stale threshold, `walking_slots` and `walk_rate`
are cleared so a frozen document cannot read as "walking now".

See [`scripts/rho_status/README.md`](../../scripts/rho_status/README.md)
for secrets, the walker hop, and what is (not) published, and
[`scripts/site/README.md`](../../scripts/site/README.md) for the published
URL layout.
