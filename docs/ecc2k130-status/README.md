# ECC2K-130 status page

Static dashboard for the `ecc2k-130` distinguished-point campaign,
published at <https://aburan28.github.io/crypto/status/>. The Action
overwrites `status.json` and `history.json`; this HTML only renders those
files, and `scripts/site/build.py` lays them out for publishing.

The operations panel above the counts is the walk total the walkers
checkpointed, as `2^n`, taken from the snapshot's `work.iterations`. Two
bars sit under it. The first is **linear in work** — the share of the
`2^60.9` expected cost of a collision — so early in a campaign it is empty
and its caption says so; nothing pads the fill to a visible sliver. The
second is the exponent on a **log scale**, legible but labelled on the page
as not being progress, because each bit of it is a doubling of the work.

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
table, is a static figure: `walk-forest.jpg`, a drawing of the forest that
deterministic walks grow when their trails merge, with a few of the meetings
picked out in colour. It is an illustration of the mechanism, and the page
says so beside it; it is not rendered from this campaign's table, whose walks
never leave the private store. `scripts/site/build.py` copies it next to the
page so the relative `src` resolves both in the working tree and once
published.

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

See [`scripts/rho_status/README.md`](../../scripts/rho_status/README.md)
for secrets, the walker hop, and what is (not) published, and
[`scripts/site/README.md`](../../scripts/site/README.md) for the published
URL layout.
