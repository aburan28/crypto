# ECC2K-130 status page

Static dashboard for the `ecc2k-130` distinguished-point campaign,
published at <https://aburan28.github.io/crypto/status/>. The Action
overwrites `status.json` and `history.json`; this HTML only renders those
files, and `scripts/site/build.py` lays them out for publishing.

The `Iterations` card is derived in the browser rather than snapshotted:
the walk total behind the point count, one point per `2^25.27` iterations
at `HW(x) <= 34`, shown against the `2^60.9` expected cost of a collision.
See `scripts/rho_status/README.md` for why it is not a published field.

Open `index.html` from the working tree next to the two JSON files and it
renders exactly as published; only the site navigation links resolve solely
on the published site.

See [`scripts/rho_status/README.md`](../../scripts/rho_status/README.md)
for secrets, the walker hop, and what is (not) published, and
[`scripts/site/README.md`](../../scripts/site/README.md) for the published
URL layout.
