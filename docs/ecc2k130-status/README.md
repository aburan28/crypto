# ECC2K-130 status page

Static dashboard for the `ecc2k-130` distinguished-point campaign,
published at <https://aburan28.github.io/crypto/status/>. The Action
overwrites `status.json` and `history.json`; this HTML only renders those
files, and `scripts/site/build.py` lays them out for publishing.

Open `index.html` from the working tree next to the two JSON files and it
renders exactly as published; only the site navigation links resolve solely
on the published site.

See [`scripts/rho_status/README.md`](../../scripts/rho_status/README.md)
for secrets, the walker hop, and what is (not) published, and
[`scripts/site/README.md`](../../scripts/site/README.md) for the published
URL layout.
