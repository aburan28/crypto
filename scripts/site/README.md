# GitHub Pages site

`build.py` assembles everything published at
<https://aburan28.github.io/crypto/>. The hourly
[`ecc2k130-status`](../../.github/workflows/ecc2k130-status.yml) workflow runs
it after refreshing the campaign snapshot and uploads the result as the Pages
artifact.

```bash
python3 scripts/site/build.py --out _site   # assemble
python3 scripts/site/test_build.py          # offline tests
python3 -m http.server --directory _site 8000
```

The build only copies and lays out files that already exist in the repository,
so no published page is ever the only copy of a figure.

## What is published where

| URL | Source |
|:--|:--|
| `/` | `docs/site/index.html` |
| `/assets/site.css` | `docs/site/assets/site.css` |
| `/scoreboard/` | `docs/index-calculus-scoreboard.html` |
| `/status/` | `docs/ecc2k130-status/index.html` |
| `/status/status.json`, `/status/history.json` | the snapshot, beside the page that reads it |
| `/status.json`, `/history.json` | the same bytes at the site root |
| `/404.html`, `/robots.txt`, `/sitemap.xml` | `docs/site/404.html`, generated |

Two constraints are worth knowing before moving anything:

- **`/history.json` must stay at the site root.** The publish job fetches the
  live `https://aburan28.github.io/crypto/history.json` to merge the previous
  snapshots forward, so moving that file discards every published snapshot
  before the move. `test_build.py` pins the URL in the workflow and the file in
  the output together.
- **The scoreboard is copied, never edited here.** AGENTS.md §7 makes
  `docs/index-calculus-scoreboard.html` canonical; the published page is a
  republish of it, and `test_build.py` asserts the two are byte-identical.

The dashboard's own JSON is published twice on purpose: the copies under
`/status/` let the page use relative fetches, so it renders the same way when
opened straight from the working tree.
