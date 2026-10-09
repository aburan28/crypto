# perf-history

Performance snapshots of `aburan28/crypto`, one per measured commit, appended
by `.github/workflows/perf-history.yml` on every push to `main` and rendered
by `scripts/perf/perfhistory.py` (both on `main`).  This branch holds data
only; never merge it into `main`.

- [`docs/perf/HISTORY.md`](docs/perf/HISTORY.md) — the tables (areas,
  snapshots, movers, every kernel) and the chart
- [`docs/perf/history.svg`](docs/perf/history.svg) — chained instruction-count
  index per area over the measured commits
- [`docs/perf/history.json`](docs/perf/history.json) — the same, machine-readable
- `docs/perf/history/*.json` — the raw snapshots

```bash
git fetch origin perf-history
git show origin/perf-history:docs/perf/HISTORY.md
```
