# WeakCurves website

This is the static source for the first evidence-graded WeakCurves atlas.
`dist/` contains the complete site and a deterministic JSON snapshot. The
snapshot and criteria registry originate in
[`research/weak_curve_leaderboard_20261012`](../../research/weak_curve_leaderboard_20261012/).
The two copies of each JSON file must remain byte-identical.

The private hosted Site is registered by `.openai/hosting.json` and currently
deploys at [weakcurves.aburan28.chatgpt.site](https://weakcurves.aburan28.chatgpt.site).
The requested `weakcurves.com` hostname has been added but awaits domain DNS
validation. The hosted Site's access remains owner-private.

To rebuild the data after a source change, run the Rust exporter from its
research directory, copy `snapshot.json` and `criteria.json` into `dist/`,
then check `node --check dist/app.js` and both file comparisons. The Sites
source workflow pushes the exact site source and packages its `dist/` tree;
the research [validation record](../../research/weak_curve_leaderboard_20261012/VALIDATION.md)
names the current site commit and deployment.
