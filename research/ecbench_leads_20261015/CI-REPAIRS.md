# CI findings while pursuing the leads

The prior head321a7a5e4 passed ecbench, Rust lint and all-target compilation,
but failed these applicable checks. These are repair obligations, not research
outcomes or permission to bypass checks.

| Check | Exact finding | Action |
|---|---|---|
|ic-leaderboard, job116123014736|Standalone native-publish Cargo.lock was ignored, so --locked fails on a fresh checkout|Commit the exact existing lock and fix ignore-rule precedence|
|curve-names, job116123016154|Two measured conductor-gap models were not in the registry|Register both exact models and subgroup/generator representations natively; preserve frozen old sessions|
|IC dashboard, job116123016516|Generated overview is stale|Regenerate the presentation from its pinned evidence and review desktop/mobile output|
|Curve cover certificates, job116123015482|Catalog display check rejects the native-augmented leaderboard|Make the native generated views and catalog consistency check agree without changing frozen measurements|
|leaf-m10 archive replay|Existing exact-file-set assertion fails|Retain the separately reproduced baseline blocker; do not rewrite its archive or disable its verifier|

PR1853 also has a merge conflict and an incomplete automated review because
its evidence volume exceeds Bugbot's review limit. Per AGENTS.md, Copilot was
asked to merge main into the head, without rebasing or changing frozen evidence,
in [comment6184186354](https://github.com/aburan28/crypto/pull/1853#issuecomment-6184186354).
The native follow-up work proceeds in new source/report paths. A conflict in
frozen evidence is reported, not resolved by choosing a side. No merge is
admitted until the exact final head satisfies the applicable CI/review gates.

After roughly50 minutes without a Copilot response, the author announced the
AGENTS.md fallback, documented in comment6184777039, to integrate main locally
without rebasing or force-pushing. Generated views preserve the incoming main
scoreboard and index its12 previously added ecbench sessions from their exact
fd8bcacf4 bytes, separately from this task's373 new attempts. No new verdict is
derived for those earlier sessions. Historical scientific source arrays remain
pinned to f9011e5; catalog additions preserve old entries and get an unchanged
native cover-certificate check. The existing catalog prose renderer only lays
out the native-generated document and diagnostic fragment.
