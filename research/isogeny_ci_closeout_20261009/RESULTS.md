# Isogeny CI closeout, 2026-10-09

Hosted checks started on PR 1606 at source `110d43c6b`. The dedicated Linux glibc,
Linux musl and Windows binary jobs passed their build, tests, installed-tool searches,
independent replay and forced-software SHA checks. Apple silicon remained queued at
that observation. These are results for that recorded head, not a later commit.

The standalone algorithm lint job reported manual divisibility predicates in the
search example and API integration test under Rust 1.98. They now use the equivalent
integer `is_multiple_of` method. The PCLMUL hardware flag is scoped to x86_64,
where it is actually used, so native ARM lint no longer reports an unread field.

Two canonical leaderboard checks reported stale JSON. The incremental native view
updater had emitted new roster-object keys alphabetically, while the canonical builder
requires insertion order. The native `canonical-roster` helper restores that order
and asserts that every parsed JSON value, including all measurements, is identical;
it also preserves measurement-board bytes. HTML and Markdown required no changes.
No experimental value, curve identity, graph edge or ledger promotion changes.

After that byte-order repair, the browser's recorded leaderboard checksum must track
the new canonical bytes. The native helper's `--browser-sources` mode refreshes the
four catalogue/cover/leaderboard digests and asserts that every other browser value
is identical. This is a downstream provenance repair, with no experiment or row
promotion. The hosted leaderboard check passed at repair head `8e15a1919`.

Run the native formatter with `canonical-roster docs/ic/leaderboard.json`, or validate
its formatting using `canonical-roster docs/ic/leaderboard.json --check`. Unknown
roster fields fail closed for review. The historical profile, source archives and
measurements remain frozen at their recorded revisions.

The root library and unrelated release jobs remain failing gates. Their failure is
separate from the passing dedicated binary jobs; the PR cannot merge or publish while
applicable checks remain unsuccessful. New CI must validate the repair commit before
any earlier-head result is treated as acceptance for it.
