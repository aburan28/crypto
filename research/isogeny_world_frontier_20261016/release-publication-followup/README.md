Official publication follow-up
==============================

The originally authorized challenge merge is main05e414f3fd908bde38f0e233b4d850270c025b07,
PR1940. Its ReleaseMe workflow38857826201 built all four isogeny targets. Attempts1
and2 failed during Git tag push with GitHub HTTP 500; their prior receipts remain.
The exact v0.1.2356 lightweight ref was then created through the refs API at that
same built commit, and only the failed release job was retried.

Attempt3 is terminal failure: tag pin/push succeeded; Publish with ReleaseMe
failed at2026-10-16T08:00:49Z with "Incorrect SemVer, missing prefix". No official
v0.1.2356 release was found at the scheduled08:42 follow-up. The raw failed log
and run identity are retained here. Successful builds and an existing Git tag
are not a published release.

The pinned ReleaseMe action cede87fd8a4829ef571dff6bde807ebbfcaeb3e9 reads prefix v
and calls SemVer.fromString(version,prefix). The override 0.1.2356 lacks that
prefix. Primary source:
https://github.com/dev-build-deploy/release-me/blob/cede87fd8a4829ef571dff6bde807ebbfcaeb3e9/src/versioning.ts
The version input must be the prefixed tag, such as v0.1.2356.

Current main already contains the correction, PR1981, commit
1bdc11680a86e89d5aaab1a6641f9417f0e7b8d9: version uses needs.meta.outputs.tag.
main-prefix-fix.json preserves that exact landed patch. The newer main Release
run38873645089 (run_number2372, source 3714752ee6a9686e0fce79cb308423d9fd09cde9)
was observed in_progress during this follow-up, after starting08:25:58Z. Its
publication remains unverified. The saved compare confirms that source contains
the authorized05e challenge merge. No duplicate release was dispatched.

The previously published /releases/latest observed at this follow-up is
v0.1.0-main.1371, published2026-10-11T15:52:54Z. It is not evidence for this new
challenge/compact publication. The corrected experimental e82 executable and
its independent control custody are documented separately under../replay-profiling/;
those are research artifacts, not an official ReleaseMe publication.

Next scheduled publication check is around 09:15 UTC for the newer run. Before
calling the tool officially available, verify terminal current-source jobs,
release/tag commit, all four isogeny archives, SHA256SUMS and a downloaded native
control. Neither cancelled jobs nor missing assets qualify as success.
