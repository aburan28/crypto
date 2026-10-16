#set page(paper: "a4", margin: 19mm)
#set text(size: 11pt)
#align(center)[#text(size: 18pt, weight: "bold")[Official isogeny publication follow-up]]
#align(center)[16 October 2026 UTC · operational observations]

The authorized challenge merge is 05e414f3fd90 (PR1940). Its four native isogeny builds succeeded. ReleaseMe run 38857826201 has no published v0.1.2356 release at this follow-up; successful build artifacts and a tag do not establish publication.

#image("flow.svg", width: 100%)

= Exact failure and correction

The first two publication attempts failed on Git tag push with GitHub HTTP 500. Creating the exact tag through the refs API at the built 05e commit allowed attempt 3's tag step to succeed. Its publication then failed at 08:00:49Z: the pinned ReleaseMe action rejected override 0.1.2356 because prefix v was configured. The terminal failed-job log and run identity are retained beside this report.

The action's pinned #link("https://github.com/dev-build-deploy/release-me/blob/cede87fd8a4829ef571dff6bde807ebbfcaeb3e9/src/versioning.ts")[primary source] parses the full version with its prefix. Current main already contains #link("https://github.com/aburan28/crypto/pull/1981")[PR1981]'s correction: the version input uses the prefixed meta tag. Its exact patch and source ancestry are preserved. No additional release was dispatched here.

= Current boundary and next check

The newer #link("https://github.com/aburan28/crypto/actions/runs/38873645089")[main run 2372], source 3714752ee6a9, started 08:25:58Z and was observed building during this follow-up. It contains the 05e challenge merge and the prefix correction. Publication remains unverified. The previously observed latest published release is the older 11 October build v0.1.0-main.1371.

The corrected e82 experimental binary, SHA256 prefix b3cbd6ac5d45, and complete control archive for 94 tests have verified research custody. That is separate from an official ReleaseMe publication. The next scheduled publication check is around 09:15 UTC. Verify terminal jobs, tag/source binding, all four isogeny archives, SHA256SUMS and a downloaded native control before reporting official availability.

This is an operational release diagram, with no mathematical finding or timing comparison. Challenge scores and canonical IC/rho graphs are unchanged.
