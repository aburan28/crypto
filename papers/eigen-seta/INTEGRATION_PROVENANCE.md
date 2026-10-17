# Historical manuscript provenance

This paper, its figures, source and PDF are preserved from PR #1680 at commit `74bec88bdc6fb43724493d452d72df43ca3ea735`. It is a historical research draft, anonymized, not peer reviewed and not submitted. Its native files were retained unchanged.

The maintained implementation in [`eigen_seta/`](../../eigen_seta/) includes later KEM–DEM PKE, public-key validation and scalar-recovery changes from previously merged PR #1672. Statements in the preserved manuscript that its specification matches the active implementation exactly, that every secret byte buffer is wiped, or that its instruction-count measurements certify the current build must be read as statements about the historical source. The current [`SIDE_CHANNELS.md`](../../eigen_seta/SIDE_CHANNELS.md) documents the active limitations, including unwiped temporary buffers.

The second verification cohort binds source `d504af28e0a262cd97cdca0a58c8168e90ce6f60`, and its public-filter example binds `72051ac67`. The manuscript's cryptographic security categories remain conjectures; no independent cryptanalysis or new security certification follows from this integration. Its mathematical cyclic-kernel argument is preserved without promoting it into maintained claims during conflict resolution.

See [integration reconciliation](../../research/eigen_seta_rust_20261011/INTEGRATION_RECONCILIATION.md) and the [byte-preserving archive manifest](../../research/eigen_seta_rust_20261011/archives/pr-1680-74bec88bd/MANIFEST.json). Updating the historical manuscript to specify and measure the active build is a separate substantive revision.

## Revision 2 (2026-10-17)

The manuscript was revised in place: it now specifies against and measures
the maintained implementation in the sibling repository `eigenx` (commit
`94fe2d9`), adds the assessment of §9 and withdraws the unqualified NIST
category conjecture. The revision-1 measurements (§8.1–8.3) are kept and
labelled as measurements of the historical `crypto:eigen_seta/` build. The
revision-1 PDF remains in git history at the commit named above.

