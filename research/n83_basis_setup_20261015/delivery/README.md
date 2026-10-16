# Evidence archive and reconstruction

The archive stores each distinct payload once under objects/<SHA256>.
RECONSTRUCTION_MANIFEST.json maps EVERY original included path, byte count,
content hash and executable mode to its object. Identical base files are
represented at every original process path. No scientific values are recomputed.

HASH_MANIFEST.json outside the archive repeats the complete path inventory.
EVIDENCE_ARCHIVE_RECEIPT.json records ZIP CRC and every-object SHA/size checks,
path counts and archive SHA/size. Bundle/receipt/hashmanifest themselves are
excluded from the payload scan to avoid a circular self-hash; the archive embeds
the reconstruction manifest. Operational completion-watch IDs/receipts and chat
transcripts are not payload. Measured raw records, fixtures, validation sidecars,
all audit receipts, reports/graphs, source/binary freezes and failed attempts are.

Reconstruct only into an empty chosen destination using the manifest's relative
path field. Repository files use their normal repository-relative paths;
external frozen artifacts use external/<absolute-path-without-leading-slash>.
Never overwrite a live campaign. Validate each object's SHA256 and size before
writing and preserve its mode. Original absolute paths are provenance, not an
instruction to write there. Sources retain version/commit and build pointers;
the archive does not assert hermetic system-library or toolchain closure.
