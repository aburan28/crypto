# Frozen source snapshots for the paired n37/n41 clean panel

The committed panel (`../../paired_fullrank_clean_evidence_20260925/`) pins the
SHA-256 of the 13 source files its six runs built, at the run's checkout,
`a045a0a80d3a8af00fedc85c0342f4e7d026aacd`. Seven of those files have since
been edited by later work. `../verify_committed.py` accepts a pin when the file
in the tree still matches it, or when a snapshot here holds exactly the pinned
bytes, and fails closed otherwise.

`SNAPSHOTS.json` uses the schema of the compact-orbit store
(`../../compact_frozen_source_replay_20260929/`), and its loader validates each
entry: the gzip's hash, the decompressed length and the content hash against
the pin. Each snapshot was taken with `git show a045a0a8:<path>` and checked
against the panel's pin before it was written.

This changes no measurement. It keeps the panel's source verifiable after
the files it pins move on.
