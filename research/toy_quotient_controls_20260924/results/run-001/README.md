# Superseded audit instrumentation run

Retained for provenance, not accepted as the canonical F5B audit.
All 1,920 decomposition results passed the point oracle. However, the original
F5B budget hook could raise an exception during generator cleanup. Its partial
F5B work receipts are therefore superseded by run-002.

`source/independent.py` is the exact audit source used here; its SHA256 matches
the corresponding entry in raw.json.gz. The other measured source files remain
identical to the paths recorded in that raw artifact. `source/contract.json`
preserves protocol version 2. No success threshold or input case was changed
when fixing the budget checkpoint. The corrected run is stored separately.
