This is the first successful frozen static preflight, run from initial draft head
`c901d44` after its hash-only GitHub CI passed. It audited PATH binaries only,
so WDSat appeared unavailable even though the separately pinned #764 fixture
binary existed outside PATH on this Mac. Its source/corpus facts and blocker
verdict remain valid; its binary inventory is incomplete for the host. The
original result, stdout/stderr and receipt are retained byte-for-byte. A
pre-outcome source amendment explicitly inventories the #764 fixture path and
re-freezes before the final static preflight. This is not a solver run.
