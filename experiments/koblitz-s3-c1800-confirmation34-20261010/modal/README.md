# Public input bundle for the N53 C1800 Modal confirmation

[`PUBLIC_BUNDLE.json`](PUBLIC_BUNDLE.json) packages the 34 frozen T213–T246
public points, exact candidate and binary hashes, factor-base identities,
algorithm parameters, and all 136 ABBA/BAAB run-order entries. It is
14,824 bytes with SHA-256
`246f308821952be73b691467bb7227e27090eaeb645e1af8de36a9d63bd9f5cd`.
The bundle is built from and bound to the existing
[`PANEL_FREEZE.json`](../PANEL_FREEZE.json) and
[`EXECUTION_FREEZE.json`](../EXECUTION_FREEZE.json). It contains no fixture
scalar, scalar derivation record, repository source, or credential. The
expected scalars remain local for independent result replay.

[`prepare_public_bundle.py`](prepare_public_bundle.py) checks every point-file
and workload-manifest hash, distinctness, exact source and manifest
identities, and the complete run schedule before creating the bundle. Its
[`receipt`](PUBLIC_BUNDLE_RECEIPT.json) binds the generator source and output
hashes. Run `python3 prepare_public_bundle.py --verify` from this directory to
recompute and verify the committed bytes without writing or contacting Modal.

This is an input artifact for a later confirmation runner; it does not create
a Sandbox or transfer any data. The bounded T205
[Modal pilot](../../koblitz-s3-c1800-combined-20261010/modal/README.md)
should first establish Linux execution, guest telemetry, and observed run
duration. A separate, reviewed resource and price freeze and transfer manifest
will then determine how to partition 136 runs into bounded Sandboxes. Modal's
guest measurements remain exploratory unless a host-wide isolation receipt
passes the repository's CPU gate.
