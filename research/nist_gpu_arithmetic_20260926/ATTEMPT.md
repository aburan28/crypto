# RTX PRO 6000 execution status

The first hosted dispatch reached GitHub Actions as run 36221654987. Checkout
and the credential gate completed, but the repository does not currently expose
MODAL_TOKEN_ID / MODAL_TOKEN_SECRET to GitHub Actions. The Install Modal,
device-validation and tournament steps were therefore skipped by design.

This is an infrastructure non-result. It contains no GPU throughput
measurement and must not be cited as performance evidence.

The benchmark remains runnable in either of two ways:

1. provide Modal credentials to the workflow and update RUN_GPU once; or
2. run gpu/nist/modal_app.py::validate and ::tune on an authorized RTX PRO
   6000 host/cloud account.

The CPU oracle is a separate required gate and is not waived by this blocker.
