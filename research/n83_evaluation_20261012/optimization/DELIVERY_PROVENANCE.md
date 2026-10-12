# Measured source and delivery integration

The executed rho and controller binaries were frozen before timing. The
configured worker source is SHA-256
`fb14072f4fa5d08e18d8b79130b6d58a220adae0f76369cc9be949691c83eaa0`.
The executable is
`21827dcc645fd4591b39dd913619f047d267a98d541854bd7c83ee7143d8bbcb`.
The controller executable is
`b1b94b24fea42a7ac02e2ce21dfd913e7ba4a2b66a536f48d65bb3cfe8044a27`.
The IC executable remains the original study's
`983f28400b7852fc6628ba97d3ea83b2ab9e32d9944658c7eea831ee89a4abf2`.

Source/protocol commit4c9da3bd and final-controller commit4ae7caa2 precede all
timed children. Those commits remain reachable in this branch's history.
The earlier source snapshot is a8194d61 plus the previously documented collector
repair and new examples. No measured binary calls that collector. The compiled
worker's revision field remains null; source and executable hashes provide the
actual binding. Raw receipts do not infer a compiled revision from their filename.

After all timings and161 independent answer replays, the delivery branch
integrated main atc2963dcf via merge010e4f69. Main had already incorporated the
original study through a squash merge. Historical study files were aligned
before integration so that the follow-up adds only the new experiment and
its current-result pointer, rather than duplicating the earlier result PR.

Strong rho, WideGf2, Gf2_128, the original IC example and the new configured
worker sources did not change in that integration. Other upstream generic
binary-curve/library code did change. Delivery correctness checks therefore
run again on the integrated tree. These checks do not relabel its newly built
binaries as the earlier measured binaries or replace any frozen timing.

The first read-only merge-tree attempt could not fetch a promised Git object
through the default network sandbox. The same operation with explicit network
access succeeded. No automatic approval review rejected an action. All actual
integration was clean; no PR conflict or CI-monitoring flow was invoked.

Original14-page/48-page PDFs, raw five-pair records and original binaries remain
unchanged. The new addendum precedes them in a separate comprehensive PDF.

Final delivery checks passed: full release library3992 passed,0 failed,
112 ignored; configured example4 passed,0 failed. The initial integrated-tree
test failure and first focused rerun are retained. They exposed stale test
assumptions: Montgomery curves can now be counted, and extension descriptions
require an explicit modulus. Only the test block was corrected; no production
code changes. Native C++ enumeration gives95 affine points plus infinity for
the F101/A3/B1 Montgomery model, and rejects roots of x^2+2 over F37. The test
now proves order96, rejects wrong order101 and retains the extension's
NotEvaluated certificate boundary. Commands, statuses and all hashes are in
checks/delivery_test_only_fix_receipt.json.
