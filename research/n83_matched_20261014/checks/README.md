# Pre-timing validation

The release library suite passes4089 tests, with0 failures and117 ignored,
in335.44s. The new strict rho worker passes4 tests; the native controller
passes9, including actual process-group timeout and memory-envelope failure
checks. The native analysis passes14 tests including missing verification,
malformed/duplicate matrices, retained failures, conventional quantiles,
non-finite modeled crossovers and target-scaled jump costs.

Controller and worker raw logs are under `native_controller/`; native analysis
logs are under `../analysis/checks/`. Earlier missing-toolchain, compile and
floating-point test-assertion failures remain visible. No earlier failed run
is substituted for a pass. Compiler/source/executable pins and the excluded
old/new singleton outcome preflight are recorded separately before timing.

The [PR-CI fix receipts](pr_ci_fixes/README.md) reproduce19 historical archive
label failures, preserve all archived bytes, and verify the repaired policy
on active and frozen paths. Seven metadata tests pass. The Lean-source gate
rejects all six forbidden construct categories and fails closed without rg;
the workflow explicitly installs rg. This is local validation, not a claim
that newly triggered hosted checks already passed.

Conductor is unavailable (exit127); exact pre-edit scope attempts are under
`conductor.*`. Raw log bytes are preserved using a local Git attribute for
final blank lines. Git whitespace validation still applies to source code.
