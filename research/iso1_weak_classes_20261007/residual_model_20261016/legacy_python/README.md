# Historical Python evaluation (superseded)

These scripts produced the first V1/V2 evaluations and the forward predictions
committed in `edc0c3ddf` and `788e11c8c`. AGENTS.md ("Implementation language:
no Python") forbids Python for research result generation, so the execution
path is now the native `iso1_residual_model` binary
(`../../census_runtime/src/residual_model.rs`). It reproduces every number these
scripts printed for p = 7–19 exactness, p = 37–61 holdouts and all six forward
predictions. The scripts are kept unchanged as provenance for those commits;
do not extend them.
