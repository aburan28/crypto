# Modal source export was denied by automatic approval review

On 2026-10-13 the command `modal run
research/ecc_arithmetic_20261013/modal_panel.py` was retried with elevated
network access after the sandboxed attempt could not connect. The tool did
not start the command. Automatic review returned:

> The Modal launch uploads the repository's study files and live
> `src/binary_ecc/f2m.rs` into a third-party Modal worker via
> `add_local_dir`/`add_local_file`; although the user authorized using Modal
> for benchmarking, they did not specifically authorize exporting this
> sensitive source payload to that destination, so the sensitive-egress rule
> requires denial.

The review also instructed not to bypass the rejection through an indirect
execution. No remote run receipt exists. The original script used
`add_local_dir` for the study directory. It has since been narrowed to the
study Cargo manifest and lock, `src/main.rs`, `src/bin/compare.rs`,
`baseline_src/f2m.rs`, and the frozen `candidate_src/f2m.rs`, as now declared
in [`modal_panel.py`](modal_panel.py). Modal also sends the launcher function
source. The exact files are listed in [`MODAL_BUNDLE.sha256`](MODAL_BUNDLE.sha256).
No source was uploaded in the denied attempt. Specific authorization for this
frozen source export would be needed
before launching that script on Modal.
