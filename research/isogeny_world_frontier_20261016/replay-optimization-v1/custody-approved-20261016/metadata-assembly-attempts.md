# Metadata assembly attempts

The native transfer completed with exit 0 before the first metadata checker ran. Its first version completed byte/hash comparisons, then incorrectly required all stderr files to be empty. Native fetch and verify commands write ordinary progress to stderr. That checker stopped at its stderr assertion. The original checker is preserved as `check-approved-custody-v1.py`.

The corrected checker accepts only the expected fetch lines and the exact `9 entries verified` line. It passed, retaining all raw stdout/stderr. Native pull and verify return plain text even with `--json`; their initially unsealed filenames were corrected to `.stdout` before sealing. No native transfer was retried or evidence bytes rewritten.

An administrative HA merge-hygiene help check initially used the standalone Python 3.13 interpreter, which lacked PyYAML and failed to import the ledger module. The existing project virtualenv provides Python 3.12.8 and PyYAML 6.0.3; its help check loaded successfully. No package installation, ledger mutation, cryptographic computation or scientific launch was performed to resolve this environment issue.

A later metadata-only progress read failed before Python executed: zsh could not create its here-document temporary file (no space left on device). The system Data volume reported 151 MiB available and 100% capacity; SSD990-2 retained 601 GiB and SSD990 retained 28 GiB. Subsequent metadata commands use bash and the task-owned temporary directory on SSD990-2. The active native root-suite log was retained on SSD990-2, and no failure line was present at diagnosis. No research artifact was deleted or rewritten as a workaround. This is a shell/filesystem failure, not an algebraic outcome.

Exact shell error, transcribed from the failed tool result (exit 1; no separate raw file was captured for that failed metadata read):

```text
zsh:1: can't create temp file for here document: no space left on device
```

A later capacity check reported 1.4 GiB free on the system Data volume. This task performed no cleanup there; its new scratch stays on SSD990-2.

A second metadata progress read also failed during bash here-document temporary-file creation, before Python executed, despite the in-command TMPDIR assignment. Exact emitted error: /bin/bash: line 1: cannot create temp file for here document: No space left on device. Diagnosis again showed system Data at 100% (138 MiB free), while SSD990 had 28 GiB and SSD990-2 had 597 GiB. Metadata readers now use stored scripts or Python -c, avoiding here-document temporaries entirely. The native log continued through the quasi-subfield tests with no failure output. This second shell failure is retained separately; the earlier claim that changing shell/scratch alone resolved all temporary-file creation was too broad.

The first JavaScript preparation for this correction failed to parse because a literal backtick was not escaped (SyntaxError: Unexpected token colon). It invoked zero nested tools and wrote no file. The corrected preparation removed that quoting conflict.
