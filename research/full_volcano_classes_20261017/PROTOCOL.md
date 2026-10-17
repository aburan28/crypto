# Three additional complete ordinary classes

The user requested a few finished graphs uploaded to Google Drive. This run
adds three tractable examples; it does not replace the unfinished fixed-input
A/C/D/E classes from the 20261016 study or complete the requested ECDLP and
index-calculus campaign.

Before execution: construct ordinary CM examples for fundamental discriminants
-7, -11 and -19 with Frobenius conductors 8192, 6561 and 1728 respectively.
Search increasing traces near 2^29, 2^29 and 2^30, retaining the first prime
field and a model with certified order. Fix seeds to 1. Store exact integer
inputs and native plans before enumeration. Fields must lie between 2^51 and
2^61. The models start at the maximal order; after complete labeling, export
an additional model at the maximum conductor for each graph.

Use the existing native isogeny-class executable and automatic prime degrees.
Run to frontier exhaustion, then replay from a copy. Acceptance requires the
full predicted vertex count, all expected conductor populations, and every
stored outgoing kernel and explicit map identity to verify. Construction and
replay share arithmetic; independent implementation replay is NOT_RUN.
No new mathematical theorem or performance comparison is asserted.

Stop on correctness failure, 2 GiB generated data, or two hours per class;
retain incomplete outputs and errors. Record commands, source/binary hashes,
inputs, actual exit status, summaries, stratum counts and a visual/PDF report.
No graph-generation wall times are admitted as solver benchmarks. ECDLP
counts/times, min/max/mean, direct-IC relations and finite-field reductions
remain NOT_RUN for these new graphs.

Upload the finished SQLite graphs, inputs, reports, executable and source as a
new archive alongside the existing Drive delivery. Preserve existing files.
Verify uploaded size and metadata; distinguish that check from a remote
whole-file hash check. Commit compact evidence and durable links, not graph
databases, into the existing crypto PR.

Conductor command attempted before edits: unavailable (exit 127), and this
checkout has no .conductor/WORKFLOW.md. No reservation was obtained.
