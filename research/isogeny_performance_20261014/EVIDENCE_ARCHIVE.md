# Frozen evidence and extraction

The full frozen study at commit
`b27fbfe8e0f2c53e2bf25221b3a59fb67d1a8482` is retained in
`FROZEN_LOCAL_RUNS.tar.gz`, with repository-relative paths. It includes every
raw JSONL, successful/failed test log, constructor comparison, repeated-run
receipt, source/binary archive, report source, PDF, graph and file inventory.

The archive was extracted to a fresh directory on the second SSD. The native
auditor checked **272 file byte counts and SHA256 digests**, all PASS. A native
recursive comparison also found no difference between the extracted evidence
and the original evidence directory.

From a fresh checkout, restore the original paths with:

```sh
tar -xzf research/isogeny_performance_20261014/FROZEN_LOCAL_RUNS.tar.gz -C .
```

The original inventories remain inside the archive. Their source context is
the frozen commit, and they exclude the inventory file being written. Historical
inventories apply to their historical source tree, not to subsequently added
publication notes. `FILES.json` is retained for direct inspection, and the
final constructor comparison and repeated-run summary remain directly readable
under `evidence/`. Links to other archived evidence work after extraction.

The initial automatic approval review rejected deleting the raw working copies.
The approved alternative changed only the Git index: **no raw file was deleted
or altered**. Every uncompressed copy remains in the owning local worktree;
earlier Git commits also retain the original blobs. The final PR diff packages
the same data as a verified archive to reduce the automated review footprint.

`RAW_EVIDENCE.tar.gz` is an additional byte-identical archive of the evidence
directory alone, extracted with `-C research/isogeny_performance_20261014`.
