# Frozen pre-repository evidence

These files are unchanged copies of evidence produced before work entered this
repository. They are not new ecbench sessions or newly qualified timings.

The complete historical case B degree-3 graph has 129,140,163 vertices,
18 conductor strata and 258,280,326 directed kernel maps. Root horizontal
kernels and strict dual maps are part of its lossless encoding; there are
additional prime-degree maps if the finite degree set is enlarged.

- [Download folder](https://drive.google.com/drive/folders/1J20JzTGAkOt594wu4eDSP0vtnu2Yp_kX)
- [96 MiB graph parts](https://drive.google.com/drive/folders/1pTN-Ov-4Q2N8nqZaTsKQC_WuRgS1cqyu)
- `drive_file_index.json`: canonical 79-part URLs, bytes and part hashes.
- `parts_manifest.json`: exact concatenation order and compressed-file hash.
- `dataset_manifest.json`: raw/compressed hashes and record schema.
- `checkpoint.txt`: completed native checkpoint, including both root horizontal kernels.
- `native_validation.json`, `validation.json`: historical audit coverage.

Preserve the extra remote retry copy of part 0038; it was renamed
`RETRY_COPY_NOT_FOR_REASSEMBLY — vertices.bin.zst.part0038`. It is not in the
canonical 79-part index and must not be included in reassembly. No remote
objects were deleted for this repository integration.

The frozen `reassemble.py` is archival download glue. It is not the new native
graph engine or a new Python cryptographic execution path. Reconstruction
requires downloading the canonical parts; the enormous graph bytes are not in Git.
The source software archive is in the download folder. New large artifacts
must use the repository store workflow and retain content-addressed provenance.
