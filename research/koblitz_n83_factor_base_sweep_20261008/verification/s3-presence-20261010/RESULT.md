# N83 factor-base S3 presence check

A read-only `ListObjectsV2` request against `s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/` found all **57 expected compressed factor-base object keys**. The listing contained 72 keys overall, including 15 non-base study objects. The 57 base objects total 109,656,159 compressed bytes in this listing. [The raw S3 listing](listing.json) and [machine receipt](receipt.json) bind the listing response hash, four committed source receipts, each expected key, its listed size, modification time and ETag. No expected key was missing.

This check proves key presence at the observation time. It does not download or rehash object contents; the original upload/download hash checks remain in the source receipts. It does not add relation, rank, individual-log or total-runtime evidence and consumes none of the exhausted 7,200-second local experiment budget.
