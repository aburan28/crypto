# First K2048 launch: producer failure before object creation

The guarded construction returned `PRODUCER_FAILURE_worker_exit` with worker exit 1 and no manifest. The worker stderr was `No such file or directory`. The supervisor's local `python:3.11-slim` tag resolved to image ID `sha256:0dd364ba7e10242f07755449e3a3d0e35f9efd987952737b90def6709ab0c5ce`; a no-network, read-only check showed that this image has no `git` executable. The exporter calls Git before constructing a v2 object to attest its clean source. A separate check found `/usr/bin/git` in the previously successful, still-local image `sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31`.

The failed run's config, outer receipt and empty/stdout/error logs are retained. No object, replay or S3 upload occurred. The first round charged 101 conservative seconds, bringing the cumulative charge to 6,440.092994294 of 7,200. It is a producer failure, not a mathematical or solver outcome.
