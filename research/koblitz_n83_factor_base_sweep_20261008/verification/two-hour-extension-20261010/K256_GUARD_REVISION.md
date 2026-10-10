# K256 guard correction before launch

The frozen K256 plan and input/caps remain unchanged. The original launch
shell copied a cleanup check that could incorrectly report success when
`docker ps` itself failed. That was exposed by the retained m6 Docker daemon
EOF. No K256 worker had been launched. Preserve both shell versions and
their SHA-256 hashes. This revision accepts cleanup only when `docker ps`
returns zero and its exact-name filtered output is empty; a failed `docker
wait` is classified as producer failure. No experiment outcome is inferred
from the m6 container loss. The m6 recovery audit supersedes the invalid
cleanup flag in its original outer receipt.
