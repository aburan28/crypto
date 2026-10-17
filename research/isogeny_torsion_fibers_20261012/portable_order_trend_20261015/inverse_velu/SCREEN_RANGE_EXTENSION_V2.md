# Fourth order-bin screen range extension, V2

Frozen after the completed V1 screen and before changing its executable
or running any extension. V1 inspected all eligible primes through
`p=10007` and independently checked 513 rows; its fourth frozen order
bin `[19200,38400)` remained empty. The first three cells are fixed as
`(p,r)=(2473,2437),(5521,5569),(9619,9679)` and will not be
reselected.

Continue the same exact curve family, base-field order calculation,
different-`j` and order-three-kernel tests, factorization rule, and 64
sampled map checks on every eligible prime in the nonoverlapping range
`10009<=p<=30011`, `p=1 mod 3`. Select the first eligible prime `p`
with a multiplicity-one, prime-to-characteristic subgroup factor
`19200<=r<38400`; ties at one `p` choose the smallest such factor.
Retain every row and rejection. Stop at the range end or a 3600-second
wall or 512 MiB RSS ceiling; a cap is censored. Independently replay
every row and accepted map witness in native Rust before the selected
cell may enter a solve panel.

This is a structural input search, not a relation solve or speed
measurement. A fourth cell still requires exact generator and fiber
audits, public target custody, the N1 and N64 five-arm solves, complete
pricing, and independent replay under `FOUR_GATE_PROTOCOL_V1.md`.
Preserve V1 source, binary, raw rows, checker and report without edits.
