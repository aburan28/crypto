# N83 four-point cold factor-base confirmation

K1152 with the dynamic one-root five-byte exact index has the lowest
descriptive cold process mean across four frozen reserve points: 164.431 s.
Its 0.480 s mean advantage over K1152/root22 is positive on every point,
and the paired two-sided 95% interval is 0.108–0.852 s. The frozen winner
rule additionally required every per-point advantage to exceed half that
point's K1024 A/A spread. On point 1 the advantage was 0.161 s against a
0.521 s threshold. The selected pipeline and candidate ID therefore remain
`null` with status `BOUNDED_UNKNOWN`.

| Pipeline | Actual base points | Folded columns | Mean cold process, s | Difference from K1152/root11, s |
| --- | ---: | ---: | ---: | ---: |
| K1152/root11 dynamic | 191,232 | 1,152 | 164.431 | 0 |
| K1152/root22 dynamic | 191,232 | 1,152 | 164.911 | +0.480 |
| K1024/root22 dynamic | 169,984 | 1,024 | 166.995 | +2.564 |
| K1536/root22 dynamic | 254,976 | 1,536 | 186.012 | +21.581 |

The [frozen protocol](PROTOCOL.md) used four distinct reserve public points,
five fresh IC processes per point in rotated order, and four same-point rho
controls. All 20 IC arms, four rho controls, four K1024 A/A controls,
cross-point base identities, execution order, and resource gates passed the
independent checker. The four A/A relative differences were 0.622%, 0.382%,
0.022%, and 0.215%. The matched rho results and the target-dependent online
phase timings are retained in the [analysis](analysis_38890911296.json) and
[20 IC1 run rows](ic1_rows_38890911296.json). This cold panel uses source
SHA-256 `d094345474605075566d646ddecef8f694cdd644b0f9190a6b1ef1551fd06888`;
the separate current-source online range screen uses a different merged source
and is analyzed as its own workload.

The [GitHub run](https://github.com/aburan28/crypto/actions/runs/38890911296)
completed successfully and uploaded artifact `12045705683`. Its exact ZIP
is archived locally at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-cold-reserve-38890911296.zip`
with SHA-256
`bcf6c55708a7e8026d760205810c47de780d97b2904d2b4a87310457b5e4b9c3`.
The cloud analysis SHA-256 is
`28253a17d18b6951843eca6451782e590b6d7d07c9b7a54ec96e9850bffea45d`.
Local independent replay against the extracted raw artifact exited zero and
recreated that analysis byte for byte. The local replay log SHA-256 is
`25b9a1f631ffed166278f80cdbdf71e5ea5c478a0043d4511a16716d8aa940e7`.

The next selection step is the current-source expanded online range screen
on the previously exposed point, followed by the sealed five-fresh-point
Tier 1 same-point IC/rho panel. The cold result remains a setup-cost and
factor-base geometry comparison under its recorded cloud allocation.
