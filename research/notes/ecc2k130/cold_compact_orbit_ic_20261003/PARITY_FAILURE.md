# Constructor parity attempt 1

Before any fresh-Q measurement, the initial `compact-orbit-scan` constructor
was run against all 3,108 archived `SOURCE42.jsonl` points with the pinned
lockfile. It used `BinaryInstance::points_with_x`, whose linear
Artin–Schreier solver may return either valid root. The release test failed at
point index **666**, the first point of column 9:

| source | x | y |
|---|---:|---:|
| initial framework constructor | 107203939703 | 63515710123 |
| frozen archived producer | 107203939703 | 95496871900 |

The two y values are negatives on this binary curve. That is a support-order
and base-hash failure, even though the abscissa agrees. The repair uses the
archived producer's half-trace lift order, counts its native field work
without borrowing the framework's different Artin–Schreier calibration,
and reruns full point, label, representative and BLAKE3 hash parity. No
public-Q outcome was inspected or measured before this repair.
