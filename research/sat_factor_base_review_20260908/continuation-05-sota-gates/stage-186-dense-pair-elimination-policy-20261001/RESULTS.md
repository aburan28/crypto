# Stage 186: elimination policy after dense-pair selection

Both arms used the selected dense pair updates, authenticated the same target, and returned exhaustive UNSAT with their exact elimination counters.

| candidate / control | wall | core | RSS | performed XORs |
|---|---:|---:|---:|---:|
| full M4RI / current BlockTables | 1.847837 | 1.018215 | 0.856223 | 0.671154 |

Full M4RI reduces actual XORs and RSS but increases both wall and CPU on the selected dense-pair stack. The screen continuation condition fails, so no confirmation panel is run.

The two screen processes charge 197.916451 wall seconds, 657.884130 core-seconds, and 3972448256 bytes maximum RSS.

Current five-column BlockTables plus dense exact pair selection remains the Phase B target-specific and repository-default configuration. This is not a full attack or SOTA result.
