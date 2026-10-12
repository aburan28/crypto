# Bounded IC autolab: three-round result

The preregistered three-attempt improvement budget is complete. No selected
challenger passed the full promotion rule, so the qualified `pairinv` incumbent
remains selected. This is an audited negative result on the registered
synthetic toy panel, not evidence that the incumbent is globally optimal.

| Attempt | Selected challenger | Verified native/profile pairs | Confirmation challenger / incumbent cold Ir | Confirmation challenger / incumbent cold native time | Decision |
| --- | --- | ---: | ---: | ---: | --- |
| [One](improvement/round1/EVIDENCE.md) | `stop6` | 3,243 / 3,243 | 0.956822 | 0.964563 | Retained incumbent |
| [Two](improvement/round2/EVIDENCE.md) | `stop5_word` | 3,480 / 3,480 | 0.922725 | 0.987022 | Retained incumbent |
| [Three](improvement/round3/EVIDENCE.md) | `stop7_word` | 3,480 / 3,480 | 0.990182 | 0.985697 | Retained incumbent |

Each row is a separate registered attempt with its own fresh confirmation
points, frozen source and archived decision. Altogether 10,203/10,203
native/profile pairs were verified across the three rounds. Confirmation
required at least 20% lower complete cold instructions **and** cold native
time (each ratio at most 0.8), no single-target online regression, nominal
familywise upper bounds below one in confirmation and replay, no per-cell
regression above 10%, and independent target certificates. No selected
challenger reached either 0.8 cold point-estimate threshold. Round two's
largest cold-instruction improvement was accompanied by online regression;
round three's modest online point-estimate improvement did not make its cold
costs competitive. Online ratios from attempt one and attempts two/three use
different qualified IC reference roles, so the table does not merge them.

The [measurement contract](../MEASUREMENT.md) and
[version-two protocol](improvement-v2/PROTOCOL.md) bound exact curve,
candidate, workload and run identities; actual subgroup-usable factor-base
sizes before orbit folding; exclusive phase costs; failures and unknowns;
one supplied point per online interval; and separate cold/online IC and rho
references. The [repo skill](../../../.agents/skills/ic-autolab/SKILL.md)
now points to all three sealed rounds and forbids redispatch. The
[scoreboard](../../../docs/index-calculus-scoreboard.html) lists complete
stage tables, `S`, applicable collector-floor ratios, paired rho context and
the retained incumbent. Each round's durable archive and replay commands are
linked from its evidence note.

The competitive bounded panels exercised optimized three-summand pair-table
point decomposition, Frobenius-orbit factor bases, relation collection,
scalar-field row kernels and direct target descent as complete pipelines.
The broader generic autolab has F4/F5/SAT and sparse relation-LA adapters and
correctness controls, but those backends were not competitively qualified by
these three rounds. Base-only memory and the physical CPU model remain
unknown. All measured curves are small synthetic instances on a virtualized
Linux runner. These results establish neither an ECC2K-130 crossover nor a
globally fastest index-calculus implementation.

The next useful campaign is a **new** preregistered, source-bound complete-DLP
comparison that first qualifies generic F4/F5/SAT decomposition and larger
sparse relation-LA paths against the optimized incumbent and strong rho on
identical public targets. It must exclude every point exposed by all three
archives, declare a new resource envelope and operation boundary, and retain
failed attempts. The existing confirmation and replay data cannot be used to
retune or replay this exhausted three-attempt budget.
