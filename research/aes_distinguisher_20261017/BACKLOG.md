# Native AES follow-up backlog

These nine cards translate the advisory hypotheses into the repository's backlog format. All remain **blocked proposals**: their extension driver, independent checker and exact frozen execution protocol are missing. No card is ready, executed, leased or published to cairn. Funding units express a proposed relative ordering only.

| Priority | Card | Next discriminating test |
| ---: | --- | --- |
| 1 | [instrument-null-panel](../../backlog/XP-20261017-aes-instrument-null-panel/README.md) | 32 keys, fixed disjoint structures, 2^17 examples/key, exact cipher fixtures and no-op/ShiftRows/final-mixing mutants |
| 2 | [toy-hull-versus-trail](../../backlog/XP-20261017-aes-toy-hull-versus-trail/README.md) | complete 16-bit toy domain, <=16 keys, <=8 frozen input differences, rounds 1..5 |
| 3 | [joint-related-quartets](../../backlog/XP-20261017-aes-joint-related-quartets/README.md) | retrieve complete algorithms, then full 16-bit toy domain on <=8 keys and <=2^18 full-byte structures/key |
| 4 | [cross-ratio-superbox](../../backlog/XP-20261017-aes-cross-ratio-superbox/README.md) | exhaustive GF(16) quartets and <=2^20 selected GF(256) quartets, then <=2^18 multi-byte structures/key on <=16 keys |
| 5 | [boomerang-dl-dependencies](../../backlog/XP-20261017-aes-boomerang-dl-dependencies/README.md) | 2^24 S-box BCT checks; full toy switch domain for <=16 frozen difference pairs and <=16 keys |
| 6 | [integral-certificates](../../backlog/XP-20261017-aes-integral-certificates/README.md) | exact 4-bit S-box transitions, <=256 frozen sets with <=2^16 points, 32 validation keys |
| 7 | [schedule-isolating-panel](../../backlog/XP-20261017-aes-schedule-isolating-panel/README.md) | 64 keys, <=2^16 structures/key, <=4 frozen statistics, several reduced round counts |
| 8 | [sat-witness-and-unsat](../../backlog/XP-20261017-aes-sat-witness-and-unsat/README.md) | all 8-bit S-box relation pairs; <=32 toy cells with <=5 rounds, 300-second watchdog/cell |
| 9 | [learned-native-proposal-ranker](../../backlog/XP-20261017-aes-learned-native-proposal-ranker/README.md) | <=2^20 training examples, 2^18 validation examples, <=2^20 final examples over >=64 final keys, <=3 frozen models |

Source mechanisms, controls, falsification and novelty limits remain in [hypotheses.json](hypotheses.json); source identifiers resolve in [sources.json](sources.json). Native `backlog validate` checks card structure, not scientific truth or readiness. Existing official AES campaign queues and permissions remain unchanged.
