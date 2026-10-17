# Native independent rank/log audit

`verify_rank.cpp` uses separately implemented polynomial GF(2^83) field and
elliptic-curve arithmetic, independently replays every modular row/rank gain,
backsolves and checks every representative log and four-point answer witness.
It never links the measured Rust implementation. Its command is:

```
g++ -std=c++17 -O2 -Wall -Wextra verify_rank.cpp -o verify_rank
./verify_rank BASE.jsonl TRACE.jsonl ANSWERS.jsonl EXPECTED_RANK_SEED NEW_RECEIPT.json
```

Meter this checker separately after releasing the benchmark reservation. A
successful receipt binds exact base/trace/answer/checker bytes; campaign-level
source/binary/matrix binding is an additional acceptance requirement. A bare
checker PASS is not a complete experiment or a performance win.

`test_rank.cpp` constructs a two-column public synthetic algebra-only fixture
with fresh known logs and checks one positive and fourteen tampered/truncated
cases. This fixture is never a performance or selection input. The first test
invocation failed while hashing a nonexistent synthetic argv[0]; it was fixed
to hash the real test executable. Both logs are retained.

The SHA and previous-fixture helper headers are unchanged copies
from the completed matched study. Assertions must remain enabled; never build
with `-DNDEBUG`. The new arithmetic-header copy additionally rejects duplicate JSON keys;
its fourteen negative tests include that case. The parser is limited to these
native artifact schemas; this is not a general untrusted JSON service.
