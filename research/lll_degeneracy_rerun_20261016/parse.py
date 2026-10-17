"""Parse the probe logs into results.json (one row per probe line).
Records only what the logs say; no computation beyond counting."""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
LINE = re.compile(r"^(?P<curve>\S+?)(?:\[(?P<path>f64|HP)\])?\s+k_bits=\s*(?P<k>\d+) m=\s*(?P<m>\d+)"
                  r".*?d_seed=0x(?P<seed>[0-9A-F]+).*?→\s+(?P<outcome>[✓✗][^(]*?)\s+\((?P<ms>\d+) ms\)")
rows = []
for log in ("probes.log", "p521_hp.log"):
    test = None
    for line in (HERE / "logs" / log).read_text(encoding="utf-8").splitlines():
        if line.startswith("test ") and " ... " in line:
            test = line.split()[1]
        m = LINE.match(line)
        if m:
            path = m["path"] or ("HP" if "reduction=LLL-HP" in line else "f64")
            rows.append({"log": log, "test": test, "curve": m["curve"], "gs": path,
                         "k_bits": int(m["k"]), "m": int(m["m"]), "d_seed": "0x" + m["seed"],
                         "recovered": m["outcome"].startswith("✓"), "outcome": m["outcome"].strip(),
                         "ms": int(m["ms"])})
summary = {}
for r in rows:
    key = f'{r["curve"]} ({r["gs"]})'
    s = summary.setdefault(key, {"recovered": 0, "probes": 0})
    s["probes"] += 1
    s["recovered"] += r["recovered"]
(HERE / "results.json").write_text(json.dumps({"rows": rows, "summary": summary}, indent=1, ensure_ascii=False) + "\n")
print(json.dumps(summary, indent=1))
