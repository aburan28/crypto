import hashlib
import json
from pathlib import Path
from survey import survey

root = Path(__file__).resolve().parent
expected = (root/'results.json').read_text()
assert json.dumps(survey(), indent=2)+'\n' == expected
receipt = json.loads((root/'receipt.json').read_text())
for name, digest in receipt['sha256'].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest() == digest
print('PASS: six structural cases, exact identities, prime factors, and replay hashes.')
