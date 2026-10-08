"""Restore exact research artifacts from reviewable UTF-8 chunks.

Run from any directory. --verify checks all hashes without writing files.
Existing files are accepted only when they already contain the exact artifact.
"""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    manifest = json.loads(Path(__file__).with_name('evidence-parts.json').read_text())
    for name, entry in manifest.items():
        target = (root / name).resolve()
        target.relative_to(root)
        blocks = []
        for name in entry['parts']:
            part = (root / name).resolve()
            part.relative_to(root)
            blocks.append(part.read_bytes())
        value = b''.join(blocks)
        if len(value) != entry['bytes'] or hashlib.sha256(value).hexdigest() != entry['sha256']:
            raise ValueError('artifact checksum mismatch: ' + str(target))
        if target.exists():
            if target.read_bytes() != value:
                raise ValueError('existing artifact differs: ' + str(target))
        elif not args.verify:
            with target.open('xb') as stream:
                stream.write(value)
    print(('Verified' if args.verify else 'Restored') + f' {len(manifest)} exact research artifacts.')


if __name__ == '__main__':
    main()
