"""Verify numbered graph parts and reconstruct vertices.bin.zst."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--directory',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--verify-only',action='store_true')
    args=ap.parse_args();manifest=json.loads((args.directory/'parts_manifest.json').read_text())
    output=args.output or args.directory/manifest['original_file']
    if output.exists() and not args.verify_only:
        raise FileExistsError(f'{output} exists; use --verify-only or choose another --output')
    temporary=output.with_name(output.name+'.assembling')
    if temporary.exists() and not args.verify_only:
        raise FileExistsError(f'{temporary} exists; inspect or remove it before retrying')
    whole=hashlib.sha256();total=0
    target=None if args.verify_only else temporary.open('xb')
    try:
        for part in manifest['parts']:
            p=args.directory/part['file'];digest=hashlib.sha256();size=0
            with p.open('rb') as source:
                while block:=source.read(8*1024*1024):
                    digest.update(block);whole.update(block);size+=len(block);total+=len(block)
                    if target is not None:target.write(block)
            if size!=part['bytes'] or digest.hexdigest()!=part['sha256']:
                raise ValueError(f'Part checksum/size mismatch: {p.name}')
            print(f'Verified {p.name}',flush=True)
        if total!=manifest['original_bytes'] or whole.hexdigest()!=manifest['original_sha256']:
            raise ValueError('Combined graph checksum/size mismatch')
    finally:
        if target is not None:target.close()
    if target is not None:temporary.rename(output)
    print(f"Verified complete compressed graph: {total:,} bytes; SHA-256 {whole.hexdigest()}")
    if not args.verify_only:print(f'Created {output}. Decompress with: zstd -d "{output}"')


if __name__=='__main__':main()
