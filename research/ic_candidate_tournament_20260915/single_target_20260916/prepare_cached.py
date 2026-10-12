"""Reuse verified builds from the pre-measurement CPU-affinity preparation failure."""
import json
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
import tournament

OLD = ROOT/'runs/round-0006-single'
cache = {}
for manifest_path in [OLD/'source-manifest.json', *sorted((OLD/'source_candidates').glob('*/source-manifest.json'))]:
    directory = manifest_path.parent
    manifest = tournament.read(manifest_path)
    assert all(tournament.digest(directory/'source'/name) == digest for name, digest in manifest.items())
    assert 'Finished `release` profile' in (directory/'build.log').read_text()
    assert (directory/'worker').is_file()
    cache[tournament.objhash(manifest)] = directory


def reuse(source, destination):
    manifest = {str(path.relative_to(source)): tournament.digest(path)
                for path in sorted(source.rglob('*')) if path.is_file()}
    directory = cache[tournament.objhash(manifest)]
    shutil.copytree(directory/'source', destination/'source')
    for name in ('worker', 'source-manifest.json', 'build.log'):
        shutil.copy2(directory/name, destination/name)
    assert tournament.digest(directory/'worker') == tournament.digest(destination/'worker')
    tournament.write(destination/'build-reuse.json', {
        'reason': 'Original prepare failed on CPU affinity before creating fixtures or running trials.',
        'original_build': str(directory), 'source_manifest_sha256': tournament.objhash(manifest),
        'worker_sha256': tournament.digest(destination/'worker'),
        'all_source_files_rehashed': True,
    }, exclusive=True)
    return destination/'worker', manifest


tournament.snapshot_build = reuse
tournament.main()
