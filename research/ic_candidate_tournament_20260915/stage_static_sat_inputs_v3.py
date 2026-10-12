"""Stage admitted native bytes and reusable geometry; never execute a solver."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

from identity import canonical
from oracle import require
from static_sat_assets_v3 import freeze_assets
from static_sat_inputs_v3 import native_admission

PUBLISHED_SHA256 = '039e535fa73e5c794aeff666bc02989076f04abbf751b5f609884a2785760678'


def stage(evidence, exporter_source, cms_build_root, output):
    data = Path(evidence).read_bytes()
    require(hashlib.sha256(data).hexdigest() == PUBLISHED_SHA256,
            'native staging requires the accepted diagnostic evidence archive')
    wanted = {'f5/stdout.json', 'f5/source-manifest.json', 'f5/root-source.tar.gz',
              'sat-v2/build/exporter', 'sat-v2/build/build-record.json',
              'sat-v2/build/build.log', 'sat-v2/build/build-exit.json',
              'sat-v2/cms-executable', 'sat-v2/cms-build-receipt.json'}
    old = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for item in tar:
            if item.name in wanted:
                require(item.isfile() and item.name not in old, 'unsafe native source archive member')
                old[item.name] = tar.extractfile(item).read()
    require(set(old) == wanted, 'accepted native evidence missing')
    report = json.loads(old['f5/stdout.json'])
    fixture = dict(report['fixture'], targets=[], target_seeds=[],
                   target_scalar_constructed=False)
    base = dict(factor_base=[[int(x), int(y)] for x, y in report['factor_base']],
                columns=29, column_convention='cofactor')
    files = {'fixture.json': canonical(fixture)+b'\n',
             'base.json': canonical(base)+b'\n',
             'bin/exporter': old['sat-v2/build/exporter'],
             'bin/cms': old['sat-v2/cms-executable'],
             'exporter/source.rs': Path(exporter_source).read_bytes(),
             'exporter/build-record.json': old['sat-v2/build/build-record.json'],
             'exporter/build.log': old['sat-v2/build/build.log'],
             'exporter/build-exit.json': old['sat-v2/build/build-exit.json'],
             'rust/source-manifest.json': old['f5/source-manifest.json'],
             'rust/root-source.tar.gz': old['f5/root-source.tar.gz'],
             'cms/receipt.json': old['sat-v2/cms-build-receipt.json']}
    for name in ('source', 'cadical', 'cadiback'):
        path = Path(cms_build_root)/(name+'.tar')
        require(path.is_file() and not path.is_symlink(), 'CMS build source archive missing or symlinked')
        files['cms/'+name+'.tar'] = path.read_bytes()
    native_admission(files)
    return freeze_assets(files, {'bin/exporter', 'bin/cms'}, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--published-evidence', type=Path, required=True)
    parser.add_argument('--exporter-source', type=Path, required=True)
    parser.add_argument('--cms-build-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    _, seal = stage(args.published_evidence, args.exporter_source, args.cms_build_root, args.out)
    print(json.dumps(seal, sort_keys=True))
