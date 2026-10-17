"""Metadata/copy/archive only. Native programs produce and verify the data.

Usage: REPO CAMPAIGN_ROOT FRESH_STAGE FRESH_ARCHIVE
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from package_checkpoint import archive_tree, digest

repo, root, stage, archive=map(Path,sys.argv[1:])
assert not stage.exists() and not archive.exists()
stage.mkdir()
for name in ['source','source-replay','prefix','prefix-replay','derived-source',
             'derived-source-replay','derived-prefix','derived-prefix-replay','build-final']:
    shutil.copytree(root/name,stage/name)
for name in ['freeze.json','source-replay-match.json','prefix-replay-match.json',
             'B-compressed-cache-eviction.json','B-duplicate-integrity.log',
             'survey-before-pause.json','survey-paused-exit.json','survey-paused-status.json']:
    shutil.copy2(root/name,stage/name)
for pattern in ['*-isolation.jsonl','*.stdout.log','render-*-final.log']:
    for p in root.glob(pattern):shutil.copy2(p,stage/p.name)
execution=repo/'research/full_volcano_graphwide_ic_20261017/execution'
shutil.copytree(execution/'full_point_order',stage/'source-code',ignore=shutil.ignore_patterns('evidence','__pycache__'))
(stage/'library-source').mkdir()
for path,name in [('src/cryptanalysis/ic_boundary.rs','ic_boundary.rs'),
                  ('src/cryptanalysis/ecbench/canonical.rs','canonical.rs'),
                  ('Cargo.toml','Cargo.toml'),('Cargo.lock','Cargo.lock')]:
    (stage/'library-source'/name).write_bytes(subprocess.check_output(['git','show',f'2bf7700b73b3a84810aa86cb29837d1359b28322:{path}'],cwd=repo))
(stage/'inputs').mkdir()
prefix=Path('/workspace/isogeny_class_cli_delivery/graphs/partial-prefixes/C.sqlite')
assert digest(prefix)=='468ce8320b6b7deae37198ca3afc07534fa4fe0e950f621e45ad322f23a2c881'
shutil.copy2(prefix,stage/'inputs/C.sqlite')
shutil.copy2('/workspace/original_local_strata_20261017T1500Z/C.json',stage/'inputs/C-local-conductors.json')
reports=repo/'research/full_volcano_graphwide_ic_20261017/reports'
for name in ['run-20261017T1909Z-C-full-point-source','run-20261017T1913Z-C-full-point-source-replay',
             'run-20261017T1911Z-C-full-point-prefix','run-20261017T1913Z-C-full-point-prefix-replay']:
    shutil.copytree(reports/name,stage/'reports'/name)
(stage/'README.txt').write_text('''Original C: complete selected-point-order diagnostic, not full C class.

Primary prefix: 1,869 models, two public one-target workloads/model, two arms.
Each arm verifies all 3,738 targets: 7,476 method records. Separate full replay
matches every non-clock row, fixture and truth field. Source-only primary and
replay are separate four-record populations. Do not pool them or the older
order-7-only diagnostic with the prefix population.

M is the exact order of each selected point, not a maximality or full-group
cyclicity certificate. All its 2/3/7 prime-power components are retained.
EC PH/BSGS and reduced Tate/Miller + finite-field PH/relation algebra use the
same public input. Scalar truth is outside solver arguments. Field images
have their exact prime-power orders checked. EC verification rejects aliases.
Relations include identity/dependent equations; per-component rank and work
are separate. Small digit groups explain dense yields, not a general oracle.

This is not NFS/integer-smoothness IC, a large-k backend, an ecbench-admitted
speedup, or an independent arithmetic proof. CPU4 wrapper journals accompany
every run; PSI is unavailable and the full-prefix replay is contended. All
clocks remain unqualified and are never pooled with replay observations.
Curve-field inversions, sqrt internals, selected scalar/integer work and
allocation remain unpriced. No combined priced cost or asymptotic claim.

Original C remains 1,934,901,248 vertices / 165 strata. The retained prefix
has only 517 locally exact labels / 24 observed strata, 1,352 unknown labels,
and 1,024 expanded vertices. Original A--E completion, the full method matrix,
transported targets and independent replay remain open. The three additional
complete graph datasets and their ongoing bounded survey are separate.

Linux x86-64 inspection: read prefix/meta.json and derived-prefix/summary.json.
Native replay, fresh output required:
  build-final/full-point run inputs/C.sqlite inputs/C-local-conductors.json /fresh/run prefix
  build-final/full-point compare prefix /fresh/run /fresh/receipt.json
  build-final/full-point-summary /fresh/run /fresh/summary
Source-only mode is 'source'. Every command refuses existing output paths.
Run no builds/tests while measuring. See source-code/PROTOCOL.md and BUILD.sh.
Linked library source is 2bf7700b73b3a84810aa86cb29837d1359b28322; exact binary,
compiler, rlib digests and source snapshots are supplied. Byte-reproducible
compilation is not asserted. All eight final report pages were inspected.

The B cache-eviction receipt concerns only a verified reproducible compressed
local copy. Its full raw 9.3-GB graph and original manifests remain retained;
Drive's prior numbered delivery remains separate. No S3 upload is claimed.
''')
lines=[f'{digest(p)}  {p.relative_to(stage).as_posix()}' for p in sorted(stage.rglob('*')) if p.is_file()]
(stage/'MANIFEST.sha256').write_text('\n'.join(lines)+'\n')
archive_tree(stage,archive)
print(json.dumps({'path':str(archive),'bytes':archive.stat().st_size,'sha256':digest(archive)}))
