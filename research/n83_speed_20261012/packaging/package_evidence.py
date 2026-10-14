"""Archive existing evidence/documents only; no research computation."""
from pathlib import Path
import json
import subprocess
import zipfile

repo = Path('/workspace/repos/crypto')
output = Path('/workspace/report_output')
paths = subprocess.check_output([
    'git', 'ls-files', '-z',
    'research/n83_evaluation_20261012', 'research/n83_speed_20261012',
    'examples/n83_strong_point_eval.rs', 'examples/n83_rho_optimization_eval.rs',
    'examples/n83_rho_speed_eval.rs',
    '.github/workflows/n83-scalar-arithmetic-lean.yml',
    'experiments/koblitz-single-target-n83-20261006',
], cwd=repo).split(b'\0')
bundle = output / 'n83_faster_rho_evidence_bundle.zip'
with zipfile.ZipFile(bundle, 'x', compression=zipfile.ZIP_DEFLATED) as archive:
    for raw in paths:
        if raw:
            relative = raw.decode()
            archive.write(repo / relative, relative)
    archive.write(output / 'unified_curve_basis_report_faster_rho.pdf',
                  'user_reports/unified_curve_basis_report_faster_rho.pdf')
    archive.writestr('BUNDLE_README.txt',
        'Start at research/n83_speed_20261012/README.md.\n'
        'New frozen rho training/selection/holdout, A/A noise, all139 native answer replays.\n'
        'Checked arithmetic lemmas and explicit machine-word verification scope.\n'
        'Both successful and unselected configurations remain in the evidence.\n'
        'Earlier studies, exact legacy fixture and comprehensive70-page PDF are included.\n'
        'Python files render/inspect/package recorded evidence only; no results recomputed.\n'
        'Archived Python drivers were not run for these native comparisons.\n')
with zipfile.ZipFile(bundle) as archive:
    bad = archive.testzip()
    if bad:
        raise RuntimeError(bad)
    receipt = {'archive_files': len(archive.infolist()), 'zip_crc': 'PASS',
               'archive_bytes': bundle.stat().st_size,
               'native_results_recomputed': False}
(output / 'n83_faster_packaging_receipt.json').write_text(
    json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt))
