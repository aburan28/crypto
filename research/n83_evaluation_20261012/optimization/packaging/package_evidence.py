"""Package existing native evidence and rendered documents; no result analysis."""
from pathlib import Path
import subprocess
import zipfile
import json
repo = Path('/workspace/repos/crypto')
output = Path('/workspace/report_output')
paths = subprocess.check_output([
    'git', 'ls-files', '-z', 'research/n83_evaluation_20261012',
    'examples/n83_strong_point_eval.rs', 'examples/n83_rho_optimization_eval.rs',
    'experiments/koblitz-single-target-n83-20261006'], cwd=repo).split(b'\0')
bundle = output / 'n83_optimized_rho_evidence_bundle.zip'
with zipfile.ZipFile(bundle, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
    for raw in paths:
        if raw:
            relative = raw.decode()
            archive.write(repo / relative, relative)
    archive.write(output / 'unified_curve_basis_report_optimized_rho.pdf',
                  'user_reports/unified_curve_basis_report_optimized_rho.pdf')
    archive.writestr('BUNDLE_README.txt',
        'Start at research/n83_evaluation_20261012/optimization/README.md.\n'
        'Native frozen training, selection, held-out records and independent replay.\n'
        'Earlier default-rho results and archived fixture are historical evidence.\n'
        'Python files render/package recorded data only; archived Python drivers were not run.\n'
        'No native results were recomputed while packaging.\n')
with zipfile.ZipFile(bundle) as archive:
    bad = archive.testzip()
    if bad:
        raise RuntimeError(bad)
    receipt = {'archive_files': len(archive.infolist()), 'zip_crc': 'PASS',
               'archive_bytes': bundle.stat().st_size, 'native_results_recomputed': False}
(output / 'n83_optimized_packaging_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps(receipt))
