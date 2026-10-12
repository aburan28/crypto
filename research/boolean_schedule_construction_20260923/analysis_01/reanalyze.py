#!/usr/bin/env python3
"""Add a corrected analysis without changing the sealed execution bundle."""
import argparse
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,default=HERE/'run_01');parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();root=args.run.resolve();out=args.out.resolve()
    for name,digest in json.loads((root/'manifest.json').read_text())['files'].items():assert sha(root/name)==digest,name
    out.mkdir(parents=True,exist_ok=False)
    for name in ['corrected_analysis.py','reanalyze.py']:shutil.copyfile(HERE/name,out/name)
    metadata={'schema_version':1,'complete':False,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'input_manifest_sha256':sha(root/'manifest.json'),'input_metadata_sha256':sha(root/'metadata.json'),
        'original_analyzer_sha256':sha(root/'analyze.py'),'corrected_analyzer_sha256':sha(out/'corrected_analysis.py'),
        'correction':'Use construction_dispatch_* for the dispatch arms and write analysis only to the new output directory. Cost formulas, guards, thresholds and observations are unchanged.'}
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    spec=importlib.util.spec_from_file_location('corrected_construction',out/'corrected_analysis.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.main(root,out)
    metadata.update(complete=True,ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (out/'manifest.json').write_text(json.dumps({'scope':'Additive correction of analysis only; original execution retained.','files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}},indent=2)+'\n')
if __name__=='__main__':main()
