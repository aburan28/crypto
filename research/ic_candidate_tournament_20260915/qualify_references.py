#!/usr/bin/env python3
"""Run the predeclared development panel through the existing tournament CLI."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from oracle import require
from tournament import read, write

HERE = Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifacts',type=Path,required=True,
                        help='Producer artifacts from the same passing workflow run')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    commands=[]

    def run(command):
        commands.append([str(value) for value in command])
        write(out/'commands.json',commands)
        with (out/f'command-{len(commands)-1}.log').open('w') as log:
            subprocess.run(commands[-1],stdout=log,stderr=subprocess.STDOUT,check=True)

    prepared={}
    for reference in ('both','scaled','pairinv'):
        matches=list(args.artifacts.resolve().glob(f'ic-producer-{reference}-*/ic-producer'))
        require(len(matches)==1,'expected one passing producer artifact for '+reference)
        prepared[reference]=matches[0]
        receipt=read(matches[0]/'preparation.json')
        require(receipt['reference']==reference and receipt['instrumented'], 'wrong prepared source')
        run(['cargo','fetch','--locked','--manifest-path',matches[0]/'source/Cargo.toml'])
    config=dict(solver='pair_table',linear_algebra='tiny_gauss',summands=3,batch_trials=1,max_trials=65536)
    registry=[dict(id='incumbent',config=config)]+[
        dict(id=reference,config=config,source_root=str(prepared[reference]/'source'))
        for reference in ('scaled','pairinv')]
    write(out/'candidates.json',registry,exclusive=True)
    shutil.copy2(HERE/'goal_20260924/reference-qualification/PROTOCOL.md',out/'PROTOCOL.md')
    campaign=out/'tournament'
    run([sys.executable,HERE/'tournament.py','prepare','--qualification',
         '--source-root',prepared['both']/'source','--out',campaign,'--candidates',out/'candidates.json',
         '--cells','17a1,19a0,23a0,23a1,31a0','--holdout-cells','29a1',
         '--profile','pilot','--seed','2026092541','--timeout','180','--max-processes','1400',
         '--qualification-widths','1','2','4','8','16','32',
         '--selection-width','3','--exploration-slots','1','--comparison-kind','factor-base-policy'])
    run([sys.executable,campaign/'evaluator/tournament.py','run','--round',campaign])
    run([sys.executable,campaign/'evaluator/tournament.py','verify','--round',campaign])
    report=read(campaign/'qualification.json')
    require(report['cases']==15 and len(report['cells'])==5 and report['repetitions']==3,
            'declared qualification panel differs')
    require(len(report['table'])==21 and report['promotion_eligible'] is False,
            'missing reference or invalid promotion')
    write(out/'summary.json',report,exclusive=True)
    print(json.dumps(report))


if __name__=='__main__':
    main()
