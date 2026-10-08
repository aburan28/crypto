"""Audit the last round and assemble a reviewable operation result."""
import difflib
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
LAST = ROOT/'runs/round-0005-batch16'


def read(path):
    return json.loads(path.read_text())


def main():
    os.sched_setaffinity(0,{2,3})
    with (LAST/'operation.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if not (LAST/'decision.json').exists():
            raise RuntimeError('Final tournament stopped without a complete decision')
        with (LAST/'publication.log').open('w') as log:
            subprocess.run([sys.executable,str(ROOT/'report.py'),'--round',str(LAST),
                            '--scoreboard','/home/ubuntu/crypto/docs/index-calculus-scoreboard.html'],
                           stdout=log,stderr=subprocess.STDOUT,check=True)
        decision=read(LAST/'decision.json')
        winner=next(a for a in read(LAST/'candidates.json') if a['id']==decision['winner'])
        (LAST/'winner-config.json').write_text(json.dumps(winner,indent=2)+'\n')
        source=LAST/winner['source_directory']
        rel=Path('src/cryptanalysis/koblitz_index_calculus.rs')
        old=(ROOT/'runs/round-0002/source'/rel).read_text()
        new=(source/rel).read_text()
        (WORK/'WINNER.patch').write_text(''.join(difflib.unified_diff(
            old.splitlines(True),new.splitlines(True),fromfile='a/'+str(rel),tofile='b/'+str(rel))))
        (WORK/'WINNER.json').write_text(json.dumps({'source_root':str(source),
            'worker':str(LAST/winner['binary_relative']),'config':winner['config'],
            'target_count':read(LAST/'contract.json')['target_count'],
            'source_manifest_sha256':winner['source_manifest_sha256'],
            'decision':str(LAST/'decision.json'),'rho_parity':decision.get('rho_parity')},indent=2)+'\n')
        rounds=[]
        for name in ['round-0002','round-0003b','round-0004','round-0005-batch16']:
            path=ROOT/'runs'/name
            d=read(path/'decision.json');m=read(path/'measurements.json');c=read(path/'contract.json')
            row=next(x for x in m['tables']['confirmation'] if x['variant']==d['winner'])
            pair=d.get('winner_over_rho',{}).get('confirmation',{})
            rounds.append({'round':name,'targets':c.get('target_count',1),'winner':d['winner'],
                           'instruction_ratio_to_rho':row['candidate_over_rho'],
                           'instruction_ci95_to_rho':pair.get('ci95'),
                           'native_ratio_to_rho':pair.get('native_wall_candidate_over_baseline'),
                           'native_ci95_to_rho':pair.get('native_wall_ci95'),
                           'rho_parity':d.get('rho_parity',False),'audit':read(path/'audit.json')})
        data={'rounds':rounds,'new_profiled_trials':sum(x['audit']['trial_receipts'] for x in rounds[1:]),
              'batch_screen_trials':36,'scope':'Separate complete cold single-target and 16-target workloads on five small Koblitz cells; existing signed-Frobenius rho reference.'}
        (WORK/'operation-result.json').write_text(json.dumps(data,indent=2)+'\n')
        text=['# Continued IC tournament results','',
              f"Batch parity verdict: **{decision.get('rho_parity',False)}**; selected implementation **{decision['winner']}**.",'',
              f"Completed three new tournaments with **{data['new_profiled_trials']:,} profiled trials**, each paired with a fresh native run. Every listed round passed its independent artifact/correctness audit. The 36-trial batch screen is separate development evidence.",'',
              '## Profiled instruction cost relative to matched rho','',
              '| Round | Targets per cold job | Selected candidate | Instructions / rho | Paired 95% interval |',
              '|---|---:|---|---:|---|']
        for x in rounds:
            text.append(f"| [{x['round']}](../runs/{x['round']}/REPORT.md) | {x['targets']} | {x['winner']} | {x['instruction_ratio_to_rho']:.4f} | {x['instruction_ci95_to_rho'] or 'not calculated in original protocol'} |")
        text += ['', '## Native process time relative to matched rho','',
                 '| Round | Targets per cold job | Selected candidate | Time / rho | Paired 95% interval |',
                 '|---|---:|---|---:|---|']
        for x in rounds:
            ratio=x['native_ratio_to_rho']
            text.append(f"| {x['round']} | {x['targets']} | {x['winner']} | {ratio:.4f} | {x['native_ci95_to_rho']} |" if ratio is not None else
                        f"| {x['round']} | {x['targets']} | {x['winner']} | diagnostic only | old polling-based timing excluded |")
        text += ['', '## Interpretation','',
                 'Every ratio uses a fresh matched rho run in the same round. The 16-target panel charges all setup once to the complete job and solves every target; it is separate from the single-target result. Rho uses the existing per-target solver API on the same constructed curve. Additional cross-target rho optimizations have not been measured here.', '',
                 'The parity rule was declared before measurement: both candidate/rho upper paired 95% limits and every curve-cell ratio must be at most 1.10, in instructions and native time, on confirmation and replay. The full decisions contain the replay evidence.', '',
                 'These are implementation improvements in fixed-compiler Valgrind amd64 guest instructions and matched native wall time. Kernel/device and external-audit work are outside the instruction count. No arithmetic-complexity, broader-family, or cryptographic-size claim follows. The K-instruction rank floor is deliberately weak.', '',
                 'The successful mechanisms are exact arithmetic substitutions, folded pair-table construction, smaller relation batches, and—where selected—cheaper descent initialization/checking. Every returned scalar still passes general worker verification and the independent Python checker.', '',
                 '## Reproduce and review','',
                 '- [Complete winner source/configuration](WINNER.json) and [cumulative source patch](WINNER.patch).',
                 '- [Operation plan](PLAN.md), [single-target successor plan](ROUND4.md), and [batch plan](ROUND5.md).',
                 '- [Operating guide](../OPERATIONS.md) and [controller tests](controller-tests.json).',
                 '- [Arithmetic/table equivalence tests](preflight-next-tests.log) and [installed skill validation](skill-validation.json).',
                 '- [Retained failed build attempt](../runs/round-0003/prepare_failure.json); no measurements came from it.', '',
                 'Production library defaults were not changed. The frozen source, worker, configurations, raw profiles, receipts and replay evidence are retained under each linked round.']
        (WORK/'RESULTS.md').write_text('\n'.join(text)+'\n')
        with (WORK/'next-proposal.json').open('w') as output:
            subprocess.run([sys.executable,str(ROOT/'tournament.py'),'propose','--from-round',str(LAST),
                            '--out',str(WORK/'next-candidates.json')],stdout=output,check=True)
        print(json.dumps({'report':str(WORK/'RESULTS.md'),'trials':data['new_profiled_trials'],
                          'winner':decision['winner'],'batch_parity':decision.get('rho_parity')}),flush=True)


if __name__=='__main__': main()
