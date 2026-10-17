#!/usr/bin/env python3
"""Resume frozen high-K arms after a post-producer import error."""
import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import run
import verify

HERE = Path(__file__).resolve().parent
ARMS = (4096, 4608, 5120)


def resume(binary, root):
    original = run.frozen(binary)
    corrected = json.loads((HERE/'resume_plan.json').read_text())
    first = root/'k3584'
    receipt_path, analysis_path = first/'run_receipt.json', first/'analysis.json'
    if (corrected['schema'] != 'n83.root11-highk-resume-plan/v1' or
            corrected['original_plan_sha256'] != run.digest(HERE/'plan.json') or
            corrected['resume_script_sha256'] != run.digest(Path(__file__)) or
            corrected['raw_root'] != str(root) or
            corrected['resume_k_values'] != list(ARMS) or
            corrected['cap_seconds_per_arm'] != run.CAP or
            corrected['orchestration_failure_sha256'] !=
            run.digest(root/'orchestration_failure.json') or
            corrected['k3584_raw_receipt_sha256'] != run.digest(receipt_path) or
            corrected['k3584_analysis_sha256'] != run.digest(analysis_path)):
        raise ValueError('sealed resume plan or first-arm evidence differs')
    prior = json.loads(analysis_path.read_text())
    checked = verify.verify_one(first, 3584, original)
    if (checked != prior or prior['status'] !=
            'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY'):
        raise ValueError('K3584 independent replay differs')
    for k in ARMS:
        folder = root/f'k{k}'
        if folder.exists():
            raise ValueError(f'K{k} folder already exists; inspect before retrying')
        folder.mkdir()
        env = run.environment(folder)
        command = [str(binary), f'construct:83:1:{k}', str(run.FIXTURE),
                   '29', str(folder/'targets_result.jsonl')]
        run.write_once(folder/'command.json', {
            'schema':'n83.root11-highk-screen-command/v1', 'k':k,
            'argv':command, 'environment':env, 'cap_seconds':run.CAP,
            'plan_sha256':run.digest(HERE/'plan.json'),
            'resume_plan_sha256':run.digest(HERE/'resume_plan.json')})
        print(f'K{k}: starting resumed arm at {time.time():.3f}',flush=True)
        start=time.monotonic()
        rss_before=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        try:
            process=subprocess.run(command,cwd=run.ROOT,env=env,
                                   capture_output=True,text=True,
                                   timeout=run.CAP,check=False)
            status='EXIT_0' if process.returncode==0 else 'PROCESS_FAILED'
            code,stdout,stderr=process.returncode,process.stdout,process.stderr
        except subprocess.TimeoutExpired as error:
            status,code='TIMEOUT',None
            stdout=(error.stdout or b'').decode(errors='replace')
            stderr=(error.stderr or b'').decode(errors='replace')
        wall=time.monotonic()-start
        rss_after=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        (folder/'stdout.jsonl').write_text(stdout)
        (folder/'stderr.txt').write_text(stderr)
        try:summary=json.loads(stdout) if status=='EXIT_0' else None
        except json.JSONDecodeError:status,summary='MALFORMED_SUMMARY',None
        run.write_once(folder/'run_receipt.json', {
            'schema':'n83.root11-highk-screen-run/v1', 'k':k,
            'producer_status':status,'exit_code':code,
            'outer_wall_seconds':wall,
            'children_maxrss_before_raw':rss_before,
            'children_maxrss_after_raw':rss_after,
            'children_maxrss_unit':'bytes' if sys.platform=='darwin' else 'kib',
            'stdout_sha256':run.digest(folder/'stdout.jsonl'),
            'stderr_sha256':run.digest(folder/'stderr.txt'),
            'base_dump_sha256':run.digest(folder/'base_dump.jsonl') if
                (folder/'base_dump.jsonl').is_file() else None,
            'target_results_sha256':run.digest(folder/'targets_result.jsonl') if
                (folder/'targets_result.jsonl').is_file() else None,
            'summary':summary,
            'resume_plan_sha256':run.digest(HERE/'resume_plan.json')})
        print(f'K{k}: {status} after {wall:.3f}s',flush=True)
        if status!='EXIT_0':break
        try:analysis=verify.verify_one(folder,k,original)
        except (OSError,ValueError,KeyError,TypeError,IndexError) as error:
            analysis={'schema':'n83.root11-highk-screen-analysis/v1',
                      'status':'VERIFICATION_FAILED_OR_CENSORED',
                      'k':k,'error':str(error)}
        run.write_once(folder/'analysis.json',analysis)
        print(f'K{k}: {analysis["status"]}',flush=True)
        if analysis['status']!='PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY':break


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    resume(args.binary.resolve(),args.root.resolve())


if __name__=='__main__':main()
