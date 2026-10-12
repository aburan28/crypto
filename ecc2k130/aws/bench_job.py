"""Run one benchmark job script on a fresh, bounded, single-GPU EC2 host without SSH.

    python3 aws/bench_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains

Launches one on-demand instance from the existing `ecc2k130-worker` launch
template (AMI, security group, root disk; the `meow34` key pair of AGENTS.md
section 9 is attached when the region has it, but nothing here connects to the
host), drives it entirely through user-data: the source tree is fetched from a
presigned S3 URL, the job runs inside nvidia/cuda:13.3.1-devel-ubuntu24.04 with
the tree at /work and /results mounted, the results tarball is uploaded to a
presigned URL, and the instance shuts itself down (terminate on shutdown), with
a systemd deadline as the backstop.  The instance gets no IAM profile and no
long-lived credential.  Results land under --out; the launch receipt is
launch.json beside them.

Like native_benchmark.py this creates no IAM policy, edits no launch template,
touches no fleet and no campaign checkpoint.  The S3 prefix is unique per run.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import tarfile
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
CUDA_IMAGE = 'nvidia/cuda:13.3.1-devel-ubuntu24.04'


def user_data(urls, source_sha, job, deadline_minutes, job_minutes, results_env):
    # All substitutions are generated identifiers, validated paths or presigned URLs.
    q = shlex.quote
    env = ' '.join('-e %s' % q('%s=%s' % kv) for kv in results_env)
    return f'''#!/bin/bash
set -uo pipefail
systemd-run --unit=bench-job-deadline --on-active={int(deadline_minutes)}m /sbin/shutdown -h now
mkdir -p /opt/bench/src /opt/bench/results
exec > /opt/bench/results/bootstrap.log 2>&1
finish() {{
    rc=$?
    trap - EXIT
    set +e
    printf '%s\\n' "$rc" > /opt/bench/results/exit-code
    tar -czf /opt/bench/results.tgz -C /opt/bench/results .
    curl --fail --silent --show-error --retry 3 --max-time 300 --upload-file /opt/bench/results.tgz {q(urls['results'])}
    curl --fail --silent --show-error --retry 3 --max-time 60 --upload-file /opt/bench/results/exit-code {q(urls['done'])}
    shutdown -h now
}}
trap finish EXIT
curl --fail --silent --show-error --retry 3 --max-time 300 {q(urls['source'])} -o /opt/bench/source.tgz
printf '%s  %s\\n' '{source_sha}' /opt/bench/source.tgz | sha256sum -c -
tar -xzf /opt/bench/source.tgz -C /opt/bench/src
for attempt in {{1..90}}; do
    if nvidia-smi -L && docker info >/dev/null 2>&1; then break; fi
    sleep 5
done
nvidia-smi -L
nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit --format=csv,noheader
docker pull {CUDA_IMAGE}
timeout {int(job_minutes) * 60} docker run --rm --gpus all {env} \\
    -v /opt/bench/src:/work -v /opt/bench/results:/results -w /work \\
    {CUDA_IMAGE} bash {q('/work/' + job)}
'''


def source_archive(path, job, extra):
    rev = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(['git', 'status', '--porcelain', '--', '.'], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    with tarfile.open(path, 'w:gz') as archive:
        for name in ('Makefile', 'src', 'include', 'codegen', 'generated', job) + tuple(extra):
            archive.add(ROOT / name, arcname=name,
                        filter=lambda info: None if '__pycache__' in info.name or info.name.endswith('.pyc') else info)
        marker = path.parent / 'SOURCE_REV'
        marker.write_text('%s%s\n' % (rev, ' (uncommitted changes)' if dirty else ''))
        archive.add(marker, arcname='SOURCE_REV')
    return hashlib.sha256(path.read_bytes()).hexdigest(), rev, bool(dirty)


def launch_request(template, subnet, data, token, instance_type, key_name, name):
    required = ('ImageId', 'SecurityGroupIds', 'BlockDeviceMappings')
    for field in required:
        if not template.get(field): raise ValueError('launch template lacks ' + field)
    disks = template['BlockDeviceMappings']
    if len(disks) != 1 or not disks[0].get('Ebs', {}).get('DeleteOnTermination'):
        raise ValueError('expected exactly one disposable root volume')
    if not 75 <= disks[0]['Ebs'].get('VolumeSize', 0) <= 200:
        raise ValueError('root volume must be 75 through 200 GiB')
    tags = [{'Key': 'Name', 'Value': name}, {'Key': 'Project', 'Value': 'ecc2k130-benchmark'},
            {'Key': 'BenchmarkToken', 'Value': token}]
    request = dict(ImageId=template['ImageId'], SecurityGroupIds=template['SecurityGroupIds'],
                   BlockDeviceMappings=disks, SubnetId=subnet, InstanceType=instance_type,
                   MinCount=1, MaxCount=1, ClientToken=token, UserData=data,
                   InstanceInitiatedShutdownBehavior='terminate',
                   MetadataOptions={'HttpTokens': 'required', 'HttpPutResponseHopLimit': 2},
                   TagSpecifications=[{'ResourceType': kind, 'Tags': tags} for kind in ('instance', 'volume')])
    if key_name: request['KeyName'] = key_name
    return request


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--job', required=True, help='job script path relative to ecc2k130/, run inside the CUDA container')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--region', default='us-east-2')
    parser.add_argument('--availability-zone')
    parser.add_argument('--instance-type', default='g7e.2xlarge',
                        choices=('g7e.2xlarge', 'g7e.4xlarge', 'g7.2xlarge', 'g7.4xlarge', 'g6e.xlarge', 'g6e.2xlarge'))
    parser.add_argument('--launch-template', default='ecc2k130-worker')
    parser.add_argument('--key-name', default='meow34', help="EC2 key pair to attach (AGENTS.md section 9); '' for none")
    parser.add_argument('--extra', action='append', default=[], help='extra tree paths to ship (relative to ecc2k130/)')
    parser.add_argument('--env', action='append', default=[], help='KEY=VALUE passed into the container')
    parser.add_argument('--job-minutes', type=int, default=80)
    parser.add_argument('--deadline-minutes', type=int, default=100, help='in-instance hard shutdown')
    parser.add_argument('--wait-minutes', type=int, default=95)
    args = parser.parse_args()
    if not re.fullmatch(r'[a-z]{2}-[a-z]+-\d', args.region): parser.error('invalid region')
    job = args.job.strip('/')
    if not (ROOT / job).is_file(): parser.error('job script not found: ' + str(ROOT / job))
    import boto3
    from botocore.config import Config
    from botocore.exceptions import ClientError, ReadTimeoutError
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    config = Config(connect_timeout=30, read_timeout=30, retries={'max_attempts': 3, 'mode': 'standard'})
    session = boto3.Session(region_name=args.region)
    ec2 = session.client('ec2', config=config)
    launch_ec2 = session.client('ec2', config=config.merge(Config(read_timeout=300)))
    account = session.client('sts', config=config).get_caller_identity()['Account']
    bucket = 'ecc2k130-' + account
    location = session.client('s3', config=config).get_bucket_location(Bucket=bucket)['LocationConstraint'] or 'us-east-1'
    s3 = session.client('s3', region_name=location, config=config.merge(Config(signature_version='s3v4')))
    token = uuid.uuid4().hex; prefix = 'benchmarks/jobs/' + token
    template = ec2.describe_launch_template_versions(LaunchTemplateName=args.launch_template,
                                                     Versions=['$Default'])['LaunchTemplateVersions'][0]
    data = template['LaunchTemplateData']
    hardware = ec2.describe_instance_types(InstanceTypes=[args.instance_type])['InstanceTypes'][0]
    if sum(gpu['Count'] for gpu in hardware.get('GpuInfo', {}).get('Gpus', [])) != 1:
        raise RuntimeError('benchmark requires exactly one GPU')
    key_name = args.key_name or None
    if key_name:
        try: ec2.describe_key_pairs(KeyNames=[key_name])
        except ClientError:
            print('key pair', key_name, 'is not in', args.region, '- launching without one', flush=True); key_name = None
    group = ec2.describe_security_groups(GroupIds=data['SecurityGroupIds'])['SecurityGroups'][0]
    subnets = ec2.describe_subnets(Filters=[{'Name': 'vpc-id', 'Values': [group['VpcId']]},
                                            {'Name': 'default-for-az', 'Values': ['true']}])['Subnets']
    zones = {item['Location'] for item in ec2.describe_instance_type_offerings(LocationType='availability-zone',
             Filters=[{'Name': 'instance-type', 'Values': [args.instance_type]}])['InstanceTypeOfferings']}
    eligible = sorted((s for s in subnets if s['AvailabilityZone'] in zones and s['MapPublicIpOnLaunch']
                       and (args.availability_zone is None or s['AvailabilityZone'] == args.availability_zone)),
                      key=lambda s: s['AvailabilityZone'])
    if not eligible: raise RuntimeError('no default subnet offers %s in %s' % (args.instance_type, args.region))
    source_sha, rev, dirty = source_archive(out / 'source.tgz', job, args.extra)
    s3.upload_file(str(out / 'source.tgz'), bucket, prefix + '/source.tgz')
    urls = {name: s3.generate_presigned_url(method, Params={'Bucket': bucket, 'Key': prefix + '/' + key}, ExpiresIn=4 * 3600)
            for name, method, key in [('source', 'get_object', 'source.tgz'),
                                      ('results', 'put_object', 'results.tgz'), ('done', 'put_object', 'done')]}
    env = [tuple(kv.split('=', 1)) for kv in args.env]
    script = user_data(urls, source_sha, job, args.deadline_minutes, args.job_minutes, env)
    receipt = dict(valid=False, region=args.region, instanceType=args.instance_type, token=token, job=job,
                   templateName=args.launch_template, templateVersion=template['VersionNumber'],
                   keyName=key_name, sourceSha256=source_sha, gitRev=rev, gitDirty=dirty,
                   resultPrefix=prefix, instanceId=None, availabilityZone=None, launchedAt=None, doneAt=None)
    def save(): (out / 'launch.json').write_text(json.dumps(receipt, indent=2) + '\n')
    def stop(signum, frame): raise KeyboardInterrupt('benchmark cancelled')
    signal.signal(signal.SIGTERM, stop)
    save()
    instance = None
    try:
        # Try each eligible zone in turn: capacity for these instance types is per zone.
        last = None
        for subnet in eligible:
            request = launch_request(data, subnet['SubnetId'], script, token + subnet['AvailabilityZone'][-1],
                                     args.instance_type, key_name, 'ecc-bench-job-' + token[:8])
            try:
                response = launch_ec2.run_instances(**request)
            except ClientError as exc:
                code = exc.response['Error']['Code']
                print('launch in', subnet['AvailabilityZone'], 'failed:', code, flush=True)
                last = exc
                if code in ('InsufficientInstanceCapacity', 'Unsupported', 'InvalidParameterValue'): continue
                raise
            launched = response['Instances'][0]
            instance = launched['InstanceId']
            receipt['instanceId'] = instance
            receipt['availabilityZone'] = launched['Placement']['AvailabilityZone']
            receipt['launchedAt'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
            save()
            break
        if instance is None: raise RuntimeError('no zone had capacity') from last
        print('Started', instance, 'in', receipt['availabilityZone'], flush=True)
        deadline = time.monotonic() + args.wait_minutes * 60
        while time.monotonic() < deadline:
            try:
                done = s3.get_object(Bucket=bucket, Key=prefix + '/done')['Body'].read().decode().strip()
            except ClientError as exc:
                if exc.response['Error']['Code'] not in ('NoSuchKey', '404'): raise
            else:
                receipt['doneAt'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
                receipt['exitCode'] = done
                s3.download_file(bucket, prefix + '/results.tgz', str(out / 'results.tgz'))
                with tarfile.open(out / 'results.tgz') as archive:
                    archive.extractall(out / 'results', filter='data')
                receipt['valid'] = done == '0'
                print('job exit code', done, '- results in', out / 'results', flush=True)
                break
            state = ec2.describe_instances(InstanceIds=[instance])['Reservations'][0]['Instances'][0]['State']['Name']
            if state in ('shutting-down', 'terminated'):
                raise RuntimeError('instance stopped before uploading a result')
            print(time.strftime('%H:%M:%S'), 'instance', state, flush=True)
            time.sleep(30)
        else:
            raise TimeoutError('wait deadline reached')
    finally:
        if instance:
            try:
                ec2.terminate_instances(InstanceIds=[instance])
                receipt['terminationRequestedFor'] = [instance]
            except ClientError as exc:
                receipt['terminationError'] = str(exc)
        save()
    return 0 if receipt['valid'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
