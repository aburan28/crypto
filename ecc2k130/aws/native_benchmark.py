"""One bounded, isolated G7e benchmark using the existing worker infrastructure.

Uses existing AWS credentials through boto3's normal chain. Does not create IAM
policies, edit a launch template, change a fleet, or use campaign checkpoints.
Results and the source archive live under a unique benchmarks/native/ prefix.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import signal
import tarfile
import time
import uuid

ROOT=Path(__file__).resolve().parents[1]


def bootstrap(bucket,prefix,region,source_sha):
    # All substitutions are generated identifiers or validated AWS names.
    return f'''#!/bin/bash
set -euo pipefail
systemd-run --unit=native-benchmark-deadline --on-active=55m /sbin/shutdown -h now
mkdir -p /opt/native/src /opt/native/results
exec > /opt/native/results/bootstrap.log 2>&1
export AWS_DEFAULT_REGION={region}
finish() {{
    rc=$?
    trap - EXIT
    printf '%s\\n' "$rc" > /opt/native/results/exit-code
    tar -czf /opt/native/results.tgz -C /opt/native/results .
    aws s3 cp /opt/native/results.tgz s3://{bucket}/{prefix}/results.tgz --only-show-errors || true
    aws s3 cp /opt/native/results/exit-code s3://{bucket}/{prefix}/done --only-show-errors || true
    shutdown -h now
}}
trap finish EXIT
aws s3 cp s3://{bucket}/{prefix}/source.tgz /opt/native/source.tgz --only-show-errors
printf '%s  %s\\n' '{source_sha}' /opt/native/source.tgz | sha256sum -c -
tar -xzf /opt/native/source.tgz -C /opt/native/src
for attempt in {{1..60}}; do
    if nvidia-smi -L && docker info >/dev/null 2>&1; then break; fi
    sleep 5
done
nvidia-smi -L
docker info >/dev/null
timeout 2700 docker run --rm --gpus all \\
    -v /opt/native/src:/work -v /opt/native/results:/results -w /work \\
    nvidia/cuda:13.3.1-devel-ubuntu24.04 \\
    bash -c 'set -euo pipefail; apt-get update -qq; DEBIAN_FRONTEND=noninteractive apt-get install -y -qq build-essential python3; python3 codegen/native_candidate_bench.py --out /results/bench'
'''


def launch_request(template,subnet,user_data,token):
    # Whitelist fields: never inherit production UserData, tags, spot/fleet
    # settings, extra disks or a live worker's network interface.
    for field in ('ImageId','IamInstanceProfile','SecurityGroupIds','BlockDeviceMappings'):
        if not template.get(field): raise ValueError('launch template lacks '+field)
    disks=template['BlockDeviceMappings']
    if len(disks)!=1 or not disks[0].get('Ebs',{}).get('DeleteOnTermination'):
        raise ValueError('expected exactly one disposable root volume')
    if not 75<=disks[0]['Ebs'].get('VolumeSize',0)<=100:
        raise ValueError('root volume must be 75 through 100 GiB')
    tags=[{'Key':'Name','Value':'ecc-native-benchmark-'+token},
          {'Key':'Project','Value':'ecc2k130-benchmark'},
          {'Key':'BenchmarkToken','Value':token}]
    return dict(ImageId=template['ImageId'],IamInstanceProfile=template['IamInstanceProfile'],
        SecurityGroupIds=template['SecurityGroupIds'],BlockDeviceMappings=disks,
        SubnetId=subnet,InstanceType='g7e.2xlarge',MinCount=1,MaxCount=1,ClientToken=token,
        # boto3 performs the base64 encoding for RunInstances.
        UserData=user_data,
        InstanceInitiatedShutdownBehavior='terminate',
        MetadataOptions={'HttpTokens':'required','HttpPutResponseHopLimit':2},
        TagSpecifications=[{'ResourceType':kind,'Tags':tags} for kind in ('instance','volume')])


def source_archive(path):
    with tarfile.open(path,'w:gz') as archive:
        for name in ('Makefile','src','include','codegen'):
            archive.add(ROOT/name,arcname=name,
                        filter=lambda info: None if '__pycache__' in info.name else info)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--region',default='us-west-2')
    args=parser.parse_args()
    if not re.fullmatch(r'[a-z]{2}-[a-z]+-\d',args.region): parser.error('invalid region')
    import boto3
    from botocore.config import Config
    from botocore.exceptions import ClientError
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    config=Config(connect_timeout=10,read_timeout=30,retries={'max_attempts':3,'mode':'standard'})
    session=boto3.Session(region_name=args.region)
    ec2=session.client('ec2',config=config);s3=session.client('s3',config=config)
    account=session.client('sts',config=config).get_caller_identity()['Account']
    if not re.fullmatch(r'\d{12}',account): raise RuntimeError('invalid account identifier')
    bucket='ecc2k130-'+account
    token=uuid.uuid4().hex;prefix='benchmarks/native/'+token
    template=ec2.describe_launch_template_versions(LaunchTemplateName='ecc2k130-worker',Versions=['$Default'])['LaunchTemplateVersions'][0]
    data=template['LaunchTemplateData']
    group=ec2.describe_security_groups(GroupIds=data['SecurityGroupIds'])['SecurityGroups'][0]
    subnets=ec2.describe_subnets(Filters=[{'Name':'vpc-id','Values':[group['VpcId']]},
                                        {'Name':'default-for-az','Values':['true']}])['Subnets']
    zones={item['Location'] for item in ec2.describe_instance_type_offerings(LocationType='availability-zone',
        Filters=[{'Name':'instance-type','Values':['g7e.2xlarge']}])['InstanceTypeOfferings']}
    eligible=sorted((s for s in subnets if s['AvailabilityZone'] in zones and s['MapPublicIpOnLaunch']),key=lambda s:s['SubnetId'])
    if not eligible: raise RuntimeError('no existing default subnet offers g7e.2xlarge in '+args.region)
    source_sha=source_archive(out/'source.tgz')
    user_data=bootstrap(bucket,prefix,args.region,source_sha)
    request=launch_request(data,eligible[0]['SubnetId'],user_data,token)
    # Check all launch permissions before uploading source or allocating a GPU.
    try:
        ec2.run_instances(**request,DryRun=True)
        raise RuntimeError('unexpected successful EC2 DryRun response')
    except ClientError as exc:
        if exc.response['Error']['Code']!='DryRunOperation': raise
    s3.upload_file(str(out/'source.tgz'),bucket,prefix+'/source.tgz')
    receipt=dict(valid=False,region=args.region,instanceType='g7e.2xlarge',token=token,
                 templateVersion=template['VersionNumber'],sourceSha256=source_sha,
                 resultPrefix=prefix,instanceId=None)
    def save(): (out/'launch.json').write_text(json.dumps(receipt,indent=2)+'\n')
    def stop(signum,frame): raise KeyboardInterrupt('benchmark cancelled')
    signal.signal(signal.SIGTERM,stop)
    try:
        save()
        # ClientToken makes a retried request idempotent. Only this returned
        # instance ID, never a Project-tag search, is used for termination.
        response=ec2.run_instances(**request)
        instance=response['Instances'][0]['InstanceId'];receipt['instanceId']=instance;save()
        print('Started isolated benchmark',instance,flush=True)
        deadline=time.monotonic()+50*60
        while time.monotonic()<deadline:
            try:
                done=s3.get_object(Bucket=bucket,Key=prefix+'/done')['Body'].read().decode().strip()
            except ClientError as exc:
                if exc.response['Error']['Code'] not in ('NoSuchKey','404'): raise
            else:
                s3.download_file(bucket,prefix+'/results.tgz',str(out/'results.tgz'))
                with tarfile.open(out/'results.tgz') as archive:
                    try:
                        member=archive.getmember('./bench/result.json')
                    except KeyError:
                        member=None
                    if member:
                        payload=archive.extractfile(member).read()
                        (out/'result.json').write_bytes(payload)
                        result=json.loads(payload)
                        receipt['valid']=done=='0' and result.get('valid') is True and not result.get('compileOnly')
                        print(json.dumps(result.get('comparisons',{}),indent=2),flush=True)
                    for name in ('./bootstrap.log','./exit-code'):
                        try: (out/Path(name).name).write_bytes(archive.extractfile(name).read())
                        except KeyError: pass
                if not receipt['valid']: raise RuntimeError('remote validation failed; see retained logs')
                break
            state=ec2.describe_instances(InstanceIds=[instance])['Reservations'][0]['Instances'][0]['State']['Name']
            if state in ('shutting-down','terminated'): raise RuntimeError('instance stopped before uploading a result')
            print('Benchmark instance',state,flush=True)
            time.sleep(20)
        else: raise TimeoutError('50-minute benchmark deadline reached')
    finally:
        # Also recover the ID after an interrupted launch response, using this
        # unique idempotency token only. The in-instance deadline is independent.
        instance=receipt.get('instanceId')
        if not instance:
            matches=ec2.describe_instances(Filters=[{'Name':'client-token','Values':[token]}])['Reservations']
            ids=[i['InstanceId'] for r in matches for i in r['Instances']]
        else: ids=[instance]
        if ids:
            ec2.terminate_instances(InstanceIds=ids)
            receipt['terminationRequestedFor']=ids
        save()
    return 0


if __name__=='__main__':
    raise SystemExit(main())
