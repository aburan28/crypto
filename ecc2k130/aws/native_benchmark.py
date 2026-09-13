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
import shlex
import signal
import tarfile
import time
import uuid

ROOT=Path(__file__).resolve().parents[1]


def bootstrap(bucket,prefix,region,source_sha,transfer_urls=None):
    # All substitutions are generated identifiers or validated AWS names.
    script = f'''#!/bin/bash
set -euo pipefail
systemd-run --unit=native-benchmark-deadline --on-active=55m /sbin/shutdown -h now
mkdir -p /opt/native/src /opt/native/results
exec > /opt/native/results/bootstrap.log 2>&1
export AWS_DEFAULT_REGION={region}
finish() {{
    rc=$?
    trap - EXIT
    set +e
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
    if transfer_urls is not None:
        # Presigned single-object capabilities expire after one hour. The VM
        # receives no long-lived credential or IAM profile. Never log URLs.
        commands={
            f'aws s3 cp s3://{bucket}/{prefix}/source.tgz /opt/native/source.tgz --only-show-errors':
                'curl --fail --silent --show-error --retry 2 --max-time 180 '
                + shlex.quote(transfer_urls['source'])+' -o /opt/native/source.tgz',
            f'aws s3 cp /opt/native/results.tgz s3://{bucket}/{prefix}/results.tgz --only-show-errors':
                'curl --fail --silent --show-error --retry 2 --max-time 180 '
                '--upload-file /opt/native/results.tgz '+shlex.quote(transfer_urls['results']),
            f'aws s3 cp /opt/native/results/exit-code s3://{bucket}/{prefix}/done --only-show-errors':
                'curl --fail --silent --show-error --retry 2 --max-time 60 '
                '--upload-file /opt/native/results/exit-code '+shlex.quote(transfer_urls['done']),
        }
        for old,new in commands.items():
            if script.count(old)!=1: raise ValueError('transfer command contract changed')
            script=script.replace(old,new)
    return script


def launch_request(template,subnet,user_data,token,require_profile=True,instance_type='g7e.2xlarge'):
    # Whitelist fields: never inherit production UserData, tags, spot/fleet
    # settings, extra disks or a live worker's network interface.
    required=('ImageId','SecurityGroupIds','BlockDeviceMappings')
    if instance_type not in ('g7e.2xlarge','g7e.4xlarge'):
        raise ValueError('only the bounded single-GPU G7e sizes are allowed')
    if require_profile: required+=('IamInstanceProfile',)
    for field in required:
        if not template.get(field): raise ValueError('launch template lacks '+field)
    disks=template['BlockDeviceMappings']
    if len(disks)!=1 or not disks[0].get('Ebs',{}).get('DeleteOnTermination'):
        raise ValueError('expected exactly one disposable root volume')
    if not 75<=disks[0]['Ebs'].get('VolumeSize',0)<=100:
        raise ValueError('root volume must be 75 through 100 GiB')
    tags=[{'Key':'Name','Value':'ecc-native-benchmark-'+token},
          {'Key':'Project','Value':'ecc2k130-benchmark'},
          {'Key':'BenchmarkToken','Value':token}]
    request=dict(ImageId=template['ImageId'],
        SecurityGroupIds=template['SecurityGroupIds'],BlockDeviceMappings=disks,
        SubnetId=subnet,InstanceType=instance_type,MinCount=1,MaxCount=1,ClientToken=token,
        # boto3 performs the base64 encoding for RunInstances.
        UserData=user_data,
        InstanceInitiatedShutdownBehavior='terminate',
        MetadataOptions={'HttpTokens':'required','HttpPutResponseHopLimit':2},
        TagSpecifications=[{'ResourceType':kind,'Tags':tags} for kind in ('instance','volume')])
    if require_profile: request['IamInstanceProfile']=template['IamInstanceProfile']
    return request


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
    parser.add_argument('--launch-template',default='ecc2k130-worker')
    parser.add_argument('--launch-config',type=Path,help='explicit AMI/security-group/root-disk JSON instead of an existing template')
    parser.add_argument('--s3-region',help='region of the existing benchmark bucket; defaults to EC2 region')
    parser.add_argument('--availability-zone',help='select an existing eligible subnet in this zone')
    parser.add_argument('--instance-type',choices=('g7e.2xlarge','g7e.4xlarge'),default='g7e.2xlarge')
    parser.add_argument('--automatic-placement',action='store_true',
                        help='let AWS select capacity in the existing default VPC')
    parser.add_argument('--presigned-transfer',action='store_true',
                        help='use expiring object URLs instead of an instance IAM profile')
    args=parser.parse_args()
    if args.automatic_placement and args.availability_zone:
        parser.error('automatic placement and an explicit availability zone are mutually exclusive')
    if not re.fullmatch(r'[a-z]{2}-[a-z]+-\d',args.region): parser.error('invalid region')
    import boto3
    from botocore.config import Config
    from botocore.exceptions import ClientError
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    config=Config(connect_timeout=30,read_timeout=30,retries={'max_attempts':2,'mode':'standard'})
    session=boto3.Session(region_name=args.region)
    ec2=session.client('ec2',config=config)
    s3=session.client('s3',region_name=args.s3_region or args.region,config=config.merge(Config(signature_version='s3v4')))
    account=session.client('sts',config=config).get_caller_identity()['Account']
    if not re.fullmatch(r'\d{12}',account): raise RuntimeError('invalid account identifier')
    bucket='ecc2k130-'+account
    token=uuid.uuid4().hex;prefix='benchmarks/native/'+token
    if args.launch_config:
        template={'VersionNumber':None,'LaunchTemplateData':json.loads(args.launch_config.read_text())}
    else:
        template=ec2.describe_launch_template_versions(LaunchTemplateName=args.launch_template,Versions=['$Default'])['LaunchTemplateVersions'][0]
    data=template['LaunchTemplateData']
    hardware=ec2.describe_instance_types(InstanceTypes=[args.instance_type])['InstanceTypes'][0]
    if sum(gpu['Count'] for gpu in hardware.get('GpuInfo',{}).get('Gpus',[]))!=1:
        raise RuntimeError('benchmark requires exactly one GPU')
    group=ec2.describe_security_groups(GroupIds=data['SecurityGroupIds'])['SecurityGroups'][0]
    subnets=ec2.describe_subnets(Filters=[{'Name':'vpc-id','Values':[group['VpcId']]},
                                        {'Name':'default-for-az','Values':['true']}])['Subnets']
    zones={item['Location'] for item in ec2.describe_instance_type_offerings(LocationType='availability-zone',
        Filters=[{'Name':'instance-type','Values':[args.instance_type]}])['InstanceTypeOfferings']}
    eligible=sorted((s for s in subnets if s['AvailabilityZone'] in zones and s['MapPublicIpOnLaunch']
                     and (args.availability_zone is None or s['AvailabilityZone']==args.availability_zone)),
                    key=lambda s:s['SubnetId'])
    if not eligible: raise RuntimeError('no existing default subnet offers '+args.instance_type+' in '+args.region)
    source_sha=source_archive(out/'source.tgz')
    transfer_urls=None
    if args.presigned_transfer:
        transfer_urls={name:s3.generate_presigned_url(method,
            Params={'Bucket':bucket,'Key':prefix+'/'+key},ExpiresIn=3600)
            for name,method,key in [('source','get_object','source.tgz'),
                                    ('results','put_object','results.tgz'),('done','put_object','done')]}
    user_data=bootstrap(bucket,prefix,args.region,source_sha,transfer_urls)
    request=launch_request(data,eligible[0]['SubnetId'],user_data,token,
                           require_profile=not args.presigned_transfer,instance_type=args.instance_type)
    if args.automatic_placement:
        vpc=ec2.describe_vpcs(VpcIds=[group['VpcId']])['Vpcs'][0]
        if not vpc.get('IsDefault') or not all(s['MapPublicIpOnLaunch'] for s in subnets):
            raise RuntimeError('automatic placement requires the existing default VPC and public default subnets')
        request.pop('SubnetId')
    # Check all launch permissions before uploading source or allocating a GPU.
    try:
        ec2.run_instances(**request,DryRun=True)
        raise RuntimeError('unexpected successful EC2 DryRun response')
    except ClientError as exc:
        if exc.response['Error']['Code']!='DryRunOperation': raise
    s3.upload_file(str(out/'source.tgz'),bucket,prefix+'/source.tgz')
    receipt=dict(valid=False,region=args.region,instanceType=args.instance_type,token=token,
                 availabilityZone=None if args.automatic_placement else eligible[0]['AvailabilityZone'],
                 automaticPlacement=args.automatic_placement,
                 templateName=None if args.launch_config else args.launch_template,presignedTransfer=args.presigned_transfer,
                 s3Region=args.s3_region or args.region,
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
        launched=response['Instances'][0]
        instance=launched['InstanceId'];receipt['instanceId']=instance
        receipt['availabilityZone']=launched['Placement']['AvailabilityZone'];save()
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
