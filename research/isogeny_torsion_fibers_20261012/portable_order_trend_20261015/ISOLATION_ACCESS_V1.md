# Isolated timing access check, 2026-10-15

The requested bare-metal L2 timing was checked after the four-order panel.
`isolab` is absent from this environment and `ISOLAB_URL` is unset; no
configured isolab worker or credential was exposed to this session. The
existing AWS CLI identity could read EC2 inventory and service quotas.
No EC2 instance was launched or altered.

Read-only commands and returned values:

```text
aws service-quotas get-service-quota --region REGION --service-code ec2 --quota-code L-1216C47A --query 'Quota.Value' --output text
us-west-2  0.0
us-east-1  0.0
us-east-2  0.0
us-west-1  0.0
eu-west-1  0.0

aws ec2 describe-instance-types --region us-west-2 --instance-types c6i.metal --query 'InstanceTypes[0].[VCpuInfo.DefaultVCpus,Hypervisor,BareMetal]' --output json
[128, null, true]

aws ec2 describe-instances --region us-west-2 --query 'Reservations[].Instances[].[InstanceType,State.Name]' --output json
No running instances; the listed stopped instances belonged to other work.
```

The Standard On-Demand quota in every checked region is below the 128 vCPUs
needed for `c6i.metal`. This is an account access limit. A quota increase or
an existing authorized bare-metal lab worker would enable an L2 attempt,
subject to actual host isolation telemetry and A/A admission. GitHub's
hosted Linux runner is used only for independent-host L1-or-lower replay.
It cannot discharge the isolated-wall gate.
