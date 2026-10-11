//! AWS EC2 spot through the `aws` CLI. Prices are live
//! (`describe-spot-price-history`, current price per availability zone).

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;

use super::{sort_offers, CatalogEntry, InstState, LaunchRequest, Offer, Provider};
use crate::record::InstanceRef;
use crate::spec::{Placement, Resources};
use crate::util::{iso8601, now, run, run_ok, s};

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct AwsConfig {
    pub regions: Vec<String>,
    pub catalog: Vec<CatalogEntry>,
    /// Instance profile whose role can read and write the store (and nothing else).
    #[serde(default)]
    pub instance_profile: Option<String>,
    #[serde(default)]
    pub security_group_ids: Vec<String>,
    /// Subnet per availability zone; without one the default VPC's subnet is used.
    #[serde(default)]
    pub subnet_by_zone: BTreeMap<String, String>,
    /// AMI per architecture; defaults to the latest Amazon Linux 2023, which ships the `aws` CLI.
    #[serde(default)]
    pub image_x86_64: Option<String>,
    #[serde(default)]
    pub image_aarch64: Option<String>,
    #[serde(default)]
    pub key_name: Option<String>,
    #[serde(default = "default_disk")]
    pub root_volume_gb: u32,
    /// Shell lines run by the bootstrap before it fetches the agent (credentials for a store on another cloud, packages, ...).
    #[serde(default)]
    pub extra_bootstrap: Vec<String>,
    /// AWS CLI profile used by the controller.
    #[serde(default)]
    pub profile: Option<String>,
}

fn default_disk() -> u32 {
    30
}

pub struct Aws {
    pub cfg: AwsConfig,
}

impl Aws {
    fn base(&self, region: &str) -> Vec<String> {
        let mut v = s(&["aws", "--output", "json", "--region", region]);
        if let Some(p) = &self.cfg.profile {
            v.push("--profile".into());
            v.push(p.clone());
        }
        v
    }

    fn image(&self, arch: &str) -> String {
        match arch {
            "aarch64" => self.cfg.image_aarch64.clone().unwrap_or_else(|| {
                "resolve:ssm:/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64".into()
            }),
            _ => self.cfg.image_x86_64.clone().unwrap_or_else(|| {
                "resolve:ssm:/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64".into()
            }),
        }
    }

    /// The `run-instances` arguments, separate from running them so tests can check them.
    pub fn launch_argv(
        &self,
        offer: &Offer,
        req: &LaunchRequest,
        user_data_path: &str,
    ) -> Result<Vec<String>> {
        let mut spot = serde_json::json!({"SpotInstanceType": "one-time", "InstanceInterruptionBehavior": "terminate"});
        if let Some(max) = req.max_hourly_usd {
            spot["MaxPrice"] = Value::String(format!("{max:.4}"));
        }
        let market = serde_json::json!({"MarketType": "spot", "SpotOptions": spot});
        let tags = format!(
            "ResourceType=instance,Tags=[{{Key=Name,Value={}}},{{Key=spotlab:job,Value={}}},{{Key=spotlab:attempt,Value={}}}]",
            req.name, req.job_id, req.attempt
        );
        let mut v = self.base(&offer.region);
        v.extend(s(&["ec2", "run-instances", "--count", "1"]));
        v.extend(["--image-id".into(), self.image(&offer.arch)]);
        v.extend(["--instance-type".into(), offer.instance_type.clone()]);
        v.extend(["--instance-market-options".into(), market.to_string()]);
        v.extend(s(&["--instance-initiated-shutdown-behavior", "terminate"]));
        v.extend(s(&[
            "--metadata-options",
            "HttpTokens=required,HttpEndpoint=enabled",
        ]));
        v.extend(["--user-data".into(), format!("file://{user_data_path}")]);
        v.extend(["--tag-specifications".into(), tags]);
        v.extend([
            "--block-device-mappings".into(),
            format!(
                "[{{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{{\"VolumeSize\":{},\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}}}]",
                self.cfg.root_volume_gb
            ),
        ]);
        match self.cfg.subnet_by_zone.get(&offer.zone) {
            Some(sn) => v.extend(["--subnet-id".into(), sn.clone()]),
            None => v.extend([
                "--placement".into(),
                format!("AvailabilityZone={}", offer.zone),
            ]),
        }
        if !self.cfg.security_group_ids.is_empty() {
            v.push("--security-group-ids".into());
            v.extend(self.cfg.security_group_ids.iter().cloned());
        }
        if let Some(p) = &self.cfg.instance_profile {
            v.extend(["--iam-instance-profile".into(), format!("Name={p}")]);
        }
        if let Some(k) = &self.cfg.key_name {
            v.extend(["--key-name".into(), k.clone()]);
        }
        Ok(v)
    }
}

/// Current spot prices from `describe-spot-price-history` output: the
/// newest entry per (type, zone).
pub fn parse_spot_prices(json: &[u8]) -> Result<BTreeMap<(String, String), f64>> {
    let v: Value = serde_json::from_slice(json).context("parse spot price history")?;
    let mut newest: BTreeMap<(String, String), (String, f64)> = BTreeMap::new();
    for e in v["SpotPriceHistory"]
        .as_array()
        .cloned()
        .unwrap_or_default()
    {
        let (Some(t), Some(z), Some(p), ts) = (
            e["InstanceType"].as_str(),
            e["AvailabilityZone"].as_str(),
            e["SpotPrice"].as_str().and_then(|x| x.parse::<f64>().ok()),
            e["Timestamp"].as_str().unwrap_or("").to_string(),
        ) else {
            continue;
        };
        let k = (t.to_string(), z.to_string());
        match newest.get(&k) {
            Some((old, _)) if old >= &ts => {}
            _ => {
                newest.insert(k, (ts, p));
            }
        }
    }
    Ok(newest.into_iter().map(|(k, (_, p))| (k, p)).collect())
}

impl Provider for Aws {
    fn name(&self) -> &str {
        "aws"
    }

    fn offers(&self, r: &Resources, p: &Placement) -> Result<Vec<Offer>> {
        let fit: Vec<&CatalogEntry> = self.cfg.catalog.iter().filter(|e| e.fits(r, p)).collect();
        if fit.is_empty() {
            return Ok(vec![]);
        }
        let mut out = Vec::new();
        for region in &self.cfg.regions {
            let mut argv = self.base(region);
            argv.extend(s(&[
                "ec2",
                "describe-spot-price-history",
                "--product-descriptions",
                "Linux/UNIX",
            ]));
            argv.extend(["--start-time".into(), iso8601(now())]);
            argv.push("--instance-types".into());
            argv.extend(fit.iter().map(|e| e.instance_type.clone()));
            let json = run_ok(&argv).with_context(|| format!("spot prices in {region}"))?;
            for ((t, z), price) in parse_spot_prices(&json)? {
                let e = fit
                    .iter()
                    .find(|e| e.instance_type == t)
                    .expect("queried type");
                out.push(Offer {
                    provider: "aws".into(),
                    region: region.clone(),
                    zone: z,
                    instance_type: t,
                    vcpus: e.vcpus,
                    memory_gb: e.memory_gb,
                    arch: e.arch.clone(),
                    hourly_usd: price,
                    price_source: "live".into(),
                });
            }
        }
        sort_offers(&mut out);
        Ok(out)
    }

    fn launch(&self, offer: &Offer, req: &LaunchRequest) -> Result<InstanceRef> {
        let path = std::env::temp_dir().join(format!(
            "spotlab-userdata-{}-{}.sh",
            req.job_id, req.attempt
        ));
        std::fs::write(&path, &req.bootstrap)?;
        let argv = self.launch_argv(offer, req, &path.display().to_string())?;
        let out = run(&argv, None);
        std::fs::remove_file(&path).ok();
        let out = out?;
        if !out.ok() {
            bail!("run-instances: {}", out.stderr.trim());
        }
        let v: Value = serde_json::from_slice(&out.stdout)?;
        let id = v["Instances"][0]["InstanceId"]
            .as_str()
            .context("run-instances returned no instance id")?;
        Ok(InstanceRef {
            provider: "aws".into(),
            region: offer.region.clone(),
            zone: offer.zone.clone(),
            instance_type: offer.instance_type.clone(),
            instance_id: id.to_string(),
        })
    }

    fn state(&self, inst: &InstanceRef) -> Result<InstState> {
        let mut argv = self.base(&inst.region);
        argv.extend(s(&[
            "ec2",
            "describe-instances",
            "--instance-ids",
            &inst.instance_id,
        ]));
        let out = run(&argv, None)?;
        if !out.ok() {
            if out.stderr.contains("InvalidInstanceID.NotFound") {
                return Ok(InstState::Gone("not found".into()));
            }
            bail!("describe-instances: {}", out.stderr.trim());
        }
        let v: Value = serde_json::from_slice(&out.stdout)?;
        let i = &v["Reservations"][0]["Instances"][0];
        let state = i["State"]["Name"].as_str().unwrap_or("unknown");
        Ok(match state {
            "pending" => InstState::Pending,
            "running" => InstState::Running,
            other => {
                let why = i["StateReason"]["Message"].as_str().unwrap_or("");
                InstState::Gone(format!("{other}: {why}"))
            }
        })
    }

    fn terminate(&self, inst: &InstanceRef) -> Result<()> {
        let mut argv = self.base(&inst.region);
        argv.extend(s(&[
            "ec2",
            "terminate-instances",
            "--instance-ids",
            &inst.instance_id,
        ]));
        let out = run(&argv, None)?;
        if !out.ok() && !out.stderr.contains("InvalidInstanceID.NotFound") {
            bail!("terminate-instances: {}", out.stderr.trim());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newest_price_per_zone_wins() {
        let j = br#"{"SpotPriceHistory":[
            {"AvailabilityZone":"us-east-1a","InstanceType":"c7a.large","SpotPrice":"0.0400","Timestamp":"2026-10-09T10:00:00+00:00"},
            {"AvailabilityZone":"us-east-1a","InstanceType":"c7a.large","SpotPrice":"0.0300","Timestamp":"2026-10-09T12:00:00+00:00"},
            {"AvailabilityZone":"us-east-1b","InstanceType":"c7a.large","SpotPrice":"0.0500","Timestamp":"2026-10-09T09:00:00+00:00"}]}"#;
        let m = parse_spot_prices(j).unwrap();
        assert_eq!(
            m[&("c7a.large".to_string(), "us-east-1a".to_string())],
            0.03
        );
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn launch_args_are_spot_and_self_terminating() {
        let aws = Aws {
            cfg: AwsConfig {
                regions: vec!["us-east-1".into()],
                instance_profile: Some("spotlab".into()),
                ..Default::default()
            },
        };
        let offer = Offer {
            provider: "aws".into(),
            region: "us-east-1".into(),
            zone: "us-east-1a".into(),
            instance_type: "c7a.large".into(),
            vcpus: 2,
            memory_gb: 4.0,
            arch: "x86_64".into(),
            hourly_usd: 0.03,
            price_source: "live".into(),
        };
        let req = LaunchRequest {
            name: "spotlab-j-1-a1".into(),
            job_id: "j-1".into(),
            attempt: 1,
            bootstrap: String::new(),
            max_run_s: 3600,
            max_hourly_usd: Some(0.05),
            store_url: "s3://b/p".into(),
        };
        let a = aws
            .launch_argv(&offer, &req, "/tmp/u.sh")
            .unwrap()
            .join(" ");
        assert!(a.contains("\"MarketType\":\"spot\""));
        assert!(a.contains("\"MaxPrice\":\"0.0500\""));
        assert!(a.contains("--instance-initiated-shutdown-behavior terminate"));
        assert!(a.contains("--placement AvailabilityZone=us-east-1a"));
        assert!(a.contains("--iam-instance-profile Name=spotlab"));
        assert!(a.contains("al2023-ami-kernel-default-x86_64"));
    }
}
