//! Google Compute Engine Spot VMs through the `gcloud` CLI. Prices come from
//! the catalog in the config (GCE changes Spot prices at most about once a
//! month); an entry without a price is never offered.

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use super::{static_offers, CatalogEntry, InstState, LaunchRequest, Offer, Provider};
use crate::record::InstanceRef;
use crate::spec::{Placement, Resources};
use crate::util::{run, s};

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct GcpConfig {
    pub project: String,
    pub zones: Vec<String>,
    pub catalog: Vec<CatalogEntry>,
    /// Service account whose role can read and write the store.
    #[serde(default)]
    pub service_account: Option<String>,
    /// Image family per architecture; the Debian images ship `gcloud`.
    #[serde(default)]
    pub image_family_x86_64: Option<String>,
    #[serde(default)]
    pub image_family_aarch64: Option<String>,
    #[serde(default)]
    pub image_project: Option<String>,
    #[serde(default = "default_disk")]
    pub boot_disk_gb: u32,
    #[serde(default)]
    pub network: Option<String>,
    #[serde(default)]
    pub subnet: Option<String>,
    #[serde(default)]
    pub extra_bootstrap: Vec<String>,
}

fn default_disk() -> u32 {
    20
}

pub struct Gcp {
    pub cfg: GcpConfig,
}

impl Gcp {
    pub fn launch_argv(
        &self,
        offer: &Offer,
        req: &LaunchRequest,
        script_path: &str,
    ) -> Vec<String> {
        let family = match offer.arch.as_str() {
            "aarch64" => self
                .cfg
                .image_family_aarch64
                .clone()
                .unwrap_or_else(|| "debian-12-arm64".into()),
            _ => self
                .cfg
                .image_family_x86_64
                .clone()
                .unwrap_or_else(|| "debian-12".into()),
        };
        let mut v = s(&[
            "gcloud",
            "compute",
            "instances",
            "create",
            &req.name,
            "--format=json",
            "--quiet",
        ]);
        v.push(format!("--project={}", self.cfg.project));
        v.push(format!("--zone={}", offer.zone));
        v.push(format!("--machine-type={}", offer.instance_type));
        v.extend(s(&[
            "--provisioning-model=SPOT",
            "--instance-termination-action=DELETE",
        ]));
        // GCE deletes the VM itself when this runs out, whatever the agent does.
        v.push(format!("--max-run-duration={}s", req.max_run_s.max(60)));
        v.push(format!("--image-family={family}"));
        v.push(format!(
            "--image-project={}",
            self.cfg
                .image_project
                .clone()
                .unwrap_or_else(|| "debian-cloud".into())
        ));
        v.push(format!("--boot-disk-size={}GB", self.cfg.boot_disk_gb));
        v.push(format!("--metadata-from-file=startup-script={script_path}"));
        v.push(format!(
            "--labels=spotlab-job={},spotlab-attempt={}",
            label(&req.job_id),
            req.attempt
        ));
        v.push("--scopes=cloud-platform".into());
        if let Some(sa) = &self.cfg.service_account {
            v.push(format!("--service-account={sa}"));
        }
        if let Some(n) = &self.cfg.network {
            v.push(format!("--network={n}"));
        }
        if let Some(sn) = &self.cfg.subnet {
            v.push(format!("--subnet={sn}"));
        }
        v
    }

    fn base(&self, verb: &str, inst: &InstanceRef) -> Vec<String> {
        let mut v = s(&[
            "gcloud",
            "compute",
            "instances",
            verb,
            &inst.instance_id,
            "--format=json",
            "--quiet",
        ]);
        v.push(format!("--project={}", self.cfg.project));
        v.push(format!("--zone={}", inst.zone));
        v
    }
}

/// GCE label values: lower case, digits, `-` and `_`, at most 63 characters.
fn label(x: &str) -> String {
    x.to_ascii_lowercase()
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '-'
            }
        })
        .take(63)
        .collect()
}

fn not_found(stderr: &str) -> bool {
    stderr.contains("was not found") || stderr.contains("notFound") || stderr.contains("404")
}

impl Provider for Gcp {
    fn name(&self) -> &str {
        "gcp"
    }

    fn offers(&self, r: &Resources, p: &Placement) -> Result<Vec<Offer>> {
        Ok(static_offers(
            "gcp",
            &self.cfg.zones,
            &self.cfg.catalog,
            r,
            p,
        ))
    }

    fn launch(&self, offer: &Offer, req: &LaunchRequest) -> Result<InstanceRef> {
        let path =
            std::env::temp_dir().join(format!("spotlab-startup-{}-{}.sh", req.job_id, req.attempt));
        std::fs::write(&path, &req.bootstrap)?;
        let out = run(
            &self.launch_argv(offer, req, &path.display().to_string()),
            None,
        );
        std::fs::remove_file(&path).ok();
        let out = out?;
        if !out.ok() {
            bail!("instances create: {}", out.stderr.trim());
        }
        Ok(InstanceRef {
            provider: "gcp".into(),
            region: offer.region.clone(),
            zone: offer.zone.clone(),
            instance_type: offer.instance_type.clone(),
            instance_id: req.name.clone(),
        })
    }

    fn state(&self, inst: &InstanceRef) -> Result<InstState> {
        let out = run(&self.base("describe", inst), None)?;
        if !out.ok() {
            if not_found(&out.stderr) {
                return Ok(InstState::Gone("deleted".into()));
            }
            bail!("instances describe: {}", out.stderr.trim());
        }
        let v: Value = serde_json::from_slice(&out.stdout).context("parse instance")?;
        Ok(match v["status"].as_str().unwrap_or("UNKNOWN") {
            "PROVISIONING" | "STAGING" => InstState::Pending,
            "RUNNING" => InstState::Running,
            other => InstState::Gone(other.to_ascii_lowercase()),
        })
    }

    fn terminate(&self, inst: &InstanceRef) -> Result<()> {
        let out = run(&self.base("delete", inst), None)?;
        if !out.ok() && !not_found(&out.stderr) {
            bail!("instances delete: {}", out.stderr.trim());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn launch_args_are_spot_with_deadline() {
        let g = Gcp {
            cfg: GcpConfig {
                project: "p".into(),
                zones: vec!["us-central1-a".into()],
                ..Default::default()
            },
        };
        let offer = Offer {
            provider: "gcp".into(),
            region: "us-central1".into(),
            zone: "us-central1-a".into(),
            instance_type: "c3d-highcpu-4".into(),
            vcpus: 4,
            memory_gb: 8.0,
            arch: "x86_64".into(),
            hourly_usd: 0.05,
            price_source: "config".into(),
        };
        let req = LaunchRequest {
            name: "spotlab-j-abc-a2".into(),
            job_id: "j-ABC".into(),
            attempt: 2,
            bootstrap: String::new(),
            max_run_s: 7200,
            max_hourly_usd: None,
            store_url: "gs://b/p".into(),
        };
        let a = g.launch_argv(&offer, &req, "/tmp/s.sh").join(" ");
        assert!(a.contains("--provisioning-model=SPOT"));
        assert!(a.contains("--instance-termination-action=DELETE"));
        assert!(a.contains("--max-run-duration=7200s"));
        assert!(a.contains("--labels=spotlab-job=j-abc,spotlab-attempt=2"));
        assert!(a.contains("--image-family=debian-12 "));
    }
}
