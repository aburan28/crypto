//! Cloud providers: price discovery, launch, state, terminate.

pub mod aws;
pub mod gcp;
pub mod local;

use anyhow::Result;
use serde::{Deserialize, Serialize};

use crate::record::InstanceRef;
use crate::spec::{Placement, Resources};

/// An instance type the config allows, with what is known about it.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CatalogEntry {
    pub instance_type: String,
    pub vcpus: u32,
    pub memory_gb: f64,
    /// `x86_64` or `aarch64`.
    pub arch: String,
    /// CPU features the type is known to have, lower case (`avx2`, `avx512`, `vpclmulqdq`, `sha`, ...).
    #[serde(default)]
    pub cpu_flags: Vec<String>,
    /// Spot price per hour in USD where no live price API is used (GCP, local).
    /// Leave unset rather than guess: an entry without a price is never offered.
    #[serde(default)]
    pub spot_hourly_usd: Option<f64>,
    #[serde(default)]
    pub note: Option<String>,
}

impl CatalogEntry {
    pub fn fits(&self, r: &Resources, p: &Placement) -> bool {
        self.vcpus >= r.vcpus
            && self.memory_gb + 1e-9 >= r.memory_gb
            && (r.arch == "any" || r.arch == self.arch)
            && r.cpu_flags
                .iter()
                .all(|f| self.cpu_flags.iter().any(|g| g.eq_ignore_ascii_case(f)))
            && (p.instance_types.is_empty()
                || p.instance_types.iter().any(|t| t == &self.instance_type))
    }
}

/// A launchable (type, zone) at a price.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Offer {
    pub provider: String,
    pub region: String,
    pub zone: String,
    pub instance_type: String,
    pub vcpus: u32,
    pub memory_gb: f64,
    pub arch: String,
    pub hourly_usd: f64,
    /// `live` (queried from the provider just now) or `config` (from the catalog).
    pub price_source: String,
}

#[derive(Clone, Debug, PartialEq)]
pub enum InstState {
    Pending,
    Running,
    /// The instance no longer exists or has stopped; the string says why, as far as the provider knows.
    Gone(String),
}

/// What the controller passes to a launch.
pub struct LaunchRequest {
    pub name: String,
    pub job_id: String,
    pub attempt: u32,
    /// Script run as root at first boot (user-data / startup-script).
    pub bootstrap: String,
    /// Hard ceiling on the instance's lifetime.
    pub max_run_s: u64,
    /// Upper bound for the spot price, where the provider supports one.
    pub max_hourly_usd: Option<f64>,
    pub store_url: String,
}

pub trait Provider: Send + Sync {
    fn name(&self) -> &str;
    /// Offers for the catalog entries that fit, cheapest first.
    fn offers(&self, r: &Resources, p: &Placement) -> Result<Vec<Offer>>;
    fn launch(&self, offer: &Offer, req: &LaunchRequest) -> Result<InstanceRef>;
    fn state(&self, inst: &InstanceRef) -> Result<InstState>;
    fn terminate(&self, inst: &InstanceRef) -> Result<()>;
}

/// Offers from a static-price catalog across zones.
pub fn static_offers(
    provider: &str,
    zones: &[String],
    catalog: &[CatalogEntry],
    r: &Resources,
    p: &Placement,
) -> Vec<Offer> {
    let mut out = Vec::new();
    for e in catalog.iter().filter(|e| e.fits(r, p)) {
        let Some(price) = e.spot_hourly_usd else {
            continue;
        };
        for z in zones {
            out.push(Offer {
                provider: provider.into(),
                region: region_of(z),
                zone: z.clone(),
                instance_type: e.instance_type.clone(),
                vcpus: e.vcpus,
                memory_gb: e.memory_gb,
                arch: e.arch.clone(),
                hourly_usd: price,
                price_source: "config".into(),
            });
        }
    }
    sort_offers(&mut out);
    out
}

pub fn sort_offers(v: &mut [Offer]) {
    v.sort_by(|a, b| {
        a.hourly_usd
            .partial_cmp(&b.hourly_usd)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.zone.cmp(&b.zone))
            .then(a.instance_type.cmp(&b.instance_type))
    });
}

/// `us-central1-a` -> `us-central1`; `us-east-1a` -> `us-east-1`.
pub fn region_of(zone: &str) -> String {
    if let Some(i) = zone.rfind('-') {
        let tail = &zone[i + 1..];
        if tail.len() == 1 && tail.chars().all(|c| c.is_ascii_lowercase()) {
            return zone[..i].to_string();
        }
    }
    let t = zone.trim_end_matches(|c: char| c.is_ascii_lowercase());
    if t.len() < zone.len() && t.ends_with(|c: char| c.is_ascii_digit()) {
        return t.to_string();
    }
    zone.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regions() {
        assert_eq!(region_of("us-central1-a"), "us-central1");
        assert_eq!(region_of("us-east-1a"), "us-east-1");
        assert_eq!(region_of("local"), "local");
    }

    #[test]
    fn fit_rules() {
        let e = CatalogEntry {
            instance_type: "c7a.xlarge".into(),
            vcpus: 4,
            memory_gb: 8.0,
            arch: "x86_64".into(),
            cpu_flags: vec!["avx512".into(), "vpclmulqdq".into()],
            spot_hourly_usd: None,
            note: None,
        };
        let mut r = Resources {
            vcpus: 4,
            memory_gb: 8.0,
            arch: "any".into(),
            cpu_flags: vec!["AVX512".into()],
        };
        let p = Placement::default();
        assert!(e.fits(&r, &p));
        r.arch = "aarch64".into();
        assert!(!e.fits(&r, &p));
        r.arch = "x86_64".into();
        r.vcpus = 8;
        assert!(!e.fits(&r, &p));
    }
}
