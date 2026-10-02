//! The host capsule: what the machine was, so a number keeps its host.
//!
//! Facts split two ways, as `docs/ic/measurement/README.md` §6 splits
//! them.  **Stable** facts (CPU model and features, topology, NUMA,
//! frequency policy, SMT, kernel, virtualisation, toolchain) are hashed
//! into `env_class_id`; two wall-clock arms must share it.  **Volatile**
//! facts (load, PSI, ticks) are sampled around every run by the runner
//! and never hashed.  Build provenance (the binary's hash, the commit)
//! is recorded beside the class but not in it: a baseline and a
//! candidate are different binaries on the same host class.
//!
//! The class id here is `ECBENV1h…`, not ICMS's `ENV1h…`: the two hash
//! different fact sets, and a shared prefix would invite a join that
//! means nothing.

use std::collections::BTreeMap;
use std::process::Command;

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::cryptanalysis::ecbench::canonical::{sha256_hex, short_id};
use crate::cryptanalysis::ecbench::isolation::{topology, CpuSettings, Topology};

/// One performance level of a heterogeneous CPU (Apple's P and E
/// cores), from `hw.perflevelN.*`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PerfLevel {
    pub name: String,
    pub physical: Option<u32>,
    pub logical: Option<u32>,
}

/// Facts that do not change between two runs on one machine.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StableFacts {
    pub os: String,
    pub arch: String,
    pub kernel_release: Option<String>,
    pub cpu_vendor: Option<String>,
    pub cpu_model: Option<String>,
    /// `family/model/stepping` on x86.
    pub cpu_signature: Option<String>,
    pub microcode: Option<String>,
    /// Features AGENTS.md §10 names, detected at run time.
    pub features: Vec<String>,
    /// SHA-256 of the full sorted flag list (`/proc/cpuinfo`), Linux only.
    pub flags_sha256: Option<String>,
    pub logical_cpus: u32,
    pub physical_cores: Option<u32>,
    pub packages: Option<u32>,
    pub numa_nodes: u32,
    pub topology: Topology,
    /// Heterogeneous cores (macOS): a run cannot choose between them.
    pub perf_levels: Vec<PerfLevel>,
    pub mem_total_kib: Option<u64>,
    /// Frequency policy, SMT, isolation flags per CPU (Linux).
    pub cpu_settings: BTreeMap<u32, CpuSettings>,
    pub transparent_hugepages: Option<String>,
    /// `hypervisor` CPU flag, DMI product, or `kern.hv_vmm_present`.
    pub virtualization: BTreeMap<String, String>,
    /// `rustc -Vv` of the toolchain on PATH.  The binary's own hash is
    /// the stronger fact and sits in [`BuildFacts`].
    pub rustc_vv: Option<String>,
}

/// Where the measured code came from.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct BuildFacts {
    pub binary_path: Option<String>,
    pub binary_sha256: Option<String>,
    pub git_commit: Option<String>,
    /// `git status --porcelain` was nonempty.
    pub git_dirty: Option<bool>,
    /// Built with debug assertions: an unoptimised build times nothing.
    pub debug_assertions: bool,
    pub target_os: String,
    pub target_arch: String,
    /// Compiler flags in the environment that change generated code.
    pub build_env: BTreeMap<String, String>,
}

/// The whole capsule.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct HostCapsule {
    pub schema: String,
    pub env_class_id: String,
    pub env_class_sha256: String,
    pub stable: StableFacts,
    pub build: BuildFacts,
    pub captured_unix_ms: u128,
}

fn run(cmd: &str, args: &[&str]) -> Option<String> {
    let out = Command::new(cmd).args(args).output().ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

#[cfg(target_os = "macos")]
fn sysctl(name: &str) -> Option<String> {
    run("sysctl", &["-n", name]).filter(|s| !s.is_empty())
}

fn uname_release() -> Option<String> {
    #[cfg(unix)]
    {
        // SAFETY: uname fills a zeroed struct.
        unsafe {
            let mut u: libc::utsname = std::mem::zeroed();
            if libc::uname(&mut u) == 0 {
                let c = std::ffi::CStr::from_ptr(u.release.as_ptr());
                return Some(c.to_string_lossy().into_owned());
            }
        }
    }
    None
}

/// Features named in AGENTS.md §10, as the running CPU reports them.
pub fn detected_features() -> Vec<String> {
    let mut f: Vec<&str> = Vec::new();
    #[cfg(target_arch = "x86_64")]
    {
        macro_rules! probe {
            ($($name:tt),*) => {$(
                if std::arch::is_x86_feature_detected!($name) { f.push($name); }
            )*};
        }
        probe!(
            "popcnt",
            "sse4.2",
            "avx",
            "avx2",
            "bmi1",
            "bmi2",
            "adx",
            "fma",
            "pclmulqdq",
            "vpclmulqdq",
            "gfni",
            "aes",
            "vaes",
            "sha",
            "avx512f",
            "avx512bw",
            "avx512dq",
            "avx512vl",
            "avx512ifma",
            "avx512vbmi"
        );
    }
    #[cfg(target_arch = "aarch64")]
    {
        macro_rules! probe {
            ($($name:tt),*) => {$(
                if std::arch::is_aarch64_feature_detected!($name) { f.push($name); }
            )*};
        }
        probe!("neon", "aes", "pmull", "sha2", "sha3", "crc", "lse", "dotprod", "sve", "sve2");
    }
    let mut v: Vec<String> = f.into_iter().map(String::from).collect();
    v.sort();
    v
}

#[cfg(target_os = "linux")]
fn cpuinfo() -> BTreeMap<String, String> {
    // First processor's block: model name, vendor, family, model,
    // stepping, microcode, flags.
    let mut out = BTreeMap::new();
    if let Ok(text) = std::fs::read_to_string("/proc/cpuinfo") {
        for line in text.lines() {
            if line.trim().is_empty() && !out.is_empty() {
                break;
            }
            if let Some((k, v)) = line.split_once(':') {
                out.entry(k.trim().to_string())
                    .or_insert_with(|| v.trim().to_string());
            }
        }
    }
    out
}

fn stable_facts() -> StableFacts {
    let topo = topology();
    let logical = std::thread::available_parallelism()
        .map(|n| n.get() as u32)
        .unwrap_or(1);
    let mut virtualization = BTreeMap::new();
    #[allow(unused_mut)]
    let mut perf_levels = Vec::new();
    #[allow(unused_mut)]
    let mut facts = StableFacts {
        os: std::env::consts::OS.into(),
        arch: std::env::consts::ARCH.into(),
        kernel_release: uname_release(),
        cpu_vendor: None,
        cpu_model: None,
        cpu_signature: None,
        microcode: None,
        features: detected_features(),
        flags_sha256: None,
        logical_cpus: logical,
        physical_cores: None,
        packages: None,
        numa_nodes: topo.nodes.len().max(1) as u32,
        topology: topo.clone(),
        perf_levels: Vec::new(),
        mem_total_kib: None,
        cpu_settings: BTreeMap::new(),
        transparent_hugepages: None,
        virtualization: BTreeMap::new(),
        rustc_vv: run("rustc", &["-Vv"]),
    };
    #[cfg(target_os = "linux")]
    {
        let ci = cpuinfo();
        facts.cpu_vendor = ci.get("vendor_id").or(ci.get("CPU implementer")).cloned();
        facts.cpu_model = ci.get("model name").or(ci.get("CPU part")).cloned();
        if let (Some(fam), Some(m), Some(s)) =
            (ci.get("cpu family"), ci.get("model"), ci.get("stepping"))
        {
            facts.cpu_signature = Some(format!("{fam}/{m}/{s}"));
        }
        facts.microcode = ci.get("microcode").cloned();
        if let Some(flags) = ci.get("flags").or(ci.get("Features")) {
            let mut fl: Vec<&str> = flags.split_whitespace().collect();
            fl.sort_unstable();
            if fl.contains(&"hypervisor") {
                virtualization.insert("cpu_flag".into(), "hypervisor".into());
            }
            facts.flags_sha256 = Some(sha256_hex(fl.join(" ").as_bytes()));
        }
        let cores: std::collections::BTreeSet<(Option<u32>, Option<u32>)> =
            topo.cpus.iter().map(|c| (c.package, c.core_id)).collect();
        facts.physical_cores = (!topo.cpus.is_empty()).then_some(cores.len() as u32);
        let pkgs: std::collections::BTreeSet<Option<u32>> =
            topo.cpus.iter().map(|c| c.package).collect();
        facts.packages = (!topo.cpus.is_empty()).then_some(pkgs.len() as u32);
        facts.mem_total_kib = std::fs::read_to_string("/proc/meminfo").ok().and_then(|m| {
            m.lines()
                .find(|l| l.starts_with("MemTotal:"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse().ok())
        });
        facts.cpu_settings = topo
            .cpus
            .iter()
            .map(|c| {
                (
                    c.cpu,
                    crate::cryptanalysis::ecbench::isolation::cpu_settings(c.cpu),
                )
            })
            .collect();
        facts.transparent_hugepages =
            std::fs::read_to_string("/sys/kernel/mm/transparent_hugepage/enabled")
                .ok()
                .map(|s| s.trim().to_string());
        for (k, path) in [
            ("dmi_product", "/sys/class/dmi/id/product_name"),
            ("dmi_vendor", "/sys/class/dmi/id/sys_vendor"),
        ] {
            if let Ok(v) = std::fs::read_to_string(path) {
                virtualization.insert(k.into(), v.trim().into());
            }
        }
        if let Some(v) = run("systemd-detect-virt", &["--vm"]) {
            virtualization.insert("systemd_detect_virt_vm".into(), v);
        }
    }
    #[cfg(target_os = "macos")]
    {
        facts.cpu_model = sysctl("machdep.cpu.brand_string");
        facts.cpu_vendor = sysctl("machdep.cpu.vendor");
        facts.physical_cores = sysctl("hw.physicalcpu").and_then(|s| s.parse().ok());
        facts.packages = sysctl("hw.packages").and_then(|s| s.parse().ok());
        facts.mem_total_kib = sysctl("hw.memsize")
            .and_then(|s| s.parse::<u64>().ok())
            .map(|b| b / 1024);
        let levels: u32 = sysctl("hw.nperflevels")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        for l in 0..levels {
            perf_levels.push(PerfLevel {
                name: sysctl(&format!("hw.perflevel{l}.name")).unwrap_or_default(),
                physical: sysctl(&format!("hw.perflevel{l}.physicalcpu"))
                    .and_then(|s| s.parse().ok()),
                logical: sysctl(&format!("hw.perflevel{l}.logicalcpu"))
                    .and_then(|s| s.parse().ok()),
            });
        }
        if let Some(v) = sysctl("kern.hv_vmm_present") {
            virtualization.insert("kern_hv_vmm_present".into(), v);
        }
    }
    facts.virtualization = virtualization;
    facts.perf_levels = perf_levels;
    facts
}

/// Whether the capsule shows a virtual machine.  `None` when nothing
/// said either way, which never passes a bare-metal gate.
pub fn is_virtual(s: &StableFacts) -> Option<bool> {
    if s.virtualization.contains_key("cpu_flag") {
        return Some(true);
    }
    if let Some(v) = s.virtualization.get("systemd_detect_virt_vm") {
        return Some(v != "none");
    }
    if let Some(v) = s.virtualization.get("kern_hv_vmm_present") {
        return Some(v == "1");
    }
    None
}

fn build_facts() -> BuildFacts {
    let exe = std::env::current_exe().ok();
    let binary_sha256 = exe
        .as_ref()
        .and_then(|p| std::fs::read(p).ok())
        .map(|b| sha256_hex(&b));
    let git_commit = run("git", &["rev-parse", "HEAD"]);
    let git_dirty =
        run("git", &["status", "--porcelain", "--untracked-files=no"]).map(|s| !s.is_empty());
    let build_env = [
        "RUSTFLAGS",
        "CARGO_ENCODED_RUSTFLAGS",
        "CARGO_PROFILE_RELEASE_LTO",
        "CC",
        "CFLAGS",
    ]
    .iter()
    .filter_map(|k| std::env::var(k).ok().map(|v| (k.to_string(), v)))
    .collect();
    BuildFacts {
        binary_path: exe.map(|p| p.display().to_string()),
        binary_sha256,
        git_commit,
        git_dirty,
        debug_assertions: cfg!(debug_assertions),
        target_os: std::env::consts::OS.into(),
        target_arch: std::env::consts::ARCH.into(),
        build_env,
    }
}

/// The stable facts that define the class.  Topology enters as its
/// shape (CPU → core, package, node, siblings), settings as values.
fn class_view(s: &StableFacts) -> Value {
    json!({
        "schema": "ecbench.env_class/v1",
        "os": s.os,
        "arch": s.arch,
        "kernel_release": s.kernel_release,
        "cpu_vendor": s.cpu_vendor,
        "cpu_model": s.cpu_model,
        "cpu_signature": s.cpu_signature,
        "microcode": s.microcode,
        "features": s.features,
        "flags_sha256": s.flags_sha256,
        "logical_cpus": s.logical_cpus,
        "physical_cores": s.physical_cores,
        "packages": s.packages,
        "numa_nodes": s.numa_nodes,
        "topology": serde_json::to_value(&s.topology).unwrap_or(Value::Null),
        "perf_levels": serde_json::to_value(&s.perf_levels).unwrap_or(Value::Null),
        "cpu_settings": serde_json::to_value(&s.cpu_settings).unwrap_or(Value::Null),
        "transparent_hugepages": s.transparent_hugepages,
        "virtualization": s.virtualization,
        "rustc_vv": s.rustc_vv,
    })
}

/// Capture the capsule now.
pub fn capture() -> Result<HostCapsule, String> {
    let stable = stable_facts();
    // NUMA distances and memory sizes are integers; the topology carries
    // no float, so the class view is hashable as is.
    let (env_class_id, env_class_sha256) = short_id("ECBENV1", &class_view(&stable))?;
    Ok(HostCapsule {
        schema: "ecbench.host/v1".into(),
        env_class_id,
        env_class_sha256,
        stable,
        build: build_facts(),
        captured_unix_ms: std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_millis())
            .unwrap_or(0),
    })
}

/// Recompute the class id from a capsule's stable facts; an audit
/// compares it with the stored one.
pub fn recompute_class(c: &HostCapsule) -> Result<String, String> {
    Ok(short_id("ECBENV1", &class_view(&c.stable))?.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capsule_class_is_stable_and_recomputable() {
        let a = capture().unwrap();
        assert!(a.env_class_id.starts_with("ECBENV1h"));
        assert_eq!(recompute_class(&a).unwrap(), a.env_class_id);
        assert!(!a.stable.os.is_empty());
    }
}
