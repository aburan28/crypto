//! `store` subcommands, mounted into the repository's binaries
//! (`crypto store …`, `ecbench store …`) so there is one way to move large
//! files wherever the work is done.

use std::path::PathBuf;

use clap::{Args, Subcommand};
use serde_json::json;

use super::{CacheState, PushMeta, RemoteState, Store};

#[derive(Args)]
pub struct StoreArgs {
    #[command(subcommand)]
    pub cmd: StoreCmd,
    /// Machine-readable output where a command has a report.
    #[arg(long, global = true)]
    pub json: bool,
}

#[derive(Subcommand)]
pub enum StoreCmd {
    /// Show the index, cache and remote this checkout resolves to.
    Status,
    /// List indexed names (optionally those starting with a prefix).
    Ls { prefix: Option<String> },
    /// Show one entry in full.
    Info {
        /// A name, `store:NAME` or `sha256:DIGEST`.
        spec: String,
    },
    /// Upload a file to the remote and record it in the index.
    Push {
        file: PathBuf,
        /// Logical name, for example `factor-bases/<curve>/<recipe>.bin`.
        #[arg(long)]
        name: String,
        #[arg(long, default_value = "")]
        description: String,
        #[arg(long = "tag")]
        tags: Vec<String>,
        /// How the file was made: the command, or what it came from.
        #[arg(long, default_value = "")]
        provenance: String,
        /// Move an existing name to this content (the old digest is kept).
        #[arg(long)]
        replace: bool,
        /// Do not keep a copy in the local cache.
        #[arg(long)]
        no_cache: bool,
    },
    /// Make an object available locally (downloading if needed) and print its path.
    Pull {
        /// A name, `store:NAME` or `sha256:DIGEST`.
        spec: String,
        /// Also copy it to this path (an existing file is not overwritten).
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Re-hash cached objects; with --remote also check the remote has each.
    Verify {
        names: Vec<String>,
        #[arg(long)]
        remote: bool,
    },
    /// Drop cached objects the index no longer names (dry run unless --apply).
    /// The remote is never touched.
    Gc {
        #[arg(long)]
        apply: bool,
    },
    /// Remove a name from the index. The remote object is kept.
    Rm { name: String },
}

fn mib(bytes: u64) -> String {
    format!("{:.1} MiB", bytes as f64 / (1024.0 * 1024.0))
}

pub fn run(args: &StoreArgs) -> Result<(), String> {
    let mut store = Store::open_default()?;
    match &args.cmd {
        StoreCmd::Status => {
            let remote = store
                .cfg
                .remote
                .as_ref()
                .map(|r| r.display())
                .unwrap_or_else(|| "none".to_string());
            let reachable = match store.remote() {
                Ok(b) => match b.list("sha256/") {
                    Ok(v) => format!("ok ({} objects)", v.len()),
                    Err(e) => format!("error: {e}"),
                },
                Err(_) => "not configured".to_string(),
            };
            let total: u64 = store.index.entries.values().map(|e| e.bytes).sum();
            if args.json {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&json!({
                        "index": store.cfg.index_path,
                        "entries": store.index.entries.len(),
                        "indexed_bytes": total,
                        "cache": store.cfg.cache_dir,
                        "remote": remote,
                        "remote_status": reachable,
                        "endpoint": store.cfg.endpoint,
                        "aws": store.cfg.aws,
                    }))
                    .map_err(|e| e.to_string())?
                );
            } else {
                println!("index    {}", store.cfg.index_path.display());
                println!(
                    "entries  {} ({} indexed)",
                    store.index.entries.len(),
                    mib(total)
                );
                println!("cache    {}", store.cfg.cache_dir.display());
                println!("remote   {remote}");
                println!("remote   {reachable}");
            }
        }
        StoreCmd::Ls { prefix } => {
            let rows: Vec<_> = store
                .index
                .entries
                .iter()
                .filter(|(n, _)| prefix.as_ref().is_none_or(|p| n.starts_with(p.as_str())))
                .collect();
            if args.json {
                let v: Vec<_> = rows
                    .iter()
                    .map(|(n, e)| {
                        json!({"name": n, "sha256": e.sha256, "bytes": e.bytes,
                        "cached": store.cache_path(&e.sha256).exists(), "tags": e.tags})
                    })
                    .collect();
                println!(
                    "{}",
                    serde_json::to_string_pretty(&v).map_err(|e| e.to_string())?
                );
            } else {
                for (n, e) in rows {
                    let cached = if store.cache_path(&e.sha256).exists() {
                        "cached"
                    } else {
                        "remote"
                    };
                    println!("{:>12}  {cached}  {}  {n}", e.bytes, &e.sha256[..12]);
                }
            }
        }
        StoreCmd::Info { spec } => {
            let (name, sha, _) = store.lookup(spec)?;
            let entry = name.as_deref().and_then(|n| store.index.get(n));
            let cached = store.cache_path(&sha).exists();
            println!(
                "{}",
                serde_json::to_string_pretty(&json!({
                    "name": name,
                    "sha256": sha,
                    "entry": entry,
                    "cached": cached,
                    "cache_path": store.cache_path(&sha),
                    "object_key": super::backend::object_key(&sha)?,
                }))
                .map_err(|e| e.to_string())?
            );
        }
        StoreCmd::Push {
            file,
            name,
            description,
            tags,
            provenance,
            replace,
            no_cache,
        } => {
            let r = store.push(
                name,
                file,
                PushMeta {
                    description: description.clone(),
                    tags: tags.clone(),
                    provenance: provenance.clone(),
                    replace: *replace,
                    no_cache: *no_cache,
                },
            )?;
            if args.json {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&json!({
                        "name": r.name, "sha256": r.sha256, "bytes": r.bytes,
                        "uploaded": r.uploaded, "cached": r.cached, "replaced": r.replaced,
                    }))
                    .map_err(|e| e.to_string())?
                );
            } else {
                eprintln!(
                    "{} {} ({}, sha256 {}){}",
                    if r.uploaded {
                        "uploaded"
                    } else {
                        "already in the remote:"
                    },
                    r.name,
                    mib(r.bytes),
                    r.sha256,
                    if r.replaced {
                        "; moved the name, old digest kept as superseded"
                    } else {
                        ""
                    }
                );
                eprintln!(
                    "index updated: {} (commit it)",
                    store.cfg.index_path.display()
                );
            }
        }
        StoreCmd::Pull { spec, out } => {
            let got = store.ensure(spec)?;
            if got.fetched {
                eprintln!(
                    "fetched {} ({}, sha256 {})",
                    got.name.as_deref().unwrap_or(spec),
                    mib(got.bytes),
                    got.sha256
                );
            }
            if let Some(out) = out {
                if std::fs::symlink_metadata(out).is_ok() {
                    return Err(format!("{} already exists; not overwriting", out.display()));
                }
                super::backend::copy_into(&got.path, out)?;
                println!("{}", out.display());
            } else {
                println!("{}", got.path.display());
            }
        }
        StoreCmd::Verify { names, remote } => {
            let rows = store.verify(names, *remote)?;
            let mut bad = 0;
            for r in &rows {
                let cache = match &r.cache {
                    CacheState::Ok => "cache ok".to_string(),
                    CacheState::Missing => "not cached".to_string(),
                    CacheState::Corrupt(why) => {
                        bad += 1;
                        format!("CACHE CORRUPT: {why}")
                    }
                };
                let remote = match &r.remote {
                    None => String::new(),
                    Some(RemoteState::Ok) => ", remote ok".to_string(),
                    Some(RemoteState::Absent) => {
                        bad += 1;
                        ", REMOTE ABSENT".to_string()
                    }
                    Some(RemoteState::Mismatch(why)) => {
                        bad += 1;
                        format!(", REMOTE MISMATCH: {why}")
                    }
                };
                println!("{}: {cache}{remote}", r.name);
            }
            if bad > 0 {
                return Err(format!("{bad} problem(s)"));
            }
            eprintln!("{} entries verified", rows.len());
        }
        StoreCmd::Gc { apply } => {
            let doomed = store.gc(*apply)?;
            let total: u64 = doomed.iter().map(|(_, b)| b).sum();
            for (p, _) in &doomed {
                println!("{}", p.display());
            }
            eprintln!(
                "{} {} file(s), {}",
                if *apply { "removed" } else { "would remove" },
                doomed.len(),
                mib(total)
            );
        }
        StoreCmd::Rm { name } => {
            if store.index.entries.remove(name).is_none() {
                return Err(format!("{name} is not in the index"));
            }
            store.save_index()?;
            eprintln!("removed {name} from the index; the remote object is kept");
        }
    }
    Ok(())
}
