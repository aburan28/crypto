//! `frozen-materialize`: build the frozen source tree of a compact-orbit
//! evaluation (see the library documentation).
//!
//! ```text
//! frozen-materialize materialize --freeze F [--freeze F2 ...] --dest DIR [--receipt FILE] [--root DIR]
//! frozen-materialize snapshot-base --freeze F --commit SHA --basis evaluation_source_commit|source_base_main [--repo DIR]
//! ```

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use frozen_materialize::{materialize, snapshot_base};

/// The checkout this binary was built from:
/// `<root>/research/notes/ecc2k130/compact_frozen_source_replay_20260929/frozen-materialize`.
fn default_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(5)
        .expect("manifest dir is five levels below the repository root")
        .to_path_buf()
}

fn replay_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("manifest dir has a parent")
        .to_path_buf()
}

struct Args {
    freezes: Vec<PathBuf>,
    dest: Option<PathBuf>,
    receipt: Option<PathBuf>,
    root: Option<PathBuf>,
    commit: Option<String>,
    basis: Option<String>,
    repo: Option<PathBuf>,
}

fn parse(mut it: impl Iterator<Item = String>) -> Result<Args, String> {
    let mut a = Args {
        freezes: vec![],
        dest: None,
        receipt: None,
        root: None,
        commit: None,
        basis: None,
        repo: None,
    };
    while let Some(flag) = it.next() {
        let mut value = |name: &str| it.next().ok_or(format!("{name} needs a value"));
        match flag.as_str() {
            "--freeze" => a.freezes.push(PathBuf::from(value("--freeze")?)),
            "--dest" => a.dest = Some(PathBuf::from(value("--dest")?)),
            "--receipt" => a.receipt = Some(PathBuf::from(value("--receipt")?)),
            "--root" => a.root = Some(PathBuf::from(value("--root")?)),
            "--commit" => a.commit = Some(value("--commit")?),
            "--basis" => a.basis = Some(value("--basis")?),
            "--repo" => a.repo = Some(PathBuf::from(value("--repo")?)),
            other => return Err(format!("unknown argument {other:?}")),
        }
    }
    Ok(a)
}

fn run() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let cmd = args
        .next()
        .ok_or("usage: frozen-materialize materialize|snapshot-base ...")?;
    let a = parse(args)?;
    match cmd.as_str() {
        "materialize" => {
            let root = a.root.unwrap_or_else(default_root);
            let freezes: Vec<PathBuf> = a
                .freezes
                .iter()
                .map(|f| std::fs::canonicalize(f).map_err(|e| format!("{}: {e}", f.display())))
                .collect::<Result<_, _>>()?;
            let root =
                std::fs::canonicalize(&root).map_err(|e| format!("{}: {e}", root.display()))?;
            let dest = a.dest.ok_or("--dest is required")?;
            let dest = if dest.is_absolute() {
                dest
            } else {
                std::env::current_dir()
                    .map_err(|e| e.to_string())?
                    .join(dest)
            };
            let receipt = materialize(&root, &freezes, &dest, &replay_dir())?;
            let text = serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())? + "\n";
            if let Some(path) = a.receipt {
                std::fs::write(&path, &text).map_err(|e| format!("{}: {e}", path.display()))?;
            }
            print!("{text}");
            Ok(())
        }
        "snapshot-base" => {
            let repo = a.repo.unwrap_or_else(default_root);
            let freeze = a.freezes.first().ok_or("--freeze is required")?;
            let freeze = std::fs::canonicalize(freeze).map_err(|e| e.to_string())?;
            let root = std::fs::canonicalize(&repo).map_err(|e| e.to_string())?;
            let rel = freeze
                .strip_prefix(&root)
                .map_err(|_| "the freeze is not under --repo".to_string())?
                .to_string_lossy()
                .into_owned();
            let entry = snapshot_base(
                &repo,
                &replay_dir(),
                &rel,
                a.commit.as_deref().ok_or("--commit is required")?,
                a.basis.as_deref().ok_or("--basis is required")?,
            )?;
            println!(
                "{}",
                serde_json::to_string_pretty(&entry).map_err(|e| e.to_string())?
            );
            Ok(())
        }
        other => Err(format!("unknown command {other:?}")),
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("frozen-materialize: {e}");
            ExitCode::FAILURE
        }
    }
}
