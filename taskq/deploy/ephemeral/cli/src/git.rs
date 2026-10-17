//! Source snapshot: the newest commit a remote already has, plus a binary patch
//! of everything after it. Same rules as `pushed_base` and `worktree_patch` in
//! `../ephemeral`.

use anyhow::{Context, Result, bail};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const MAX_UNTRACKED: u64 = 5 << 20;

pub struct Repo {
    pub root: PathBuf,
}

impl Repo {
    pub fn discover() -> Result<Repo> {
        let out = Command::new("git")
            .args(["rev-parse", "--show-toplevel"])
            .output()?;
        if !out.status.success() {
            bail!("not in a git repository (use --no-source)");
        }
        Ok(Repo {
            root: PathBuf::from(String::from_utf8(out.stdout)?.trim()),
        })
    }

    fn cmd(&self) -> Command {
        let mut c = Command::new("git");
        c.current_dir(&self.root);
        c
    }

    /// stdout of a git command, trimmed; None when it fails.
    pub fn try_git(&self, args: &[&str]) -> Option<String> {
        let out = self.cmd().args(args).stderr(Stdio::null()).output().ok()?;
        out.status
            .success()
            .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
    }

    pub fn git(&self, args: &[&str]) -> Result<String> {
        self.try_git(args)
            .with_context(|| format!("git {} failed", args.join(" ")))
    }

    fn on_a_remote(&self, commit: &str) -> bool {
        self.try_git(&["for-each-ref", "--contains", commit, "refs/remotes"])
            .is_some_and(|s| !s.is_empty())
    }

    /// `ref` resolved to a commit that some remote has.
    pub fn pushed_ref(&self, r: &str) -> Result<String> {
        let commit = self
            .try_git(&["rev-parse", "--verify", &format!("{r}^{{commit}}")])
            .with_context(|| format!("unknown ref {r}"))?;
        if !self.on_a_remote(&commit) {
            bail!(
                "{r} ({commit}) is not on any remote: push it, or drop --ref to ship the working tree"
            );
        }
        Ok(commit)
    }

    /// The newest commit of HEAD's history that a remote already has.
    pub fn pushed_base(&self) -> Result<String> {
        let unpushed = self.git(&["rev-list", "--topo-order", "HEAD", "--not", "--remotes"])?;
        let Some(oldest) = unpushed.lines().last() else {
            return self.git(&["rev-parse", "HEAD"]);
        };
        let base = self
            .try_git(&["rev-parse", "--verify", "-q", &format!("{oldest}^")])
            .context("no commit of HEAD's history is on a remote: push something first")?;
        if self.on_a_remote(&base) {
            return Ok(base);
        }
        let upstream = self
            .try_git(&["rev-parse", "--abbrev-ref", "origin/HEAD"])
            .unwrap_or_else(|| "origin/main".into());
        self.try_git(&["merge-base", "HEAD", &upstream])
            .context("cannot find a pushed ancestor of HEAD")
    }

    /// Binary diff from `base` to the working tree, built in a copy of the
    /// index so the real one is untouched. Skipped untracked files are named
    /// on stderr, so nothing is dropped silently.
    pub fn worktree_patch(&self, base: &str, untracked: bool, tmp: &Path) -> Result<Vec<u8>> {
        let real = self.git(&["rev-parse", "--git-path", "index"])?;
        let idx = tmp.join("index");
        // A copy keeps the stat cache, so `add` hashes only what changed.
        std::fs::copy(self.root.join(real), &idx).context("copying the git index")?;
        let with_idx = |args: &[&str]| {
            let mut c = self.cmd();
            c.env("GIT_INDEX_FILE", &idx).args(args);
            c
        };
        if !with_idx(&["add", "-u"]).status()?.success() {
            bail!("git add -u failed");
        }
        if untracked {
            let out = with_idx(&["ls-files", "-z", "--others", "--exclude-standard"]).output()?;
            let mut keep = Vec::new();
            for f in out.stdout.split(|b| *b == 0).filter(|f| !f.is_empty()) {
                let name = String::from_utf8_lossy(f);
                let p = format!("/{name}");
                if ["/.venv/", "/venv/", "/node_modules/", "/__pycache__/"]
                    .iter()
                    .any(|d| p.contains(d))
                {
                    continue;
                }
                let size = std::fs::symlink_metadata(self.root.join(&*name))
                    .map(|m| m.len())
                    .unwrap_or(0);
                if size > MAX_UNTRACKED {
                    eprintln!(
                        "ephemeral: not shipping untracked {name} ({} MiB); commit it or use --artifact on the box",
                        size >> 20
                    );
                    continue;
                }
                keep.extend_from_slice(f);
                keep.push(0);
            }
            if !keep.is_empty() {
                let mut child = with_idx(&["add", "--pathspec-from-file=-", "--pathspec-file-nul"])
                    .stdin(Stdio::piped())
                    .spawn()?;
                child.stdin.take().unwrap().write_all(&keep)?;
                if !child.wait()?.success() {
                    bail!("git add of untracked files failed");
                }
            }
        }
        let out = with_idx(&["diff", "--cached", "--binary", base]).output()?;
        if !out.status.success() {
            bail!("git diff failed: {}", String::from_utf8_lossy(&out.stderr));
        }
        Ok(out.stdout)
    }
}
