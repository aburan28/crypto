use aes_discover::cipher::{Aes, Config};
use aes_discover::search::{self, DiscoveryOptions, FrozenSpec, OracleKind, ValidationPlan};
use aes_discover::{block, diagnostics, digest, hex, stream, trail, unhex};
use clap::{Args, Parser, Subcommand, ValueEnum};
use rand::seq::SliceRandom;
use rand::RngCore;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::fs::{self, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    version,
    about = "AES synthetic-input differential discovery laboratory; candidates are not cipher breaks"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Args, Clone)]
struct CipherArgs {
    #[arg(long, default_value_t = 128)]
    key_bits: usize,
    #[arg(long, default_value_t = 3)]
    rounds: usize,
    /// Retain MixColumns in the last round (a different permutation from FIPS).
    #[arg(long)]
    final_mix: bool,
}
impl CipherArgs {
    fn config(&self) -> Config {
        Config {
            key_bits: self.key_bits,
            rounds: self.rounds,
            final_mix_columns: self.final_mix,
        }
    }
}

#[derive(Args)]
struct Output {
    /// Write once; an existing file is an error. Omit for stdout.
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Subcommand)]
enum Command {
    /// Enumerate exact AES S-box DDT, LAT and BCT conventions, including BCT diagonal.
    Tables {
        #[command(flatten)]
        output: Output,
    },
    /// Bounded value-level beam proposals; DDT products are model scores only.
    Trail {
        #[command(flatten)]
        cipher: CipherArgs,
        #[arg(long, default_value = "01000000000000000000000000000000")]
        difference: String,
        #[arg(long, default_value_t = 16)]
        beam: usize,
        #[arg(long, default_value_t = 4)]
        branch: usize,
        #[command(flatten)]
        output: Output,
    },
    /// Search empirical byte/parity projections and freeze the complete validation plan.
    Discover {
        #[command(flatten)]
        cipher: CipherArgs,
        #[arg(long,value_enum,default_value_t=OracleKind::Aes)]
        oracle: OracleKind,
        #[arg(long, default_value_t = 20261017)]
        seed: u64,
        #[arg(long, default_value_t = 512)]
        samples: u64,
        #[arg(long, default_value_t = 2)]
        keys: usize,
        #[arg(long, default_value_t = 48)]
        candidates: usize,
        #[arg(long, default_value_t = 8)]
        top: usize,
        #[arg(long, default_value_t = 20261018)]
        validation_seed: u64,
        #[arg(long, default_value_t = 8192)]
        validation_samples: u64,
        #[arg(long, default_value_t = 8)]
        validation_keys: usize,
        #[arg(long, default_value_t = 0.01)]
        alpha: f64,
        #[arg(long, default_value_t = 10_000_000)]
        max_pairs: u64,
        #[command(flatten)]
        output: Output,
    },
    /// Test exactly a frozen spec; no outcome-dependent predicate/seed/sample changes.
    Validate {
        #[arg(long)]
        spec: PathBuf,
        #[arg(long)]
        expect_sha256: String,
        #[arg(long, default_value_t = 10_000_000)]
        max_pairs: u64,
        #[command(flatten)]
        output: Output,
    },
    /// Reproduce 256-text integrals against a PRP and a one-member sibling mutant.
    Integral {
        #[command(flatten)]
        cipher: CipherArgs,
        #[arg(long, default_value_t = 20261017)]
        seed: u64,
        #[arg(long, default_value_t = 8)]
        keys: usize,
        #[arg(long, default_value_t = 4)]
        structures: usize,
        #[arg(long, default_value_t = 0)]
        active_byte: usize,
        #[command(flatten)]
        output: Output,
    },
    /// Cross-ratio through inversion, actual AES affine map, and matched ablations.
    Quartets {
        #[arg(long, default_value_t = 20261017)]
        seed: u64,
        #[arg(long, default_value_t = 65536)]
        trials: u64,
        #[command(flatten)]
        output: Output,
    },
    /// One-bit key-schedule difference diagnostics for all three key lengths.
    Schedule {
        #[arg(long, default_value_t = 20261017)]
        seed: u64,
        #[arg(long, default_value_t = 8)]
        keys: usize,
        #[arg(long, default_value_t = 0)]
        key_bit: usize,
        #[command(flatten)]
        output: Output,
    },
    /// Fixed-delta forward/inverse-oracle diagnostic, with consistent PRP inverse.
    Boomerang {
        #[command(flatten)]
        cipher: CipherArgs,
        #[arg(long, default_value_t = 20261017)]
        seed: u64,
        #[arg(long, default_value_t = 4096)]
        samples: u64,
        #[arg(long)]
        alpha: String,
        #[arg(long)]
        delta: String,
        #[command(flatten)]
        output: Output,
    },
    /// Export a reproducible key-disjoint JSONL dataset for later learned proposals.
    Dataset {
        #[command(flatten)]
        cipher: CipherArgs,
        #[arg(long, default_value_t = 20261017)]
        seed: u64,
        #[arg(long, value_enum)]
        split: Split,
        #[arg(long, default_value_t = 8)]
        keys: usize,
        #[arg(long, default_value_t = 4096)]
        samples: u64,
        #[arg(long, default_value = "01000000000000000000000000000000")]
        difference: String,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value_t = 10_000_000)]
        max_pairs: u64,
    },
    /// Known-answer/reference interop; table-based, not for production encryption.
    Encrypt {
        #[command(flatten)]
        cipher: CipherArgs,
        #[arg(long)]
        key: String,
        #[arg(long)]
        plaintext: String,
        #[command(flatten)]
        output: Output,
    },
}

#[derive(Clone, Copy, ValueEnum)]
enum Split {
    Train,
    Validation,
    Test,
}
impl Split {
    fn name(self) -> &'static str {
        match self {
            Self::Train => "train",
            Self::Validation => "validation",
            Self::Test => "test",
        }
    }
}

fn open_new(path: &Path) -> Result<fs::File, String> {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| format!("{}: {e}", path.display()))
}

fn emit(value: &impl Serialize, output: Output) -> Result<(), String> {
    let mut data = serde_json::to_vec_pretty(value).map_err(|e| e.to_string())?;
    data.push(b'\n');
    match output.out {
        Some(path) => {
            let mut f = open_new(&path)?;
            f.write_all(&data).map_err(|e| e.to_string())?;
            f.sync_all().map_err(|e| e.to_string())?;
        }
        None => {
            std::io::stdout()
                .write_all(&data)
                .map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

struct DatasetOptions {
    config: Config,
    seed: u64,
    split: Split,
    keys: usize,
    samples: u64,
    difference: String,
    out: PathBuf,
    max_pairs: u64,
}

fn dataset(options: DatasetOptions) -> Result<(), String> {
    let DatasetOptions {
        config,
        seed,
        split,
        keys,
        samples,
        difference,
        out,
        max_pairs,
    } = options;
    config.check()?;
    search::budget(samples, keys, 1, max_pairs)?;
    if samples < 2 || !samples.is_multiple_of(2) {
        return Err("dataset samples per key must be positive and even".into());
    }
    let d = block(&difference)?;
    if d == [0; 16] {
        return Err("dataset input difference must be nonzero".into());
    }
    let sidecar = out.with_extension("manifest.json");
    if out == sidecar {
        return Err("dataset output must not end with .manifest.json".into());
    }
    // Reserve both outputs before generating any examples. A partial failure is
    // retained rather than overwritten or silently reused by a later attempt.
    let raw = open_new(&out)?;
    let manifest_file = open_new(&sidecar)?;
    let mut writer = BufWriter::new(raw);
    let mut hash = Sha256::new();
    let mut fingerprints = Vec::new();
    for k in 0..keys {
        let domain = format!("dataset/{}/keys", split.name());
        let key = search::public_key(seed, &domain, k as u64, config.key_bits);
        let key_id = digest(&key);
        fingerprints.push(key_id.clone());
        let aes = Aes::new(&key, config)?;
        let mut prp = search::RandomPermutation::new(
            seed,
            &format!("dataset/{}/prp", split.name()),
            k as u64,
        );
        let mut rng = stream(
            seed,
            &format!("dataset/{}/plaintexts", split.name()),
            k as u64,
        );
        let mut label_rng = stream(seed, &format!("dataset/{}/labels", split.name()), k as u64);
        let mut labels: Vec<_> = (0..samples).map(|i| i % 2).collect();
        labels.shuffle(&mut label_rng);
        let mut used = HashSet::new();
        for label in labels {
            let (p, q) = loop {
                let mut p = [0; 16];
                rng.fill_bytes(&mut p);
                let q = aes_discover::cipher::xor(p, d);
                if !used.contains(&p) && !used.contains(&q) {
                    used.insert(p);
                    used.insert(q);
                    break (p, q);
                }
            };
            let (a, b) = if label == 1 {
                (aes.encrypt(p), aes.encrypt(q))
            } else {
                (prp.encrypt(p), prp.encrypt(q))
            };
            let row = serde_json::json!({"split":split.name(),"public_key_sha256":key_id,"label_aes":label,
                "ciphertext_left":hex(&a),"ciphertext_right":hex(&b),"ciphertext_difference":hex(&aes_discover::cipher::xor(a,b))});
            let mut bytes = serde_json::to_vec(&row).map_err(|e| e.to_string())?;
            bytes.push(b'\n');
            hash.update(&bytes);
            writer.write_all(&bytes).map_err(|e| e.to_string())?;
        }
    }
    writer.flush().map_err(|e| e.to_string())?;
    writer.get_ref().sync_all().map_err(|e| e.to_string())?;
    let manifest = serde_json::json!({"schema":"aes-discover/dataset/v1","cipher":config,"seed":seed,"split":split.name(),
        "keys":keys,"public_key_sha256":fingerprints,"samples_per_key":samples,"examples":samples*keys as u64,
        "aes_examples":samples*keys as u64/2,"permutation_examples":samples*keys as u64/2,
        "input_difference":difference,"jsonl_sha256":hex(&hash.finalize()),"prng":"ChaCha20 with SHA-256 domain separation",
        "scope":"fixed public key within each group, disjoint keys and streams across splits; no training or evaluation performed; metadata and labels are never model features"});
    let mut f = manifest_file;
    f.write_all(&serde_json::to_vec_pretty(&manifest).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    f.sync_all().map_err(|e| e.to_string())?;
    Ok(())
}

fn run(cli: Cli) -> Result<(), String> {
    match cli.command {
        Command::Tables { output } => emit(&diagnostics::tables(), output),
        Command::Trail {
            cipher,
            difference,
            beam,
            branch,
            output,
        } => emit(
            &trail::search(cipher.config(), block(&difference)?, beam, branch)?,
            output,
        ),
        Command::Discover {
            cipher,
            oracle,
            seed,
            samples,
            keys,
            candidates,
            top,
            validation_seed,
            validation_samples,
            validation_keys,
            alpha,
            max_pairs,
            output,
        } => {
            let spec = search::discover(DiscoveryOptions {
                config: cipher.config(),
                oracle,
                seed,
                samples,
                keys,
                candidates,
                top,
                max_pairs,
                validation: ValidationPlan {
                    seed: validation_seed,
                    samples_per_key_per_event: validation_samples,
                    keys: validation_keys,
                    alpha,
                    primary_oracle: oracle,
                    stopping_rule: "fixed_samples_no_optional_stopping".into(),
                },
            })?;
            emit(&spec, output)
        }
        Command::Validate {
            spec,
            expect_sha256,
            max_pairs,
            output,
        } => {
            let bytes = fs::read(&spec).map_err(|e| e.to_string())?;
            if digest(&bytes) != expect_sha256 {
                return Err("frozen specification SHA-256 mismatch".into());
            }
            let spec: FrozenSpec = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
            emit(&search::validate(&spec, &bytes, max_pairs)?, output)
        }
        Command::Integral {
            cipher,
            seed,
            keys,
            structures,
            active_byte,
            output,
        } => emit(
            &diagnostics::integral(cipher.config(), seed, keys, structures, active_byte)?,
            output,
        ),
        Command::Quartets {
            seed,
            trials,
            output,
        } => emit(&diagnostics::quartets(seed, trials)?, output),
        Command::Schedule {
            seed,
            keys,
            key_bit,
            output,
        } => emit(&diagnostics::schedule(seed, keys, key_bit)?, output),
        Command::Boomerang {
            cipher,
            seed,
            samples,
            alpha,
            delta,
            output,
        } => emit(
            &diagnostics::boomerang(
                cipher.config(),
                seed,
                samples,
                block(&alpha)?,
                block(&delta)?,
            )?,
            output,
        ),
        Command::Dataset {
            cipher,
            seed,
            split,
            keys,
            samples,
            difference,
            out,
            max_pairs,
        } => dataset(DatasetOptions {
            config: cipher.config(),
            seed,
            split,
            keys,
            samples,
            difference,
            out,
            max_pairs,
        }),
        Command::Encrypt {
            cipher,
            key,
            plaintext,
            output,
        } => {
            let config = cipher.config();
            let aes = Aes::new(&unhex(&key)?, config)?;
            emit(
                &serde_json::json!({"schema":"aes-discover/encryption-fixture/v1","cipher":config,
                "ciphertext":hex(&aes.encrypt(block(&plaintext)?))}),
                output,
            )
        }
    }
}

fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("aes-discover: {e}");
        std::process::exit(2);
    }
}
