//! CLI demo for the crypto library.
//!
//! Commands:
//!   hash    <algorithm> <message>      — Compute a hash
//!   ecc     keygen                     — Generate an ECC key pair
//!   ecdsa   sign <privkey-hex> <msg>   — ECDSA sign
//!   ecdsa   verify <pubkey-hex> <msg> <r-hex> <s-hex>
//!   ecdh    <privkey-hex> <pubkey-hex> — ECDH shared secret
//!   rsa     keygen                     — Generate an RSA key pair (slow)
//!   rsa     sign <priv-json> <msg>     — RSA sign
//!   aes     encrypt <key-hex> <nonce-hex> <msg>
//!   chacha  encrypt <key-hex> <nonce-hex> <msg>
//!   hkdf    <ikm-hex> <salt-hex> <info> <length>
//!   pqc     kyber                      — Kyber KEM demo

use clap::{Parser, Subcommand};
use crypto_lib::{
    asymmetric::rsa::{rsa_sign, rsa_verify, RsaKeyPair},
    ecc::{
        curve::CurveParams,
        ecdh::ecdh,
        ecdsa::{sign, verify, EcdsaSignature},
        field::FieldElement,
        keys::{EccKeyPair, EccPublicKey},
        point::Point,
    },
    hash::{blake3_hash, sha256, sha3_256},
    kdf::hkdf::hkdf,
    pqc::kyber::{kyber_decapsulate, kyber_encapsulate, kyber_keygen},
    symmetric::{
        aes::{aes_gcm_decrypt, aes_gcm_encrypt, AesKey},
        chacha20::{chacha20_poly1305_decrypt, chacha20_poly1305_encrypt},
    },
    utils::encoding::{bigint_to_bytes_be, from_hex, to_hex},
};
use num_bigint::BigUint;

#[derive(Parser)]
#[command(name = "crypto", about = "Educational cryptography library demo")]
struct Cli {
    #[command(subcommand)]
    command: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Hash a message with the given algorithm (sha256 | sha3 | blake3)
    Hash { algorithm: String, message: String },
    /// ECC operations (keygen)
    Ecc {
        #[command(subcommand)]
        op: EccOp,
    },
    /// ECDSA sign / verify
    Ecdsa {
        #[command(subcommand)]
        op: EcdsaOp,
    },
    /// ECDH shared secret demo
    Ecdh,
    /// RSA operations
    Rsa {
        #[command(subcommand)]
        op: RsaOp,
    },
    /// AES-GCM encrypt / decrypt demo
    Aes,
    /// ChaCha20-Poly1305 encrypt / decrypt demo
    Chacha,
    /// HKDF key derivation demo
    Hkdf,
    /// Post-quantum Kyber KEM demo
    Pqc,
    /// Run all demos
    Demo,
    /// Cryptanalysis: run cipher-attack techniques (S-box analysis,
    /// boomerang distinguisher, rectangle attack, trail search, …)
    Cryptanalysis {
        #[command(subcommand)]
        op: CryptanalysisOp,
    },
    /// **Visual demos** for every part of the library — symmetric
    /// ciphers, hash constructions, ECC.
    Visual {
        #[command(subcommand)]
        op: VisualOp,
    },
}

#[derive(Subcommand)]
enum VisualOp {
    /// Symmetric ciphers: ECB-penguin, mode dataflows, AES round
    /// trace, ChaCha state, avalanche heat-map.
    Symmetric {
        /// Sub-target: `ecb-penguin`, `modes`, `aes-trace`, `chacha`,
        /// `avalanche`, `all`.
        #[arg(long, default_value = "all")]
        target: String,
    },
    /// Hash functions: Merkle-Damgård / sponge / BLAKE3 tree diagrams,
    /// SHA-256 round trace, hash avalanche.
    Hash {
        #[arg(long, default_value = "all")]
        target: String,
    },
    /// Elliptic curves: toy-curve scatter, point operations,
    /// scalar-multiplication ladder, ECDSA dataflow.
    Ecc {
        #[arg(long, default_value = "all")]
        target: String,
    },
    /// Everything at once — emits the full visual bundle.
    All,
}

#[derive(Subcommand)]
enum CryptanalysisOp {
    /// List all registered ciphers known to the auto-attack runner.
    ListCiphers,
    /// Run every applicable attack against a registered cipher and
    /// emit a Markdown report.
    ///
    /// Example: `cargo run -- cryptanalysis auto --cipher toyspn-2r`
    Auto {
        /// Cipher identifier (use `list-ciphers` to see options).
        #[arg(long)]
        cipher: String,
        /// Number of plaintext pairs for the boomerang.
        #[arg(long, default_value_t = 65_536)]
        boomerang_pairs: usize,
        /// Pool size for the rectangle attack.
        #[arg(long, default_value_t = 2048)]
        rectangle_pool: usize,
        /// Negative log₂ of the differential-trail probability
        /// threshold.  `12` here means threshold = `2⁻¹²`.
        #[arg(long, default_value_t = 12)]
        trail_threshold_neg_log2: u32,
        /// Max number of trails to render.
        #[arg(long, default_value_t = 5)]
        trail_top_k: usize,
    },
    /// Run just the boomerang distinguisher against a registered cipher.
    Boomerang {
        #[arg(long)]
        cipher: String,
        /// α and δ in hex (low byte first), e.g. `0100`.  Defaults to
        /// the cipher's canonical pair.
        #[arg(long)]
        alpha: Option<String>,
        #[arg(long)]
        delta: Option<String>,
        #[arg(long, default_value_t = 65_536)]
        pairs: usize,
    },
    /// Run just the rectangle attack against a registered cipher.
    Rectangle {
        #[arg(long)]
        cipher: String,
        #[arg(long)]
        alpha: Option<String>,
        #[arg(long)]
        delta: Option<String>,
        #[arg(long, default_value_t = 2048)]
        pool: usize,
    },
    /// Print the S-box differential / linear / boomerang report for a
    /// registered cipher.
    Sbox {
        #[arg(long)]
        cipher: String,
    },
    /// Run the full research bench and emit the Markdown report.
    /// (Equivalent to `cargo test --lib --release bench_demo --
    /// --ignored --nocapture` but accessible from the CLI.)
    Bench,
    /// Auto-run every applicable hash-function attack (length
    /// extension, birthday collision, Joux multicollision,
    /// differential bias) against the named hash.  Targets: md4, md5, sha1.
    HashAuto {
        /// Hash identifier: `md4`, `md5`, or `sha1`.
        #[arg(long)]
        hash: String,
    },
    /// Length-extension attack demo against a Merkle-Damgård hash.
    LengthExtension {
        /// Hash identifier: `md4`, `md5`, or `sha1`.
        #[arg(long)]
        hash: String,
        /// Length (in bytes) of the unknown secret prefix.
        #[arg(long)]
        secret_len: usize,
        /// Original message body (hex-encoded).
        #[arg(long)]
        message_hex: String,
        /// Suffix to append (hex-encoded).
        #[arg(long)]
        suffix_hex: String,
        /// Digest of `secret || message_hex` (hex-encoded).
        #[arg(long)]
        digest_hex: String,
    },
    /// **AES related-key attack** — propagate a key-difference
    /// through the schedule + run avalanche + run a 4-round local
    /// collision demo.  Emits a Markdown report.
    AesRelatedKey {
        /// Key size in bits: 128 or 256.
        #[arg(long, default_value_t = 128)]
        key_bits: u32,
        /// Key-difference ΔK in hex (length = key_bits / 8 bytes).
        /// Defaults to a single-active-byte difference at position 0.
        #[arg(long)]
        delta_k_hex: Option<String>,
        /// Number of avalanche trials.
        #[arg(long, default_value_t = 1024)]
        avalanche_trials: usize,
        /// Number of rounds to study (default = full).  Sets the
        /// avalanche and local-collision round count.
        #[arg(long, default_value_t = 0)]
        rounds: usize,
    },
    /// **AES visual demos** — emit ASCII diagrams for the major
    /// AES attack visualizations: truncated-diff trail, 3-round
    /// integral distinguisher, DFA single-byte-fault propagation,
    /// related-key schedule diffusion.
    AesVisualDemo {
        /// Which demo to run: `truncated-diff`, `integral`, `dfa`,
        /// `key-schedule`, or `all` (default).
        #[arg(long, default_value = "all")]
        demo: String,
    },
    /// **Visual demos for every non-AES attack** — DDT/LAT heat-maps
    /// on S-boxes, Walsh-Hadamard spectra, Pollard ρ trajectories,
    /// HNP recovery curves, Bleichenbacher bias histograms,
    /// length-extension diagrams, Joux multicollision trees, j=0
    /// twist factorisations, birthday-paradox curves.
    VisualAll {
        /// Optional single-demo target: `sbox-ddt`, `sbox-lat`,
        /// `walsh`, `pollard-rho`, `hnp`, `bleichenbacher`,
        /// `length-extension`, `joux`, `j0-twists`, `birthday`,
        /// or `all` (default).
        #[arg(long, default_value = "all")]
        target: String,
    },
}

#[derive(Subcommand)]
enum EccOp {
    Keygen,
}

#[derive(Subcommand)]
enum EcdsaOp {
    Sign {
        message: String,
    },
    Verify {
        message: String,
        pubkey_hex: String,
        r_hex: String,
        s_hex: String,
    },
}

#[derive(Subcommand)]
enum RsaOp {
    /// Generate a 1024-bit RSA key pair (faster for demo purposes)
    Keygen,
    /// Sign + verify demo
    Demo,
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Cmd::Hash { algorithm, message } => cmd_hash(&algorithm, &message),
        Cmd::Ecc { op: EccOp::Keygen } => cmd_ecc_keygen(),
        Cmd::Ecdsa { op } => cmd_ecdsa(op),
        Cmd::Ecdh => cmd_ecdh(),
        Cmd::Rsa { op } => cmd_rsa(op),
        Cmd::Aes => cmd_aes(),
        Cmd::Chacha => cmd_chacha(),
        Cmd::Hkdf => cmd_hkdf(),
        Cmd::Pqc => cmd_pqc(),
        Cmd::Demo => cmd_demo(),
        Cmd::Cryptanalysis { op } => cmd_cryptanalysis(op),
        Cmd::Visual { op } => cmd_visual(op),
    }
}

// ── Cryptanalysis subcommands ─────────────────────────────────────────────────

fn cmd_cryptanalysis(op: CryptanalysisOp) {
    use crypto_lib::cryptanalysis::auto_attack::{auto_attack_with, AutoAttackOptions};
    use crypto_lib::cryptanalysis::boomerang::{boomerang_distinguisher, rectangle_attack};
    use crypto_lib::cryptanalysis::cipher_registry::{list_ciphers, RegisteredCipher};
    use crypto_lib::cryptanalysis::research_bench::run_full_bench;
    match op {
        CryptanalysisOp::ListCiphers => {
            println!("Registered ciphers:");
            for name in list_ciphers() {
                if let Some(c) = RegisteredCipher::from_name(name) {
                    let e = c.entry();
                    println!(
                        "  {:<22} {} block={}B rounds={}",
                        name, e.description, e.block_bytes, e.rounds,
                    );
                }
            }
        }
        CryptanalysisOp::Auto {
            cipher,
            boomerang_pairs,
            rectangle_pool,
            trail_threshold_neg_log2,
            trail_top_k,
        } => {
            let opts = AutoAttackOptions {
                boomerang_pairs,
                rectangle_pool,
                trail_threshold: 2f64.powi(-(trail_threshold_neg_log2 as i32)),
                trail_top_k,
            };
            match auto_attack_with(&cipher, opts) {
                Ok(r) => {
                    println!("{}", r.markdown);
                    eprintln!(
                        "[auto-attack] {} sections in {} ms",
                        r.sections_run, r.elapsed_ms,
                    );
                }
                Err(e) => {
                    eprintln!("error: {}", e);
                    eprintln!("Try one of: {}", list_ciphers().join(", "));
                    std::process::exit(1);
                }
            }
        }
        CryptanalysisOp::Boomerang {
            cipher,
            alpha,
            delta,
            pairs,
        } => {
            let c = match RegisteredCipher::from_name(&cipher) {
                Some(c) => c,
                None => {
                    eprintln!("error: unknown cipher: {}", cipher);
                    std::process::exit(1);
                }
            };
            let entry = c.entry();
            let alpha_bytes = match alpha {
                Some(h) => from_hex(&h).expect("bad alpha hex"),
                None => entry.canonical_alpha.clone(),
            };
            let delta_bytes = match delta {
                Some(h) => from_hex(&h).expect("bad delta hex"),
                None => entry.canonical_delta.clone(),
            };
            if alpha_bytes.len() != entry.block_bytes || delta_bytes.len() != entry.block_bytes {
                eprintln!(
                    "error: alpha/delta must be {} bytes for this cipher",
                    entry.block_bytes
                );
                std::process::exit(1);
            }
            let t0 = std::time::Instant::now();
            let r = boomerang_distinguisher(&c, &alpha_bytes, &delta_bytes, pairs, None);
            let elapsed = t0.elapsed().as_millis();
            println!("# Boomerang distinguisher on `{}`", cipher);
            println!();
            println!("- elapsed:       {} ms", elapsed);
            println!("- right quartets: {} / {}", r.right_quartets, r.n_pairs);
            println!("- empirical p:   {:.4e}", r.empirical_probability);
            println!("- random baseline: {:.4e}", r.random_baseline);
            println!(
                "- distinguishes from random (10×): {}",
                if r.distinguishes_from_random(10.0) {
                    "yes"
                } else {
                    "no"
                },
            );
        }
        CryptanalysisOp::Rectangle {
            cipher,
            alpha,
            delta,
            pool,
        } => {
            let c = match RegisteredCipher::from_name(&cipher) {
                Some(c) => c,
                None => {
                    eprintln!("error: unknown cipher: {}", cipher);
                    std::process::exit(1);
                }
            };
            let entry = c.entry();
            let alpha_bytes = match alpha {
                Some(h) => from_hex(&h).expect("bad alpha hex"),
                None => entry.canonical_alpha.clone(),
            };
            let delta_bytes = match delta {
                Some(h) => from_hex(&h).expect("bad delta hex"),
                None => entry.canonical_delta.clone(),
            };
            let t0 = std::time::Instant::now();
            let r = rectangle_attack(&c, &alpha_bytes, &delta_bytes, pool, None);
            let elapsed = t0.elapsed().as_millis();
            println!("# Rectangle attack on `{}`", cipher);
            println!();
            println!("- elapsed:           {} ms", elapsed);
            println!("- pool size:         {}", r.pool_size);
            println!("- right rectangles:  {}", r.right_quartets);
        }
        CryptanalysisOp::Sbox { cipher } => {
            let c = match RegisteredCipher::from_name(&cipher) {
                Some(c) => c,
                None => {
                    eprintln!("error: unknown cipher: {}", cipher);
                    std::process::exit(1);
                }
            };
            match c.sbox() {
                Some(sbox) => {
                    let r = sbox.report();
                    println!("# S-box report for `{}`", cipher);
                    println!();
                    println!("- bits in / out:                 {} / {}", r.n_in, r.n_out);
                    println!("- bijective:                     {}", r.bijective);
                    println!("- balanced:                      {}", r.balanced);
                    println!(
                        "- differential uniformity:       {}",
                        r.differential_uniformity
                    );
                    println!(
                        "- max differential probability:  {:.4}",
                        r.max_differential_probability
                    );
                    println!("- max linear bias:               {:.4}", r.max_linear_bias);
                    println!("- nonlinearity:                  {}", r.nonlinearity);
                    println!("- algebraic degree:              {}", r.algebraic_degree);
                    println!(
                        "- boomerang uniformity:          {}",
                        r.boomerang_uniformity
                            .map(|v| v.to_string())
                            .unwrap_or_else(|| "—".into()),
                    );
                    println!("- max DLCT bias:                 {:.4}", r.max_dlct_bias);
                }
                None => {
                    eprintln!("error: cipher {} has no exposed S-box", cipher);
                    std::process::exit(1);
                }
            }
        }
        CryptanalysisOp::Bench => {
            eprintln!("Running full research bench (this takes ~15-30 s in release mode)...");
            let report = run_full_bench(1);
            println!("{}", report);
        }
        CryptanalysisOp::HashAuto { hash } => {
            use crypto_lib::cryptanalysis::hash_attacks::auto_hash_attack;
            match auto_hash_attack(&hash) {
                Ok(md) => println!("{}", md),
                Err(e) => {
                    eprintln!("error: {}", e);
                    std::process::exit(1);
                }
            }
        }
        CryptanalysisOp::LengthExtension {
            hash,
            secret_len,
            message_hex,
            suffix_hex,
            digest_hex,
        } => {
            use crypto_lib::cryptanalysis::hash_attacks::{
                length_extension_attack, Md4, Md5, MerkleDamgardHash, Sha1,
            };
            let message = from_hex(&message_hex).expect("bad message hex");
            let suffix = from_hex(&suffix_hex).expect("bad suffix hex");
            let digest = from_hex(&digest_hex).expect("bad digest hex");
            let prefix_len = secret_len + message.len();
            fn report_lex<H: MerkleDamgardHash>(
                hash: &H,
                digest: &[u8],
                prefix_len: usize,
                message: &[u8],
                suffix: &[u8],
            ) {
                let (forged, glue) = length_extension_attack(hash, digest, prefix_len, suffix);
                println!("# Length-extension attack on `{}`", hash.name());
                println!();
                println!("Forged digest: {}", to_hex(&forged));
                println!("Glue padding bytes ({} B): {}", glue.len(), to_hex(&glue));
                println!();
                println!("Effective extended message (without secret):");
                println!(
                    "  message ({}B) || glue ({}B) || suffix ({}B)",
                    message.len(),
                    glue.len(),
                    suffix.len()
                );
                let mut combined = Vec::new();
                combined.extend_from_slice(message);
                combined.extend_from_slice(&glue);
                combined.extend_from_slice(suffix);
                println!("Hex: {}", to_hex(&combined));
            }
            match hash.as_str() {
                "md4" => report_lex(&Md4, &digest, prefix_len, &message, &suffix),
                "md5" => report_lex(&Md5, &digest, prefix_len, &message, &suffix),
                "sha1" => report_lex(&Sha1, &digest, prefix_len, &message, &suffix),
                other => {
                    eprintln!("error: unknown hash '{}' (try md4, md5, sha1)", other);
                    std::process::exit(1);
                }
            }
        }
        CryptanalysisOp::AesRelatedKey {
            key_bits,
            delta_k_hex,
            avalanche_trials,
            rounds,
        } => {
            use crypto_lib::cryptanalysis::aes::related_key::{
                biryukov_khovratovich_4round_demo, format_key_schedule_diff,
                key_schedule_difference, related_key_avalanche,
            };
            use crypto_lib::symmetric::aes::AesKey;
            let key_bytes = match key_bits {
                128 => 16,
                256 => 32,
                _ => {
                    eprintln!("error: --key-bits must be 128 or 256");
                    std::process::exit(1);
                }
            };
            let default_rounds = if key_bits == 128 { 10 } else { 14 };
            let rounds = if rounds == 0 { default_rounds } else { rounds };
            let delta_k: Vec<u8> = match &delta_k_hex {
                Some(h) => {
                    let v = from_hex(h).expect("bad --delta-k-hex");
                    if v.len() != key_bytes {
                        eprintln!(
                            "error: --delta-k-hex must be {} bytes (got {})",
                            key_bytes,
                            v.len()
                        );
                        std::process::exit(1);
                    }
                    v
                }
                None => {
                    let mut d = vec![0u8; key_bytes];
                    d[0] = 0x01;
                    d
                }
            };
            // 1. Key-schedule difference propagation.
            let zero_key_bytes = vec![0u8; key_bytes];
            let zero_key = AesKey::new(&zero_key_bytes).unwrap();
            let diff = key_schedule_difference(&zero_key, &delta_k);
            println!(
                "# AES-{} related-key attack against ΔK = 0x{}",
                key_bits,
                to_hex(&delta_k)
            );
            println!();
            println!("## Key-schedule difference propagation\n");
            println!("{}", format_key_schedule_diff(&diff));
            // 2. Related-key avalanche (AES-128 only — the avalanche
            // function uses ReducedAes128 which is 16-byte-key only).
            if key_bits == 128 {
                println!(
                    "## Related-key avalanche ({} rounds, {} trials)\n",
                    rounds, avalanche_trials
                );
                let zero_key_16 = [0u8; 16];
                let mut dk16 = [0u8; 16];
                dk16.copy_from_slice(&delta_k);
                let r = related_key_avalanche(&zero_key_16, &dk16, rounds, avalanche_trials);
                println!("- mean Hamming distance: {:.2} bits", r.mean_distance_bits);
                println!("- ideal-cipher baseline: {:.2} bits", r.ideal_baseline);
                println!("- bias from ideal:       {:.2} bits", r.bias);
                println!();
                println!("## Biryukov-Khovratovich-style local collision\n");
                println!("{} trials at {} rounds:", avalanche_trials, rounds);
                let r = biryukov_khovratovich_4round_demo(avalanche_trials, rounds);
                println!(
                    "- low-Hamming-distance hits: {} / {}",
                    r.collisions, r.n_trials
                );
                println!(
                    "- empirical hit rate:        {:.4}",
                    r.empirical_probability
                );
                println!("- random baseline (Hd≤16 in 128 bits): ≈ 2⁻⁷⁹");
            } else {
                println!("(Related-key avalanche + local-collision demo only available for AES-128; key-schedule diffusion table above already shows the AES-256 weakness.)");
            }
        }
        CryptanalysisOp::AesVisualDemo { demo } => {
            use crypto_lib::cryptanalysis::aes::{
                dfa::{dfa_per_column_candidates, format_dfa_visualization, inject_fault_round_9},
                higher_order::{integral_distinguisher_3_round, render_integral_visualization},
                reduced::ReducedAes128,
                related_key::{format_key_schedule_diff, key_schedule_difference},
                truncated_diff::{
                    propagate_truncated_round, render_trail_diagram, TruncatedPattern,
                },
            };
            use crypto_lib::symmetric::aes::AesKey;
            let run_truncated = |x: &str| -> bool { x == "all" || x == "truncated-diff" };
            let run_integral = |x: &str| -> bool { x == "all" || x == "integral" };
            let run_dfa = |x: &str| -> bool { x == "all" || x == "dfa" };
            let run_ks = |x: &str| -> bool { x == "all" || x == "key-schedule" };
            if run_truncated(&demo) {
                println!("\n## (1) Truncated-differential trail on AES (single active byte)\n");
                let mut current = TruncatedPattern(1);
                let mut steps = Vec::new();
                for _ in 0..4 {
                    let step = propagate_truncated_round(current);
                    let mut next = 0u16;
                    for c in 0..4 {
                        let n = *step.mc_output_counts[c].iter().min().unwrap();
                        for r in 0..n as usize {
                            next |= 1u16 << (4 * c + r);
                        }
                    }
                    current = TruncatedPattern(next);
                    steps.push(step);
                }
                println!("{}", render_trail_diagram(&steps));
            }
            if run_integral(&demo) {
                println!("\n## (2) Square / integral distinguisher on 3-round AES\n");
                let key = [0u8; 16];
                let other = [0u8; 15];
                let xor_sum = integral_distinguisher_3_round(&key, &other);
                println!("{}", render_integral_visualization(&xor_sum));
            }
            if run_dfa(&demo) {
                println!("\n## (3) DFA / Piret-Quisquater single-byte fault on AES-128\n");
                let key_bytes = [
                    0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09,
                    0xcf, 0x4f, 0x3c,
                ];
                let cipher = ReducedAes128::new(&key_bytes, 10, false);
                let pt = [0x6bu8; 16];
                let (correct, faulted) = inject_fault_round_9(&cipher, &pt, 0, 0xCC);
                let cands = dfa_per_column_candidates(&correct, &faulted);
                println!("{}", format_dfa_visualization(&correct, &faulted, &cands));
            }
            if run_ks(&demo) {
                println!("\n## (4) Related-key schedule diffusion (AES-256 vs AES-128)\n");
                for &(bits, label) in &[(128u32, "AES-128"), (256u32, "AES-256")] {
                    let key_bytes = vec![0u8; (bits / 8) as usize];
                    let key = AesKey::new(&key_bytes).unwrap();
                    let mut delta = vec![0u8; (bits / 8) as usize];
                    delta[0] = 0x01;
                    let diff = key_schedule_difference(&key, &delta);
                    println!("\n### {} (ΔK = single byte at position 0)\n", label);
                    println!("{}", format_key_schedule_diff(&diff));
                }
            }
        }
        CryptanalysisOp::VisualAll { target } => {
            use crypto_lib::cryptanalysis::visual_demos::{
                demo_bleichenbacher_bias, demo_birthday_paradox, demo_hnp_recovery_curve,
                demo_j0_twist_factors, demo_joux_multicollision,
                demo_length_extension_diagram, demo_pollard_rho_path, demo_sbox_ddt_heatmap,
                demo_sbox_lat_heatmap, demo_walsh_spectrum, run_all_visual_demos,
            };
            let out: String = match target.as_str() {
                "sbox-ddt" => demo_sbox_ddt_heatmap(),
                "sbox-lat" => demo_sbox_lat_heatmap(),
                "walsh" => demo_walsh_spectrum(),
                "pollard-rho" => demo_pollard_rho_path(),
                "hnp" => demo_hnp_recovery_curve(),
                "bleichenbacher" => demo_bleichenbacher_bias(),
                "length-extension" => demo_length_extension_diagram(),
                "joux" => demo_joux_multicollision(),
                "j0-twists" => demo_j0_twist_factors(),
                "birthday" => demo_birthday_paradox(),
                "all" => run_all_visual_demos(),
                other => {
                    eprintln!("error: unknown target '{}'; try `all` or one of: sbox-ddt, sbox-lat, walsh, pollard-rho, hnp, bleichenbacher, length-extension, joux, j0-twists, birthday", other);
                    std::process::exit(1);
                }
            };
            println!("{}", out);
        }
    }
}

// ── Command implementations ───────────────────────────────────────────────────

fn cmd_hash(alg: &str, msg: &str) {
    match alg {
        "sha256" => println!("{}", to_hex(&sha256(msg.as_bytes()))),
        "sha3" => println!("{}", to_hex(&sha3_256(msg.as_bytes()))),
        "blake3" => println!("{}", to_hex(&blake3_hash(msg.as_bytes()))),
        _ => eprintln!("Unknown algorithm. Use: sha256, sha3, blake3"),
    }
}

fn cmd_ecc_keygen() {
    let curve = CurveParams::secp256k1();
    let kp = EccKeyPair::generate(&curve);

    println!("=== secp256k1 Key Pair ===");
    println!(
        "Private key: {}",
        to_hex(&bigint_to_bytes_be(&kp.private.scalar, 32))
    );
    if let Some(pub_bytes) = kp.public.to_sec1_uncompressed() {
        println!("Public key (uncompressed): {}", to_hex(&pub_bytes));
    }
}

fn cmd_ecdsa(op: EcdsaOp) {
    let curve = CurveParams::secp256k1();
    match op {
        EcdsaOp::Sign { message } => {
            let kp = EccKeyPair::generate(&curve);
            let sig = sign(message.as_bytes(), &kp.private, &curve);
            println!("Message: {message}");
            println!("r: {}", to_hex(&bigint_to_bytes_be(&sig.r, 32)));
            println!("s: {}", to_hex(&bigint_to_bytes_be(&sig.s, 32)));
            if let Some(pub_bytes) = kp.public.to_sec1_uncompressed() {
                println!("Public key: {}", to_hex(&pub_bytes));
            }
            let valid = verify(message.as_bytes(), &kp.public, &sig, &curve);
            println!("Verified: {valid}");
        }
        EcdsaOp::Verify {
            message,
            pubkey_hex,
            r_hex,
            s_hex,
        } => {
            let pub_bytes = from_hex(&pubkey_hex).expect("invalid pubkey hex");
            if pub_bytes.len() != 65 || pub_bytes[0] != 0x04 {
                eprintln!("Expected 65-byte uncompressed public key (04 || X || Y)");
                return;
            }
            let x = BigUint::from_bytes_be(&pub_bytes[1..33]);
            let y = BigUint::from_bytes_be(&pub_bytes[33..65]);
            let p = CurveParams::secp256k1().p;
            let pub_key = EccPublicKey {
                point: Point::Affine {
                    x: FieldElement::new(x, p.clone()),
                    y: FieldElement::new(y, p),
                },
                curve_name: "secp256k1".to_string(),
            };
            let r = BigUint::from_bytes_be(&from_hex(&r_hex).expect("bad r"));
            let s = BigUint::from_bytes_be(&from_hex(&s_hex).expect("bad s"));
            let sig = EcdsaSignature { r, s };
            let valid = verify(message.as_bytes(), &pub_key, &sig, &curve);
            println!("Signature valid: {valid}");
        }
    }
}

fn cmd_ecdh() {
    let curve = CurveParams::secp256k1();
    println!("=== ECDH Key Exchange (secp256k1) ===");

    let alice = EccKeyPair::generate(&curve);
    let bob = EccKeyPair::generate(&curve);

    let sa = ecdh(&alice.private, &bob.public, &curve).unwrap();
    let sb = ecdh(&bob.private, &alice.public, &curve).unwrap();

    println!(
        "Alice private: {}",
        to_hex(&bigint_to_bytes_be(&alice.private.scalar, 32))
    );
    println!(
        "Bob   private: {}",
        to_hex(&bigint_to_bytes_be(&bob.private.scalar, 32))
    );
    println!("Alice shared:  {}", to_hex(&sa));
    println!("Bob   shared:  {}", to_hex(&sb));
    println!("Match: {}", sa == sb);
}

fn cmd_rsa(op: RsaOp) {
    match op {
        RsaOp::Keygen => {
            println!("Generating 1024-bit RSA key pair…");
            let kp = RsaKeyPair::generate(1024);
            println!(
                "n (256-bit prefix): {}…",
                to_hex(&bigint_to_bytes_be(&kp.public.n, 128))[..64].to_string()
            );
            println!("e: {}", kp.public.e);
        }
        RsaOp::Demo => {
            println!("Generating 1024-bit RSA key pair for demo…");
            let kp = RsaKeyPair::generate(1024);
            let msg = b"RSA demo message";
            let sig = rsa_sign(msg, &kp.private);
            let ok = rsa_verify(msg, &sig, &kp.public);
            println!("Message: {:?}", std::str::from_utf8(msg).unwrap());
            println!(
                "Signature (first 16 bytes): {}",
                to_hex(&bigint_to_bytes_be(&sig, 16))
            );
            println!("Verified: {ok}");
        }
    }
}

fn cmd_aes() {
    println!("=== AES-256-GCM ===");
    let key_bytes = [0x42u8; 32];
    let key = AesKey::Aes256(key_bytes);
    let nonce = [0x00u8; 12];
    let pt = b"Hello, AES-GCM!";
    let aad = b"authenticated header";

    let ct = aes_gcm_encrypt(pt, &key, &nonce, aad);
    println!("Plaintext:  {}", std::str::from_utf8(pt).unwrap());
    println!("Ciphertext: {}", to_hex(&ct));

    let recovered = aes_gcm_decrypt(&ct, &key, &nonce, aad).unwrap();
    println!("Decrypted:  {}", std::str::from_utf8(&recovered).unwrap());
}

fn cmd_chacha() {
    println!("=== ChaCha20-Poly1305 ===");
    let key = [0x42u8; 32];
    let nonce = [0x00u8; 12];
    let pt = b"Hello, ChaCha20-Poly1305!";
    let aad = b"additional data";

    let ct = chacha20_poly1305_encrypt(pt, &key, &nonce, aad);
    println!("Plaintext:  {}", std::str::from_utf8(pt).unwrap());
    println!("Ciphertext: {}", to_hex(&ct));

    let recovered = chacha20_poly1305_decrypt(&ct, &key, &nonce, aad).unwrap();
    println!("Decrypted:  {}", std::str::from_utf8(&recovered).unwrap());
}

fn cmd_hkdf() {
    println!("=== HKDF-SHA256 ===");
    let ikm = b"input keying material";
    let salt = b"random salt";
    let info = b"application context";
    let okm = hkdf(Some(salt), ikm, info, 32).expect("demo HKDF parameters must be valid");
    println!("IKM:  {:?}", std::str::from_utf8(ikm).unwrap());
    println!("OKM:  {}", to_hex(&okm));
}

fn cmd_pqc() {
    println!("=== Kyber KEM (simplified ML-KEM educational demo) ===");
    let sk = kyber_keygen();
    println!("Key pair generated.");

    let (ct, ss_enc) = kyber_encapsulate(&sk.public);
    println!("Encapsulated. Shared secret (encaps): {}", to_hex(&ss_enc));

    let ss_dec = kyber_decapsulate(&sk, &ct);
    println!("Decapsulated. Shared secret (decaps): {}", to_hex(&ss_dec));
    println!("Secrets match: {}", ss_enc == ss_dec);
}

fn cmd_demo() {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  Comprehensive Cryptography Library Demo");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // Hashes
    let msg = b"The quick brown fox jumps over the lazy dog";
    println!("[Hash functions]");
    println!("  SHA-256:  {}", to_hex(&sha256(msg)));
    println!("  SHA3-256: {}", to_hex(&sha3_256(msg)));
    println!("  BLAKE3:   {}", to_hex(&blake3_hash(msg)));

    // ECC
    println!("\n[ECC — secp256k1]");
    let curve = CurveParams::secp256k1();
    let kp = EccKeyPair::generate(&curve);
    let sig = sign(msg, &kp.private, &curve);
    println!("  Sign + verify: {}", verify(msg, &kp.public, &sig, &curve));
    let bob = EccKeyPair::generate(&curve);
    let ss = ecdh(&kp.private, &bob.public, &curve).unwrap();
    println!("  ECDH secret:   {}", to_hex(&ss));

    // AES
    println!("\n[AES-256-GCM]");
    let aes_key = AesKey::Aes256([0x13u8; 32]);
    let nonce = [0u8; 12];
    let ct = aes_gcm_encrypt(msg, &aes_key, &nonce, b"");
    let ok = aes_gcm_decrypt(&ct, &aes_key, &nonce, b"").is_ok();
    println!("  Encrypt + decrypt: {ok}");

    // ChaCha20
    println!("\n[ChaCha20-Poly1305]");
    let cha_key = [0x37u8; 32];
    let ct2 = chacha20_poly1305_encrypt(msg, &cha_key, &nonce, b"");
    let ok2 = chacha20_poly1305_decrypt(&ct2, &cha_key, &nonce, b"").is_ok();
    println!("  Encrypt + decrypt: {ok2}");

    // HKDF
    println!("\n[HKDF-SHA256]");
    let okm =
        hkdf(Some(b"salt"), b"secret", b"ctx", 32).expect("demo HKDF parameters must be valid");
    println!("  Derived key: {}", to_hex(&okm));

    // Kyber
    println!("\n[Kyber KEM (simplified)]");
    let sk = kyber_keygen();
    let (ct3, ss3_enc) = kyber_encapsulate(&sk.public);
    let ss3_dec = kyber_decapsulate(&sk, &ct3);
    println!("  KEM secrets match: {}", ss3_enc == ss3_dec);

    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
}

// ── Visual subcommands ────────────────────────────────────────────────────────

fn cmd_visual(op: VisualOp) {
    use crypto_lib::ecc::visualize as ecc_v;
    use crypto_lib::hash::visualize as hash_v;
    use crypto_lib::symmetric::visualize as sym_v;
    match op {
        VisualOp::Symmetric { target } => match target.as_str() {
            "ecb-penguin" => println!("{}", sym_v::demo_ecb_penguin(32, 16)),
            "modes" => println!("{}", sym_v::demo_mode_dataflows()),
            "aes-trace" => {
                println!(
                    "{}",
                    sym_v::demo_aes_state_trace(b"YELLOW SUBMARINE", &[0x42u8; 16])
                )
            }
            "chacha" => println!("{}", sym_v::demo_chacha20_state()),
            "avalanche" => {
                for nr in [1usize, 2, 3, 5, 10] {
                    println!("{}", sym_v::demo_avalanche_matrix(nr));
                }
            }
            "all" => println!("{}", sym_v::run_all_symmetric_demos()),
            other => {
                eprintln!(
                    "error: unknown symmetric target '{}'; try one of: ecb-penguin, modes, aes-trace, chacha, avalanche, all",
                    other
                );
                std::process::exit(1);
            }
        },
        VisualOp::Hash { target } => match target.as_str() {
            "md" => println!("{}", hash_v::demo_merkle_damgard_diagram()),
            "sponge" => println!("{}", hash_v::demo_sponge_diagram()),
            "blake3" => println!("{}", hash_v::demo_blake3_tree_diagram()),
            "sha256-trace" => println!("{}", hash_v::demo_sha256_round_trace(b"abc")),
            "avalanche" => println!("{}", hash_v::demo_hash_avalanche()),
            "all" => println!("{}", hash_v::run_all_hash_demos()),
            other => {
                eprintln!(
                    "error: unknown hash target '{}'; try: md, sponge, blake3, sha256-trace, avalanche, all",
                    other
                );
                std::process::exit(1);
            }
        },
        VisualOp::Ecc { target } => match target.as_str() {
            "curve" => println!("{}", ecc_v::demo_toy_curve_scatter(31, 1, 1)),
            "ops" => println!("{}", ecc_v::demo_point_operations()),
            "ladder" => println!("{}", ecc_v::demo_scalar_mul_ladder(13)),
            "ecdsa" => println!("{}", ecc_v::demo_ecdsa_flow()),
            "all" => println!("{}", ecc_v::run_all_ecc_demos()),
            other => {
                eprintln!(
                    "error: unknown ecc target '{}'; try: curve, ops, ladder, ecdsa, all",
                    other
                );
                std::process::exit(1);
            }
        },
        VisualOp::All => {
            println!("{}", sym_v::run_all_symmetric_demos());
            println!("\n\n========================================\n\n");
            println!("{}", hash_v::run_all_hash_demos());
            println!("\n\n========================================\n\n");
            println!("{}", ecc_v::run_all_ecc_demos());
        }
    }
}
