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
    asymmetric::rsa::{RsaKeyPair, rsa_sign, rsa_verify},
    ecc::{
        curve::CurveParams,
        ecdh::ecdh,
        ecdsa::{sign, verify, EcdsaSignature},
        keys::{EccKeyPair, EccPublicKey},
        point::Point,
        field::FieldElement,
    },
    hash::{sha256, sha3_256, blake3_hash},
    kdf::hkdf::hkdf,
    pqc::kyber::{kyber_keygen, kyber_encapsulate, kyber_decapsulate},
    symmetric::{
        aes::{AesKey, aes_gcm_encrypt, aes_gcm_decrypt},
        chacha20::{chacha20_poly1305_encrypt, chacha20_poly1305_decrypt},
    },
    utils::encoding::{to_hex, from_hex, bigint_to_bytes_be},
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
    Hash {
        algorithm: String,
        message: String,
    },
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
}

#[derive(Subcommand)]
enum EccOp {
    Keygen,
}

#[derive(Subcommand)]
enum EcdsaOp {
    Sign   { message: String },
    Verify { message: String, pubkey_hex: String, r_hex: String, s_hex: String },
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
        Cmd::Ecc { op: EccOp::Keygen }  => cmd_ecc_keygen(),
        Cmd::Ecdsa { op }               => cmd_ecdsa(op),
        Cmd::Ecdh                       => cmd_ecdh(),
        Cmd::Rsa { op }                 => cmd_rsa(op),
        Cmd::Aes                        => cmd_aes(),
        Cmd::Chacha                     => cmd_chacha(),
        Cmd::Hkdf                       => cmd_hkdf(),
        Cmd::Pqc                        => cmd_pqc(),
        Cmd::Demo                       => cmd_demo(),
        Cmd::Cryptanalysis { op }       => cmd_cryptanalysis(op),
    }
}

// ── Cryptanalysis subcommands ─────────────────────────────────────────────────

fn cmd_cryptanalysis(op: CryptanalysisOp) {
    use crypto_lib::cryptanalysis::auto_attack::{auto_attack_with, AutoAttackOptions};
    use crypto_lib::cryptanalysis::boomerang::{
        boomerang_distinguisher, rectangle_attack,
    };
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
            if alpha_bytes.len() != entry.block_bytes
                || delta_bytes.len() != entry.block_bytes
            {
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
                    println!("- differential uniformity:       {}", r.differential_uniformity);
                    println!("- max differential probability:  {:.4}", r.max_differential_probability);
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
    }
}

// ── Command implementations ───────────────────────────────────────────────────

fn cmd_hash(alg: &str, msg: &str) {
    match alg {
        "sha256" => println!("{}", to_hex(&sha256(msg.as_bytes()))),
        "sha3"   => println!("{}", to_hex(&sha3_256(msg.as_bytes()))),
        "blake3" => println!("{}", to_hex(&blake3_hash(msg.as_bytes()))),
        _        => eprintln!("Unknown algorithm. Use: sha256, sha3, blake3"),
    }
}

fn cmd_ecc_keygen() {
    let curve = CurveParams::secp256k1();
    let kp = EccKeyPair::generate(&curve);

    println!("=== secp256k1 Key Pair ===");
    println!("Private key: {}", to_hex(&bigint_to_bytes_be(&kp.private.scalar, 32)));
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
        EcdsaOp::Verify { message, pubkey_hex, r_hex, s_hex } => {
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
    let bob   = EccKeyPair::generate(&curve);

    let sa = ecdh(&alice.private, &bob.public, &curve).unwrap();
    let sb = ecdh(&bob.private, &alice.public, &curve).unwrap();

    println!("Alice private: {}", to_hex(&bigint_to_bytes_be(&alice.private.scalar, 32)));
    println!("Bob   private: {}", to_hex(&bigint_to_bytes_be(&bob.private.scalar, 32)));
    println!("Alice shared:  {}", to_hex(&sa));
    println!("Bob   shared:  {}", to_hex(&sb));
    println!("Match: {}", sa == sb);
}

fn cmd_rsa(op: RsaOp) {
    match op {
        RsaOp::Keygen => {
            println!("Generating 1024-bit RSA key pair…");
            let kp = RsaKeyPair::generate(1024);
            println!("n (256-bit prefix): {}…", to_hex(&bigint_to_bytes_be(&kp.public.n, 128))[..64].to_string());
            println!("e: {}", kp.public.e);
        }
        RsaOp::Demo => {
            println!("Generating 1024-bit RSA key pair for demo…");
            let kp = RsaKeyPair::generate(1024);
            let msg = b"RSA demo message";
            let sig = rsa_sign(msg, &kp.private);
            let ok  = rsa_verify(msg, &sig, &kp.public);
            println!("Message: {:?}", std::str::from_utf8(msg).unwrap());
            println!("Signature (first 16 bytes): {}", to_hex(&bigint_to_bytes_be(&sig, 16)));
            println!("Verified: {ok}");
        }
    }
}

fn cmd_aes() {
    println!("=== AES-256-GCM ===");
    let key_bytes = [0x42u8; 32];
    let key   = AesKey::Aes256(key_bytes);
    let nonce = [0x00u8; 12];
    let pt    = b"Hello, AES-GCM!";
    let aad   = b"authenticated header";

    let ct = aes_gcm_encrypt(pt, &key, &nonce, aad);
    println!("Plaintext:  {}", std::str::from_utf8(pt).unwrap());
    println!("Ciphertext: {}", to_hex(&ct));

    let recovered = aes_gcm_decrypt(&ct, &key, &nonce, aad).unwrap();
    println!("Decrypted:  {}", std::str::from_utf8(&recovered).unwrap());
}

fn cmd_chacha() {
    println!("=== ChaCha20-Poly1305 ===");
    let key   = [0x42u8; 32];
    let nonce = [0x00u8; 12];
    let pt    = b"Hello, ChaCha20-Poly1305!";
    let aad   = b"additional data";

    let ct = chacha20_poly1305_encrypt(pt, &key, &nonce, aad);
    println!("Plaintext:  {}", std::str::from_utf8(pt).unwrap());
    println!("Ciphertext: {}", to_hex(&ct));

    let recovered = chacha20_poly1305_decrypt(&ct, &key, &nonce, aad).unwrap();
    println!("Decrypted:  {}", std::str::from_utf8(&recovered).unwrap());
}

fn cmd_hkdf() {
    println!("=== HKDF-SHA256 ===");
    let ikm  = b"input keying material";
    let salt = b"random salt";
    let info = b"application context";
    let okm  = hkdf(Some(salt), ikm, info, 32);
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
    let nonce   = [0u8; 12];
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
    let okm = hkdf(Some(b"salt"), b"secret", b"ctx", 32);
    println!("  Derived key: {}", to_hex(&okm));

    // Kyber
    println!("\n[Kyber KEM (simplified)]");
    let sk = kyber_keygen();
    let (ct3, ss3_enc) = kyber_encapsulate(&sk.public);
    let ss3_dec = kyber_decapsulate(&sk, &ct3);
    println!("  KEM secrets match: {}", ss3_enc == ss3_dec);

    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
}
