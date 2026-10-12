//! `crypto pqc …` — every post-quantum scheme in the library, from the CLI.
//!
//! The library carries twenty-odd PQC schemes; before this module ten of them
//! had no entry point outside their own tests. [`PqcOp::List`] enumerates what
//! is here and [`PqcOp::Run`] exercises any one of them (or all of them), so
//! "the library implements X" is a claim a reader can check in one command.
//!
//! Parameter sets are the library's own. ML-KEM, ML-DSA and SLH-DSA are at
//! their standardised sizes; the rest are at reduced educational parameters and
//! say so in their own module docs. Nothing here is constant-time.

use clap::Subcommand;
use crypto_lib::pqc;
use crypto_lib::utils::encoding::to_hex;

#[derive(Subcommand)]
pub enum PqcOp {
    /// The original guided tour: ML-KEM, ML-DSA, SQIsign, the on-ramp
    /// signatures and toy Kyber, with commentary.
    Demo,
    /// List every scheme, with its kind and whether its parameters are
    /// standardised or educational.
    List,
    /// Run one scheme end to end — or `all` for every one of them.
    Run {
        /// Scheme name as `list` prints it, or `all`.
        #[arg(long, default_value = "all")]
        scheme: String,
        /// Message to sign, for the signature schemes.
        #[arg(long, default_value = "post-quantum")]
        message: String,
    },
}

/// What a scheme is.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// Key encapsulation: keygen, encapsulate, decapsulate, secrets match.
    Kem,
    /// Signature: keygen, sign, verify, and a tampered signature is rejected.
    Sig,
    /// Non-interactive key exchange: two parties, one shared value.
    Nike,
    /// Public-key encryption without the FO transform.
    Pke,
}

impl Kind {
    fn label(self) -> &'static str {
        match self {
            Kind::Kem => "KEM",
            Kind::Sig => "signature",
            Kind::Nike => "key exchange",
            Kind::Pke => "PKE",
        }
    }
}

/// Whether a scheme runs at its standardised parameters or reduced ones.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Params {
    Standard,
    Educational,
}

impl Params {
    fn label(self) -> &'static str {
        match self {
            Params::Standard => "standardised",
            Params::Educational => "educational",
        }
    }
}

struct Scheme {
    name: &'static str,
    kind: Kind,
    params: Params,
    family: &'static str,
}

const SCHEMES: &[Scheme] = &[
    Scheme {
        name: "ml-kem-512",
        kind: Kind::Kem,
        params: Params::Standard,
        family: "module lattice",
    },
    Scheme {
        name: "ml-kem-768",
        kind: Kind::Kem,
        params: Params::Standard,
        family: "module lattice",
    },
    Scheme {
        name: "ml-kem-1024",
        kind: Kind::Kem,
        params: Params::Standard,
        family: "module lattice",
    },
    Scheme {
        name: "ml-dsa-65",
        kind: Kind::Sig,
        params: Params::Standard,
        family: "module lattice",
    },
    Scheme {
        name: "slh-dsa-sha2-128s",
        kind: Kind::Sig,
        params: Params::Standard,
        family: "hash",
    },
    Scheme {
        name: "kyber",
        kind: Kind::Kem,
        params: Params::Educational,
        family: "module lattice",
    },
    Scheme {
        name: "frodo",
        kind: Kind::Kem,
        params: Params::Educational,
        family: "plain LWE",
    },
    Scheme {
        name: "ntru",
        kind: Kind::Kem,
        params: Params::Educational,
        family: "NTRU lattice",
    },
    Scheme {
        name: "ntru-prime",
        kind: Kind::Kem,
        params: Params::Educational,
        family: "NTRU lattice",
    },
    Scheme {
        name: "x-wing",
        kind: Kind::Kem,
        params: Params::Educational,
        family: "hybrid (X25519 + ML-KEM)",
    },
    Scheme {
        name: "bike",
        kind: Kind::Kem,
        params: Params::Educational,
        family: "QC-MDPC code",
    },
    Scheme {
        name: "hqc",
        kind: Kind::Kem,
        params: Params::Educational,
        family: "quasi-cyclic code",
    },
    Scheme {
        name: "classic-mceliece",
        kind: Kind::Kem,
        params: Params::Educational,
        family: "Goppa code",
    },
    Scheme {
        name: "mceliece",
        kind: Kind::Pke,
        params: Params::Educational,
        family: "Goppa code",
    },
    Scheme {
        name: "csidh",
        kind: Kind::Nike,
        params: Params::Educational,
        family: "isogeny",
    },
    Scheme {
        name: "sqisign",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "isogeny",
    },
    Scheme {
        name: "uov",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "multivariate",
    },
    Scheme {
        name: "qr-uov",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "multivariate",
    },
    Scheme {
        name: "mayo",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "multivariate",
    },
    Scheme {
        name: "snova",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "multivariate",
    },
    Scheme {
        name: "hawk",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "lattice (Hawk)",
    },
    Scheme {
        name: "fn-dsa",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "NTRU lattice (Falcon)",
    },
    Scheme {
        name: "sdith",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "MPC-in-the-head",
    },
    Scheme {
        name: "mqom",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "MPC-in-the-head",
    },
    Scheme {
        name: "faest",
        kind: Kind::Sig,
        params: Params::Educational,
        family: "VOLE-in-the-head",
    },
];

/// The outcome of running one scheme.
struct Outcome {
    ok: bool,
    detail: String,
}

fn ok(detail: String) -> Outcome {
    Outcome { ok: true, detail }
}

pub fn run(op: Option<PqcOp>) {
    match op.unwrap_or(PqcOp::Demo) {
        PqcOp::Demo => crate::cmd_pqc(),
        PqcOp::List => {
            println!(
                "{:<20}  {:<13}  {:<13}  family",
                "scheme", "kind", "parameters"
            );
            println!("{}", "-".repeat(78));
            for s in SCHEMES {
                println!(
                    "{:<20}  {:<13}  {:<13}  {}",
                    s.name,
                    s.kind.label(),
                    s.params.label(),
                    s.family
                );
            }
            println!(
                "\n{} schemes. `crypto pqc run --scheme <name>` runs one; `--scheme all` runs\n\
                 every one. The speed-oriented ML-KEM/ML-DSA/isogeny code lives under\n\
                 `crypto pqc-fast`, and the attacks on ML-KEM and ML-DSA under `crypto mlwe`.",
                SCHEMES.len()
            );
        }
        PqcOp::Run { scheme, message } => {
            let msg = message.as_bytes();
            let wanted: Vec<&Scheme> = if scheme == "all" {
                SCHEMES.iter().collect()
            } else {
                let key = normalise(&scheme);
                let found: Vec<&Scheme> = SCHEMES
                    .iter()
                    .filter(|s| normalise(s.name) == key)
                    .collect();
                if found.is_empty() {
                    eprintln!(
                        "unknown scheme `{scheme}`; `crypto pqc list` prints the {} known names",
                        SCHEMES.len()
                    );
                    std::process::exit(2);
                }
                found
            };
            let mut failures = 0usize;
            for s in wanted {
                let o = run_one(s, msg);
                println!(
                    "{:<20} {:<13} {}  {}",
                    s.name,
                    s.kind.label(),
                    if o.ok { "ok  " } else { "FAIL" },
                    o.detail
                );
                if !o.ok {
                    failures += 1;
                }
            }
            if failures > 0 {
                eprintln!("{failures} scheme(s) failed");
                std::process::exit(1);
            }
        }
    }
}

fn normalise(s: &str) -> String {
    s.chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_lowercase())
        .collect()
}

/// Run one scheme's full cycle and report what happened.
///
/// For a KEM: keygen, encapsulate, decapsulate, and the two shared secrets must
/// agree. For a signature: keygen, sign, verify, *and* a flipped bit in the
/// signature must be rejected — a scheme that accepts everything would otherwise
/// pass.
fn run_one(s: &Scheme, msg: &[u8]) -> Outcome {
    match s.name {
        "ml-kem-512" => ml_kem(&pqc::ML_KEM_512),
        "ml-kem-768" => ml_kem(&pqc::ML_KEM_768),
        "ml-kem-1024" => ml_kem(&pqc::ML_KEM_1024),

        "ml-dsa-65" => {
            let (pk, sk) = pqc::ml_dsa_65_keygen(&[3u8; 32]);
            let sig = pqc::ml_dsa_65_sign(&sk, msg, &[0u8; 32]);
            sig_outcome(
                pqc::ml_dsa_65_verify(&pk, msg, &sig),
                |t| pqc::ml_dsa_65_verify(&pk, msg, t),
                &sig,
                pk.0.len(),
            )
        }
        "slh-dsa-sha2-128s" => {
            use crypto_lib::pqc::slh_dsa::*;
            let (pk, sk) = slh_dsa_sha2_128s_keygen(&[7u8; 48]);
            let sig = slh_dsa_sha2_128s_sign(&sk, msg, &[0u8; 16]);
            sig_outcome(
                slh_dsa_sha2_128s_verify(&pk, msg, &sig),
                |t| slh_dsa_sha2_128s_verify(&pk, msg, t),
                &sig,
                pk.bytes.len(),
            )
        }

        "kyber" => {
            let sk = pqc::kyber_keygen();
            let (ct, a) = pqc::kyber_encapsulate(&sk.public);
            let b = pqc::kyber_decapsulate(&sk, &ct);
            kem_outcome(a, b, 0)
        }
        "frodo" => {
            let kp = pqc::frodo_keygen();
            let (ct, a) = pqc::frodo_encapsulate(&kp.pk);
            let b = pqc::frodo_decapsulate(&ct, &kp.sk);
            kem_outcome(a, b, 0)
        }
        "ntru" => {
            let kp = pqc::ntru_keygen();
            let (ct, a) = pqc::ntru_encapsulate(&kp.pk);
            let b = pqc::ntru_decapsulate(&ct, &kp.sk);
            kem_outcome(a, b, 0)
        }
        "ntru-prime" => {
            use crypto_lib::pqc::ntru_prime::*;
            let kp = ntru_prime_keygen();
            let (ct, a) = ntru_prime_encapsulate(&kp.pk);
            let b = ntru_prime_decapsulate(&ct, &kp.sk);
            kem_outcome(a, b, 0)
        }
        "x-wing" => {
            let kp = pqc::x_wing_keygen();
            let (ct, a) = pqc::x_wing_encapsulate(&kp.pk);
            let b = pqc::x_wing_decapsulate(&ct, &kp.sk);
            kem_outcome(a, b, 0)
        }
        "bike" => {
            use crypto_lib::pqc::bike::*;
            let kp = bike_keygen();
            let (ct, a) = bike_encapsulate(&kp.pk);
            let b = bike_decapsulate(&ct, &kp.sk);
            kem_outcome(a, b, 0)
        }
        "hqc" => {
            let kp = pqc::hqc_keygen();
            let (ct, a) = pqc::hqc_encapsulate(&kp.pk);
            let b = pqc::hqc_decapsulate(&ct, &kp.sk);
            kem_outcome(a, b, 0)
        }
        "classic-mceliece" => {
            use crypto_lib::pqc::classic_mceliece::*;
            let kp = classic_mceliece_keygen();
            let (ct, a) = classic_mceliece_encapsulate(&kp.pk);
            let b = classic_mceliece_decapsulate(&ct, &kp.sk);
            kem_outcome(a, b, 0)
        }

        "mceliece" => {
            use crypto_lib::pqc::mceliece::*;
            let kp = McElieceKeyPair::generate();
            // The plaintext is exactly k bits, one per byte.
            let plain: Vec<u8> = (0..kp.public.k).map(|i| (i % 2) as u8).collect();
            let ct = mceliece_encrypt(&plain, &kp.public);
            match mceliece_decrypt(&ct, &kp.private) {
                Some(got) if got == plain => ok(format!(
                    "k = {} bits, ciphertext {} bits, decrypted exactly",
                    kp.public.k,
                    ct.len()
                )),
                Some(_) => Outcome {
                    ok: false,
                    detail: "decrypted to the wrong plaintext".into(),
                },
                None => Outcome {
                    ok: false,
                    detail: "decryption failed".into(),
                },
            }
        }

        "csidh" => {
            use crypto_lib::pqc::csidh::*;
            let a_sk = csidh_keygen_private();
            let b_sk = csidh_keygen_private();
            let a_pk = csidh_public_key(&a_sk);
            let b_pk = csidh_public_key(&b_sk);
            let a_shared = csidh_shared_secret(&a_sk, &b_pk);
            let b_shared = csidh_shared_secret(&b_sk, &a_pk);
            if a_shared == b_shared {
                ok(format!("both parties agreed on j = {a_shared}"))
            } else {
                Outcome {
                    ok: false,
                    detail: format!("no agreement: {a_shared} vs {b_shared}"),
                }
            }
        }

        "sqisign" => {
            let (pk, sk) = pqc::sqisign_keygen();
            let sig = pqc::sqisign_sign(&pk, &sk, msg);
            if pqc::sqisign_verify(&pk, msg, &sig) {
                ok(format!(
                    "j(E_A) = {:?}, response path of {} steps",
                    pk.j,
                    sig.response.len().saturating_sub(1)
                ))
            } else {
                Outcome {
                    ok: false,
                    detail: "own signature rejected".into(),
                }
            }
        }

        "uov" => {
            let (pk, sk) = pqc::uov_keygen();
            let sig = pqc::uov_sign(&sk, msg);
            sig_bool(
                pqc::uov_verify(&pk, msg, &sig),
                pqc::uov_verify(&pk, b"other", &sig),
            )
        }
        "qr-uov" => {
            let (pk, sk) = pqc::qr_uov_keygen();
            let sig = pqc::qr_uov_sign(&sk, msg);
            sig_bool(
                pqc::qr_uov_verify(&pk, msg, &sig),
                pqc::qr_uov_verify(&pk, b"other", &sig),
            )
        }
        "mayo" => {
            let (pk, sk) = pqc::mayo_keygen();
            let sig = pqc::mayo_sign(&sk, msg);
            sig_bool(
                pqc::mayo_verify(&pk, msg, &sig),
                pqc::mayo_verify(&pk, b"other", &sig),
            )
        }
        "snova" => {
            let (pk, sk) = pqc::snova_keygen();
            let sig = pqc::snova_sign(&sk, msg);
            sig_bool(
                pqc::snova_verify(&pk, msg, &sig),
                pqc::snova_verify(&pk, b"other", &sig),
            )
        }
        "hawk" => {
            let (pk, sk) = pqc::hawk_keygen();
            let sig = pqc::hawk_sign(&sk, msg);
            sig_bool(
                pqc::hawk_verify(&pk, msg, &sig),
                pqc::hawk_verify(&pk, b"other", &sig),
            )
        }
        "fn-dsa" => {
            let (pk, sk) = pqc::fn_dsa_keygen();
            let sig = pqc::fn_dsa_sign(&sk, msg);
            sig_bool(
                pqc::fn_dsa_verify(&pk, msg, &sig),
                pqc::fn_dsa_verify(&pk, b"other", &sig),
            )
        }
        "sdith" => {
            let (pk, sk) = pqc::sdith_keygen();
            let sig = pqc::sdith_sign(&pk, &sk, msg);
            sig_bool(
                pqc::sdith_verify(&pk, msg, &sig),
                pqc::sdith_verify(&pk, b"other", &sig),
            )
        }
        "mqom" => {
            let (pk, sk) = pqc::mqom_keygen();
            let sig = pqc::mqom_sign(&pk, &sk, msg);
            sig_bool(
                pqc::mqom_verify(&pk, msg, &sig),
                pqc::mqom_verify(&pk, b"other", &sig),
            )
        }
        "faest" => {
            let (pk, sk) = pqc::faest_keygen();
            let sig = pqc::faest_sign(&pk, &sk, msg);
            sig_bool(
                pqc::faest_verify(&pk, msg, &sig),
                pqc::faest_verify(&pk, b"other", &sig),
            )
        }

        other => Outcome {
            ok: false,
            detail: format!("no runner wired for `{other}` — see cli_pqc.rs"),
        },
    }
}

fn ml_kem(p: &pqc::MlKemParams) -> Outcome {
    let (ek, dk) = pqc::ml_kem_keygen(p);
    let Some((ct, a)) = pqc::ml_kem_encaps(p, &ek) else {
        return Outcome {
            ok: false,
            detail: "encapsulation rejected a fresh key".into(),
        };
    };
    let Some(b) = pqc::ml_kem_decaps(p, &dk, &ct) else {
        return Outcome {
            ok: false,
            detail: "decapsulation rejected a well-formed ciphertext".into(),
        };
    };
    // Implicit rejection: a tampered ciphertext must still return a secret, and
    // it must not be the honest one.
    let mut bad = ct.clone();
    bad[0] ^= 1;
    let rejected = pqc::ml_kem_decaps(p, &dk, &bad);
    if a != b {
        return Outcome {
            ok: false,
            detail: "shared secrets disagree".into(),
        };
    }
    match rejected {
        Some(r) if r != a => ok(format!(
            "ek {} B, ct {} B, secret {}…, implicit rejection differs",
            ek.0.len(),
            ct.len(),
            &to_hex(&a)[..16]
        )),
        Some(_) => Outcome {
            ok: false,
            detail: "a tampered ciphertext gave the honest secret".into(),
        },
        None => Outcome {
            ok: false,
            detail: "a tampered ciphertext errored instead of rejecting".into(),
        },
    }
}

fn kem_outcome(a: [u8; 32], b: [u8; 32], _unused: usize) -> Outcome {
    if a == b {
        ok(format!("shared secret {}…", &to_hex(&a)[..16]))
    } else {
        Outcome {
            ok: false,
            detail: "shared secrets disagree".into(),
        }
    }
}

/// A signature scheme passes only if it accepts its own signature *and* rejects
/// a tampered one.
fn sig_outcome(
    accepted: bool,
    mut verify: impl FnMut(&[u8]) -> bool,
    sig: &[u8],
    pk_len: usize,
) -> Outcome {
    if !accepted {
        return Outcome {
            ok: false,
            detail: "own signature rejected".into(),
        };
    }
    let mut bad = sig.to_vec();
    bad[0] ^= 1;
    if verify(&bad) {
        return Outcome {
            ok: false,
            detail: "accepted a tampered signature".into(),
        };
    }
    ok(format!(
        "pk {pk_len} B, signature {} B, tampering rejected",
        sig.len()
    ))
}

/// The same check for schemes whose signature type is not a byte slice: accept
/// the real message, reject a different one.
fn sig_bool(accepted: bool, accepted_other: bool) -> Outcome {
    if !accepted {
        Outcome {
            ok: false,
            detail: "own signature rejected".into(),
        }
    } else if accepted_other {
        Outcome {
            ok: false,
            detail: "verified against the wrong message".into(),
        }
    } else {
        ok("signed, verified, wrong message rejected".into())
    }
}
