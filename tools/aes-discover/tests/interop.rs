use aes_discover::cipher::{Aes, Config};
use aes_discover::{hex, stream};
use rand::RngCore;
use std::io::Write;
use std::process::{Command, Stdio};

/// Independent established implementation; fail visibly if OpenSSL is absent.
#[test]
fn openssl_differential_reference_all_key_sizes() {
    for bits in [128, 192, 256] {
        let mut rng = stream(20261017, "openssl/fixtures", bits as u64);
        for _ in 0..8 {
            let mut key = vec![0; bits / 8];
            rng.fill_bytes(&mut key);
            let mut pt = [0; 16];
            rng.fill_bytes(&mut pt);
            let aes = Aes::new(
                &key,
                Config {
                    key_bits: bits,
                    rounds: bits / 32 + 6,
                    final_mix_columns: false,
                },
            )
            .unwrap();
            let mut child = Command::new("openssl")
                .args([
                    "enc",
                    &format!("-aes-{bits}-ecb"),
                    "-K",
                    &hex(&key),
                    "-nopad",
                    "-nosalt",
                ])
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .expect("OpenSSL is a required test reference");
            child.stdin.take().unwrap().write_all(&pt).unwrap();
            let out = child.wait_with_output().unwrap();
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
            assert_eq!(out.stdout, aes.encrypt(pt));
        }
    }
}

#[test]
fn cli_rejects_bad_keys_and_round_counts() {
    let exe = env!("CARGO_BIN_EXE_aes-discover");
    let o = Command::new(exe)
        .args([
            "encrypt",
            "--key-bits",
            "192",
            "--rounds",
            "13",
            "--key",
            &"00".repeat(24),
            "--plaintext",
            &"00".repeat(16),
        ])
        .output()
        .unwrap();
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("rounds must be"));
    let o = Command::new(exe)
        .args(["encrypt", "--key", "00", "--plaintext", &"00".repeat(16)])
        .output()
        .unwrap();
    assert!(!o.status.success());
}
