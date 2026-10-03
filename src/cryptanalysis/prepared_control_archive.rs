//! Strict bounded USTAR inventory verification for retained control capsules.
//! Parses data only; it never extracts files or executes archived programs.
use super::prepared_sat_control::sha256;
use serde_json::{json, Value};
use std::collections::BTreeSet;
fn require(ok: bool, message: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(message.into())
    }
}
fn desc(bytes: &[u8]) -> Value {
    json!({"bytes":bytes.len(),"sha256":sha256(bytes)})
}
fn text(bytes: &[u8]) -> Result<String, String> {
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    require(
        bytes[end..].iter().all(|&b| b == 0 || b == b' '),
        "invalid tar header text",
    )?;
    Ok(std::str::from_utf8(&bytes[..end])
        .map_err(|e| e.to_string())?
        .trim_end()
        .into())
}
fn octal(bytes: &[u8]) -> Result<usize, String> {
    let value = std::str::from_utf8(bytes)
        .map_err(|e| e.to_string())?
        .trim_matches(['\0', ' ']);
    require(
        !value.is_empty() && value.bytes().all(|b| (b'0'..=b'7').contains(&b)),
        "invalid tar octal field",
    )?;
    usize::from_str_radix(value, 8).map_err(|e| e.to_string())
}
pub fn safe_path(name: &str) -> bool {
    !name.is_empty()
        && !name.starts_with('/')
        && !name.contains('\\')
        && name
            .trim_end_matches('/')
            .split('/')
            .all(|c| !c.is_empty() && c != "." && c != "..")
}
pub fn verify_tar_inventory(tar: &[u8], expected: &Value) -> Result<usize, String> {
    let expected = expected.as_object().ok_or("missing archive inventory")?;
    let mut seen = BTreeSet::new();
    let mut members = BTreeSet::new();
    let mut offset = 0;
    while offset + 512 <= tar.len() {
        let h = &tar[offset..offset + 512];
        if h.iter().all(|&b| b == 0) {
            require(
                tar.len() - offset >= 1024
                    && (tar.len() - offset).is_multiple_of(512)
                    && tar[offset..].iter().all(|&b| b == 0),
                "invalid tar terminal blocks",
            )?;
            require(
                seen.len() == expected.len(),
                "archive omits registered files",
            )?;
            return Ok(seen.len());
        }
        let sum: usize = h
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                if (148..156).contains(&i) {
                    32
                } else {
                    b as usize
                }
            })
            .sum();
        require(sum == octal(&h[148..156])?, "tar checksum differs")?;
        require(&h[257..265] == b"ustar\x0000", "archive is not plain USTAR")?;
        let name = text(&h[..100])?;
        let prefix = text(&h[345..500])?;
        let name = if prefix.is_empty() {
            name
        } else {
            format!("{prefix}/{name}")
        };
        require(safe_path(&name), "unsafe archive path")?;
        require(members.insert(name.clone()), "duplicate archive member")?;
        let size = octal(&h[124..136])?;
        offset += 512;
        let end = offset.checked_add(size).ok_or("archive size overflow")?;
        require(end <= tar.len(), "truncated archive member")?;
        match h[156] {
            b'5' => {
                require(
                    size == 0 && name.ends_with('/'),
                    "invalid archive directory",
                )?;
                require(
                    expected.keys().any(|k| k.starts_with(&name)),
                    "unregistered archive directory",
                )?;
            }
            0 | b'0' => {
                require(seen.insert(name.clone()), "duplicate archive member")?;
                require(
                    expected.get(&name) == Some(&desc(&tar[offset..end])),
                    "archive file differs or is unregistered",
                )?;
            }
            _ => return Err("archive contains a link or special member".into()),
        }
        offset = end
            .checked_add((512 - size % 512) % 512)
            .ok_or("archive padding overflow")?;
        require(offset <= tar.len(), "truncated archive padding")?;
    }
    Err("archive has no terminal blocks".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn tar_file(name: &str, bytes: &[u8], kind: u8) -> Vec<u8> {
        let mut h = vec![0; 512];
        h[..name.len()].copy_from_slice(name.as_bytes());
        h[124..136].copy_from_slice(format!("{:011o}\0", bytes.len()).as_bytes());
        h[156] = kind;
        h[257..265].copy_from_slice(b"ustar\x0000");
        let sum: usize = h
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                if (148..156).contains(&i) {
                    32
                } else {
                    b as usize
                }
            })
            .sum();
        h[148..156].copy_from_slice(format!("{sum:06o}\0 ").as_bytes());
        h.extend(bytes);
        h.resize(h.len().div_ceil(512) * 512, 0);
        h.extend(vec![0; 1024]);
        h
    }
    #[test]
    fn mutation_omission_duplicate_links_and_unsafe_paths_reject() {
        let expected = json!({"immutable/a":desc(b"frozen")});
        let good = tar_file("immutable/a", b"frozen", b'0');
        assert_eq!(verify_tar_inventory(&good, &expected).unwrap(), 1);
        let mut bad = good.clone();
        bad[512] ^= 1;
        assert!(verify_tar_inventory(&bad, &expected).is_err());
        assert!(verify_tar_inventory(&[0; 1024], &expected).is_err());
        let mut duplicate = good[..1024].to_vec();
        duplicate.extend(&good);
        assert!(verify_tar_inventory(&duplicate, &expected).is_err());
        for name in ["../escape", "/absolute", "immutable//a", "immutable/./a"] {
            assert!(verify_tar_inventory(&tar_file(name, b"frozen", b'0'), &expected).is_err());
        }
        assert!(
            verify_tar_inventory(&tar_file("immutable/a", b"frozen", b'2'), &expected).is_err()
        );
        assert!(verify_tar_inventory(&good[..good.len() - 1], &expected).is_err());
        let mut trailer = good;
        trailer.push(1);
        assert!(verify_tar_inventory(&trailer, &expected).is_err());
    }
    #[test]
    fn directories_must_be_registered_and_unique() {
        let expected = json!({"immutable/a":desc(b"frozen")});
        let dir = tar_file("immutable/", b"", b'5');
        let file = tar_file("immutable/a", b"frozen", b'0');
        let mut valid = dir[..512].to_vec();
        valid.extend(&file);
        assert_eq!(verify_tar_inventory(&valid, &expected).unwrap(), 1);
        let mut duplicate = dir[..512].to_vec();
        duplicate.extend(&valid);
        assert!(verify_tar_inventory(&duplicate, &expected).is_err());
        let mut unregistered = tar_file("unregistered/", b"", b'5')[..512].to_vec();
        unregistered.extend(&file);
        assert!(verify_tar_inventory(&unregistered, &expected).is_err());
    }
}
