use std::collections::BTreeSet;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use super::arithmetic::{
    centered_u, d23_control, norm, p192_parameters, p192_sqrt_discriminant_control, primitive,
};
use super::certificate::{self, Disposition, RationalFactor, SmallEntry};
use super::digest::sha256;
use super::factor_base::{mod_u64, mul_mod, Entry};
use super::Result;

pub const REFERENCE_SHELL: &str = "REF-0";
pub const REFERENCE_V_MAX: u64 = 2;
pub const REFERENCE_X_MAX_INCLUSIVE: u64 = 1024;
pub const RECORD_COUNT: u64 = 1025;
pub const RECORD_BYTES: u64 = RECORD_COUNT * 16;
pub const SHARD_COUNT: u64 = 1;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReferenceRecord {
    pub index: u64,
    pub v: u64,
    pub x: u64,
}

impl ReferenceRecord {
    pub fn encode(&self) -> [u8; 16] {
        let mut encoded = [0u8; 16];
        encoded[..8].copy_from_slice(&self.v.to_be_bytes());
        encoded[8..].copy_from_slice(&self.x.to_be_bytes());
        encoded
    }
}

pub fn record(index: u64) -> Result<ReferenceRecord> {
    match index {
        0..=511 => Ok(ReferenceRecord {
            index,
            v: 1,
            x: index
                .checked_mul(2)
                .and_then(|value| value.checked_add(1))
                .ok_or_else(|| "REF-0 x overflow".to_owned())?,
        }),
        512..=1024 => Ok(ReferenceRecord {
            index,
            v: 2,
            x: (index - 512)
                .checked_mul(2)
                .ok_or_else(|| "REF-0 x overflow".to_owned())?,
        }),
        _ => Err(format!("REF-0 index out of range: {index}")),
    }
}

pub fn records() -> Result<Vec<ReferenceRecord>> {
    (0..RECORD_COUNT).map(record).collect()
}

pub fn record_bytes() -> Result<Vec<u8>> {
    let capacity = usize::try_from(RECORD_BYTES)
        .map_err(|_| "REF-0 byte count does not fit usize".to_owned())?;
    let mut encoded = Vec::with_capacity(capacity);
    for item in records()? {
        encoded.extend_from_slice(&item.encode());
    }
    if encoded.len() as u64 != RECORD_BYTES {
        return Err("REF-0 encoded byte count mismatch".to_owned());
    }
    Ok(encoded)
}

pub fn shell_bytes() -> Vec<u8> {
    let mut encoded = Vec::with_capacity(9);
    encoded.extend_from_slice(&(REFERENCE_SHELL.len() as u32).to_be_bytes());
    encoded.extend_from_slice(REFERENCE_SHELL.as_bytes());
    encoded
}

pub fn hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(DIGITS[(byte >> 4) as usize] as char);
        output.push(DIGITS[(byte & 0x0f) as usize] as char);
    }
    output
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Ref0Arithmetic {
    primitive_count: u64,
    duplicate_indices: Vec<u64>,
    endpoint_rows: Vec<(u64, BigInt, BigInt)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StatusCounts {
    pub nonprimitive_duplicate: u64,
    pub complete: u64,
    pub one_large_prime: u64,
    pub two_large_prime: u64,
    pub rejected: u64,
    pub invalid: u64,
    pub unresolved: u64,
}

impl StatusCounts {
    fn add(&mut self, status: u8) -> Result<()> {
        let slot = match status {
            0 => &mut self.nonprimitive_duplicate,
            1 => &mut self.complete,
            2 => &mut self.one_large_prime,
            3 => &mut self.two_large_prime,
            4 => &mut self.rejected,
            5 => &mut self.invalid,
            6 => &mut self.unresolved,
            _ => return Err(format!("unknown REF-0 status {status}")),
        };
        *slot = slot
            .checked_add(1)
            .ok_or_else(|| "REF-0 status counter overflow".to_owned())?;
        Ok(())
    }

    pub fn total(&self) -> u64 {
        self.nonprimitive_duplicate
            + self.complete
            + self.one_large_prime
            + self.two_large_prime
            + self.rejected
            + self.invalid
            + self.unresolved
    }

    pub fn retained(&self) -> u64 {
        self.complete + self.one_large_prime + self.two_large_prime
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReferenceOutput {
    pub candidate_records: Vec<u8>,
    pub candidate_records_sha256: [u8; 32],
    pub candidate_shard_sha256: [u8; 32],
    pub candidate_shell_sha256: [u8; 32],
    pub disposition_records: Vec<u8>,
    pub disposition_records_sha256: [u8; 32],
    pub disposition_shard_sha256: [u8; 32],
    pub disposition_shell_sha256: [u8; 32],
    pub certificates: Vec<u8>,
    pub status_counts: StatusCounts,
}

fn arithmetic() -> Result<Ref0Arithmetic> {
    let (p, t, _) = p192_parameters()?;
    let mut primitive_count = 0u64;
    let mut seen_primitive = BTreeSet::new();
    let mut duplicate_indices = Vec::new();
    let mut endpoint_rows = Vec::new();
    for item in records()? {
        let v = BigInt::from(item.v);
        let u = centered_u(item.x, item.v, &t)?;
        let content = u.gcd(&v).abs();
        if primitive(&u, &v) {
            primitive_count += 1;
            if !seen_primitive.insert((u.to_string(), v.to_string())) {
                return Err(format!(
                    "primitive REF-0 record {} duplicates an earlier primitive",
                    item.index
                ));
            }
        } else {
            if content.is_zero() {
                return Err("zero REF-0 content".to_owned());
            }
            let normalized = (&u / &content, &v / &content);
            if !seen_primitive.contains(&(normalized.0.to_string(), normalized.1.to_string())) {
                return Err(format!(
                    "nonprimitive REF-0 record {} has no earlier primitive representative",
                    item.index
                ));
            }
            duplicate_indices.push(item.index);
        }
        if matches!(item.index, 0 | 511 | 512 | 1024) {
            let alpha_norm = norm(&u, &v, &t, &p);
            endpoint_rows.push((item.index, u, alpha_norm));
        }
    }
    Ok(Ref0Arithmetic {
        primitive_count,
        duplicate_indices,
        endpoint_rows,
    })
}

pub(crate) fn candidate_roots(bytes: &[u8]) -> ([u8; 32], [u8; 32]) {
    let shell = shell_bytes();
    let mut shard_preimage = b"P192-WCM-CANDIDATE-SHARD-v1\0".to_vec();
    shard_preimage.extend_from_slice(&shell);
    shard_preimage.extend_from_slice(&0u64.to_be_bytes());
    shard_preimage.extend_from_slice(&RECORD_COUNT.to_be_bytes());
    shard_preimage.extend_from_slice(bytes);
    let shard = sha256(&shard_preimage);
    let mut shell_preimage = b"P192-WCM-SHELL-ROOT-v1\0".to_vec();
    shell_preimage.extend_from_slice(&shell);
    shell_preimage.extend_from_slice(&SHARD_COUNT.to_be_bytes());
    shell_preimage.extend_from_slice(&0u64.to_be_bytes());
    shell_preimage.extend_from_slice(&RECORD_COUNT.to_be_bytes());
    shell_preimage.extend_from_slice(&shard);
    (shard, sha256(&shell_preimage))
}

pub(crate) fn disposition_roots(bytes: &[u8]) -> ([u8; 32], [u8; 32]) {
    let shell = shell_bytes();
    let mut shard_preimage = b"P192-WCM-DISPOSITION-SHARD-v1\0".to_vec();
    shard_preimage.extend_from_slice(&shell);
    shard_preimage.extend_from_slice(&0u64.to_be_bytes());
    shard_preimage.extend_from_slice(&0u64.to_be_bytes());
    shard_preimage.extend_from_slice(&RECORD_COUNT.to_be_bytes());
    shard_preimage.extend_from_slice(bytes);
    let shard = sha256(&shard_preimage);
    let mut shell_preimage = b"P192-WCM-DISPOSITION-SHELL-ROOT-v1\0".to_vec();
    shell_preimage.extend_from_slice(&shell);
    shell_preimage.extend_from_slice(&SHARD_COUNT.to_be_bytes());
    shell_preimage.extend_from_slice(&0u64.to_be_bytes());
    shell_preimage.extend_from_slice(&0u64.to_be_bytes());
    shell_preimage.extend_from_slice(&RECORD_COUNT.to_be_bytes());
    shell_preimage.extend_from_slice(&shard);
    (shard, sha256(&shell_preimage))
}

fn validate_fixed_reference() -> Result<()> {
    let summary = arithmetic()?;
    if summary.primitive_count != 769
        || summary.duplicate_indices != (513..=1023).step_by(2).collect::<Vec<_>>()
    {
        return Err("REF-0 primitive/duplicate inventory mismatch".to_owned());
    }
    // These are mandatory algebraic controls, not a scientific BOX census.
    d23_control()?;
    p192_sqrt_discriminant_control()?;
    Ok(())
}

fn assemble_reference(mut classified: Vec<Option<Disposition>>) -> Result<ReferenceOutput> {
    let candidate_records = record_bytes()?;
    let candidate_records_sha256 = sha256(&candidate_records);
    let (candidate_shard_sha256, candidate_shell_sha256) = candidate_roots(&candidate_records);
    let mut dispositions = Vec::new();
    let mut certificate_entries = Vec::new();
    let mut status_counts = StatusCounts {
        nonprimitive_duplicate: 0,
        complete: 0,
        one_large_prime: 0,
        two_large_prime: 0,
        rejected: 0,
        invalid: 0,
        unresolved: 0,
    };
    for item in records()? {
        let disposition = classified[item.index as usize]
            .take()
            .ok_or_else(|| format!("missing REF-0 processing index: {}", item.index))?;
        let status = disposition.status();
        status_counts.add(status)?;
        dispositions.push(status);
        dispositions.extend_from_slice(&item.index.to_be_bytes());
        if let Some(certificate) = disposition.retained() {
            dispositions.push(1);
            dispositions.extend_from_slice(&certificate.candidate_sha256);
            dispositions.extend_from_slice(&certificate.factorization_sha256);
            certificate_entries.push((item.index, certificate.clone()));
        } else {
            dispositions.push(0);
        }
    }
    if status_counts.total() != RECORD_COUNT {
        return Err("REF-0 status counts do not cover every record".to_owned());
    }
    if status_counts.invalid != 0 || status_counts.unresolved != 0 {
        return Err(format!(
            "REF-0 cannot complete with invalid={} unresolved={}",
            status_counts.invalid, status_counts.unresolved
        ));
    }
    if status_counts.retained() != certificate_entries.len() as u64 {
        return Err("REF-0 certificate inventory mismatch".to_owned());
    }
    let disposition_records_sha256 = sha256(&dispositions);
    let (disposition_shard_sha256, disposition_shell_sha256) = disposition_roots(&dispositions);

    let mut certificates = b"P192-WCM-REF-CERT-STORE-v1\0".to_vec();
    debug_assert_eq!(certificates.len(), 27);
    certificates.extend_from_slice(&(certificate_entries.len() as u64).to_be_bytes());
    for (index, certificate) in certificate_entries {
        certificates.extend_from_slice(&index.to_be_bytes());
        certificates.extend_from_slice(
            &u32::try_from(certificate.candidate_bytes.len())
                .map_err(|_| "candidate certificate exceeds u32".to_owned())?
                .to_be_bytes(),
        );
        certificates.extend_from_slice(&certificate.candidate_bytes);
        certificates.extend_from_slice(
            &u32::try_from(certificate.factorization_bytes.len())
                .map_err(|_| "factorization certificate exceeds u32".to_owned())?
                .to_be_bytes(),
        );
        certificates.extend_from_slice(&certificate.factorization_bytes);
    }
    Ok(ReferenceOutput {
        candidate_records,
        candidate_records_sha256,
        candidate_shard_sha256,
        candidate_shell_sha256,
        disposition_records: dispositions,
        disposition_records_sha256,
        disposition_shard_sha256,
        disposition_shell_sha256,
        certificates,
        status_counts,
    })
}

pub fn build_reference(
    factor_base: &[Entry],
    factor_base_sha256: &[u8; 32],
    algebraic_bound: u64,
    large_prime_bound: u64,
) -> Result<ReferenceOutput> {
    let processing_order = (0..RECORD_COUNT).collect::<Vec<_>>();
    build_reference_with_order(
        factor_base,
        factor_base_sha256,
        algebraic_bound,
        large_prime_bound,
        &processing_order,
    )
}

pub(crate) fn build_reference_with_order(
    factor_base: &[Entry],
    factor_base_sha256: &[u8; 32],
    algebraic_bound: u64,
    large_prime_bound: u64,
    processing_order: &[u64],
) -> Result<ReferenceOutput> {
    if processing_order.len() as u64 != RECORD_COUNT {
        return Err("REF-0 processing order has the wrong length".to_owned());
    }
    validate_fixed_reference()?;
    let (field_norm, trace, discriminant) = p192_parameters()?;
    let mut classified = vec![None; RECORD_COUNT as usize];
    for &index in processing_order {
        let slot = classified
            .get_mut(index as usize)
            .ok_or_else(|| format!("REF-0 processing index out of range: {index}"))?;
        if slot.is_some() {
            return Err(format!("duplicate REF-0 processing index: {index}"));
        }
        let item = record(index)?;
        *slot = Some(certificate::classify(
            &item,
            factor_base,
            factor_base_sha256,
            &trace,
            &field_norm,
            &discriminant,
            algebraic_bound,
            large_prime_bound,
        )?);
    }
    assemble_reference(classified)
}

/// Distinct prime-major segmented reference implementation.
///
/// The scalar path above completely factors one candidate before advancing.
/// This path initializes all 1025 residuals, traverses the factor base in
/// fixed segments, and applies each prime across the whole candidate segment.
pub(crate) fn build_reference_segmented(
    factor_base: &[Entry],
    factor_base_sha256: &[u8; 32],
    algebraic_bound: u64,
    large_prime_bound: u64,
) -> Result<ReferenceOutput> {
    validate_fixed_reference()?;
    let (field_norm, trace, discriminant) = p192_parameters()?;
    let mut prepared = Vec::with_capacity(RECORD_COUNT as usize);
    let mut classified = vec![None; RECORD_COUNT as usize];
    for item in records()? {
        match certificate::prepare_candidate(&item, &trace, &field_norm)? {
            Some(candidate) => prepared.push(Some(candidate)),
            None => {
                prepared.push(None);
                classified[item.index as usize] = Some(Disposition::NonprimitiveDuplicate);
            }
        }
    }
    const FACTOR_SEGMENT_LENGTH: usize = 257;
    for segment in factor_base.chunks(FACTOR_SEGMENT_LENGTH) {
        for entry in segment {
            let mut marks = Vec::new();
            for (index, candidate) in prepared.iter_mut().enumerate() {
                let Some(candidate) = candidate else {
                    continue;
                };
                let matching_roots = entry
                    .roots
                    .iter()
                    .enumerate()
                    .filter(|(_, root)| {
                        (mod_u64(&candidate.u, entry.ell)
                            + mul_mod(mod_u64(&candidate.v, entry.ell), **root, entry.ell))
                        .is_multiple_of(entry.ell)
                    })
                    .map(|(root_index, _)| root_index)
                    .collect::<Vec<_>>();
                match matching_roots.as_slice() {
                    [] => {}
                    [root_index] => marks.push((index, *root_index)),
                    _ => candidate.invalid = true,
                }
            }
            for (index, root_index) in marks {
                let candidate = prepared[index]
                    .as_mut()
                    .ok_or_else(|| "segmented mark references a duplicate row".to_owned())?;
                let mut exponent = 0u32;
                while (&candidate.residual % entry.ell).is_zero() {
                    candidate.residual /= entry.ell;
                    exponent = exponent
                        .checked_add(1)
                        .ok_or_else(|| "segmented small-prime exponent overflow".to_owned())?;
                }
                if exponent == 0
                    || (entry.kind == 1 && (entry.roots.len() != 1 || exponent != 1))
                    || (entry.kind == 0 && entry.roots.len() != 2)
                    || !matches!(entry.kind, 0 | 1)
                {
                    candidate.invalid = true;
                    continue;
                }
                let signed_exponent = match (entry.kind, root_index) {
                    (1, 0) => BigInt::one(),
                    (0, 0) => BigInt::from(exponent),
                    (0, 1) => -BigInt::from(exponent),
                    _ => {
                        candidate.invalid = true;
                        continue;
                    }
                };
                candidate.rational_factors.push(RationalFactor {
                    prime: entry.ell,
                    exponent,
                });
                candidate.small_entries.push(SmallEntry {
                    ell: entry.ell,
                    kind: entry.kind,
                    exponent: signed_exponent,
                    root: entry.roots[root_index],
                    positive_root: entry.roots[0],
                    negative_root: *entry.roots.last().ok_or_else(|| {
                        "segmented factor-base entry has no characteristic root".to_owned()
                    })?,
                });
            }
        }
    }
    for index in 0..RECORD_COUNT {
        if let Some(candidate) = prepared[index as usize].take() {
            classified[index as usize] = Some(certificate::finish_prepared(
                &record(index)?,
                candidate,
                factor_base_sha256,
                &trace,
                &field_norm,
                &discriminant,
                algebraic_bound,
                large_prime_bound,
            )?);
        }
    }
    assemble_reference(classified)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encoded(index: u64) -> String {
        hex(&record(index).unwrap().encode())
    }

    #[test]
    fn exact_record_and_shell_bytes() {
        assert_eq!(shell_bytes(), b"\0\0\0\x05REF-0");
        assert_eq!(hex(&shell_bytes()), "000000055245462d30");
        assert_eq!(encoded(0), "00000000000000010000000000000001");
        assert_eq!(encoded(1), "00000000000000010000000000000003");
        assert_eq!(encoded(511), "000000000000000100000000000003ff");
        assert_eq!(encoded(512), "00000000000000020000000000000000");
        assert_eq!(encoded(513), "00000000000000020000000000000002");
        assert_eq!(encoded(1024), "00000000000000020000000000000400");
        assert_eq!(record_bytes().unwrap().len(), 16_400);
    }

    #[test]
    fn exact_centered_endpoints_and_norms() {
        let summary = arithmetic().unwrap();
        let got = summary
            .endpoint_rows
            .into_iter()
            .map(|(index, u, alpha_norm)| (index, u.to_string(), alpha_norm.to_string()))
            .collect::<Vec<_>>();
        assert_eq!(
            got,
            vec![
                (
                    0,
                    "-15803701158356963603741338599".to_owned(),
                    "6027344765084027530636040308278493916237318129474216339879".to_owned(),
                ),
                (
                    511,
                    "-15803701158356963603741338088".to_owned(),
                    "6027344765084027530636040308278493916237318129474216601511".to_owned(),
                ),
                (
                    512,
                    "-31607402316713927207482677199".to_owned(),
                    "24109379060336110122544161233113975664949272517896865359515".to_owned(),
                ),
                (
                    1024,
                    "-31607402316713927207482676687".to_owned(),
                    "24109379060336110122544161233113975664949272517896865621659".to_owned(),
                ),
            ]
        );
    }

    #[test]
    fn primitive_and_duplicate_inventory_is_exact() {
        let summary = arithmetic().unwrap();
        assert_eq!(summary.primitive_count, 769);
        assert_eq!(summary.duplicate_indices.len(), 256);
        assert_eq!(
            summary.duplicate_indices,
            (513..=1023).step_by(2).collect::<Vec<_>>()
        );
    }

    #[test]
    fn candidate_root_preimages_are_domain_separated() {
        let bytes = record_bytes().unwrap();
        let (shard, shell) = candidate_roots(&bytes);
        assert_ne!(shard, sha256(&bytes));
        assert_ne!(shell, shard);
        assert_eq!(shard.len(), 32);
    }

    #[test]
    fn full_ref0_has_complete_canonical_dispositions() {
        let base = super::super::factor_base::build(65_521).unwrap();
        let output = build_reference(&base, &[0x5au8; 32], 65_521, 2_147_483_647).unwrap();
        let segmented =
            build_reference_segmented(&base, &[0x5au8; 32], 65_521, 2_147_483_647).unwrap();
        assert_eq!(segmented, output);
        assert_eq!(output.status_counts.total(), RECORD_COUNT);
        assert_eq!(output.status_counts.nonprimitive_duplicate, 256);
        assert_eq!(output.status_counts.invalid, 0);
        assert_eq!(output.status_counts.unresolved, 0);
        assert_eq!(output.candidate_records.len(), RECORD_BYTES as usize);
        assert_eq!(
            certificate_store_count_for_test(&output.certificates),
            output.status_counts.retained()
        );
    }

    fn certificate_store_count_for_test(bytes: &[u8]) -> u64 {
        let offset = b"P192-WCM-REF-CERT-STORE-v1\0".len();
        u64::from_be_bytes(bytes[offset..offset + 8].try_into().unwrap())
    }
}
