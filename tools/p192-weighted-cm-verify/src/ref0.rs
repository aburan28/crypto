use std::{collections::BTreeMap, fs, path::Path};

use serde::{Deserialize, Serialize};

use crate::{
    encoding::{
        candidate_shard_digest, candidate_shell_digest, disposition_shard_digest,
        disposition_shell_digest, encode_candidate_record, encode_disposition_record,
        parse_canonical_json, sha256, ShardDigest,
    },
    factor_base::VerifiedFactorBase,
    ideal::evaluate_candidate,
    params::{CURVE_UID, REF0_COUNT, REF0_SHELL_ID},
    Result,
};

const CERTIFICATE_STORE_DOMAIN: &[u8] = b"P192-WCM-REF-CERT-STORE-v1\0";
pub const REFERENCE_SCHEMA: &str = "p192-wcm-reference-box-v1";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArtifactRecord {
    pub path: String,
    pub byte_length: u64,
    pub sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bounds {
    pub v_min: u64,
    pub v_max: u64,
    pub x_min: u64,
    pub x_max_inclusive: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CandidateShardRecord {
    pub shard_index: u64,
    pub first_candidate_index: u64,
    pub record_count: u64,
    pub records_sha256: String,
    pub candidate_shard_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DispositionShardRecord {
    pub shard_index: u64,
    pub first_candidate_index: u64,
    pub record_count: u64,
    pub byte_offset: u64,
    pub byte_length: u64,
    pub records_sha256: String,
    pub disposition_shard_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CertificateRecord {
    pub path: String,
    pub byte_length: u64,
    pub sha256: String,
    pub record_count: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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
    fn increment(&mut self, status: u8) -> Result<()> {
        let slot = usize::from(status);
        if slot >= 7 {
            return Err("unknown disposition status".to_owned());
        }
        let mut values = self.values();
        values[slot] = values[slot]
            .checked_add(1)
            .ok_or_else(|| "status count overflow".to_owned())?;
        *self = Self::from_values(values);
        Ok(())
    }

    pub fn total(&self) -> u64 {
        self.values().into_iter().sum()
    }

    fn values(&self) -> [u64; 7] {
        [
            self.nonprimitive_duplicate,
            self.complete,
            self.one_large_prime,
            self.two_large_prime,
            self.rejected,
            self.invalid,
            self.unresolved,
        ]
    }

    fn from_values(values: [u64; 7]) -> Self {
        let [nonprimitive_duplicate, complete, one_large_prime, two_large_prime, rejected, invalid, unresolved] =
            values;
        Self {
            nonprimitive_duplicate,
            complete,
            one_large_prime,
            two_large_prime,
            rejected,
            invalid,
            unresolved,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceBox {
    pub schema: String,
    pub experiment_id: String,
    pub protocol_version: u64,
    pub protocol_commit: String,
    pub source_commit: String,
    pub curve_uid: String,
    pub factor_base_sha256: String,
    pub shell_id: String,
    pub bounds: Bounds,
    pub candidate_count: u64,
    pub shard_count: u64,
    pub candidate_records: ArtifactRecord,
    pub candidate_shards: Vec<CandidateShardRecord>,
    pub candidate_shell_sha256: String,
    pub disposition_records: ArtifactRecord,
    pub disposition_shards: Vec<DispositionShardRecord>,
    pub disposition_shell_sha256: String,
    pub certificates: CertificateRecord,
    pub status_counts: StatusCounts,
    pub complete: bool,
}

#[derive(Clone, Debug)]
pub struct RegeneratedReference {
    pub candidate_records: Vec<u8>,
    pub disposition_records: Vec<u8>,
    pub certificate_store: Vec<u8>,
    pub candidate_shard_sha256: [u8; 32],
    pub candidate_shell_sha256: [u8; 32],
    pub disposition_shard_sha256: [u8; 32],
    pub disposition_shell_sha256: [u8; 32],
    pub status_counts: StatusCounts,
    pub retained_count: u64,
}

#[derive(Clone, Debug)]
pub struct VerifiedReference {
    pub manifest: ReferenceBox,
    pub regenerated: RegeneratedReference,
    pub manifest_bytes: Vec<u8>,
}

pub fn pair_at(index: u64) -> Result<(u64, u64)> {
    match index {
        0..=511 => Ok((1, 2 * index + 1)),
        512..=1024 => Ok((2, 2 * (index - 512))),
        _ => Err(format!("REF-0 candidate index {index} is out of range")),
    }
}

fn append_store_entry(
    output: &mut Vec<u8>,
    index: u64,
    candidate: &[u8],
    factorization: &[u8],
) -> Result<()> {
    output.extend_from_slice(&index.to_be_bytes());
    output.extend_from_slice(
        &u32::try_from(candidate.len())
            .map_err(|_| "candidate certificate exceeds u32 length".to_owned())?
            .to_be_bytes(),
    );
    output.extend_from_slice(candidate);
    output.extend_from_slice(
        &u32::try_from(factorization.len())
            .map_err(|_| "factorization certificate exceeds u32 length".to_owned())?
            .to_be_bytes(),
    );
    output.extend_from_slice(factorization);
    Ok(())
}

pub fn regenerate(factor_base: &VerifiedFactorBase) -> Result<RegeneratedReference> {
    let mut candidate_records = Vec::with_capacity(REF0_COUNT as usize * 16);
    let mut disposition_records = Vec::new();
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
    let mut retained_count = 0u64;

    for index in 0..REF0_COUNT {
        let (v_box, x_box) = pair_at(index)?;
        candidate_records.extend_from_slice(&encode_candidate_record(v_box, x_box));
        let evaluation = evaluate_candidate(index, v_box, x_box, factor_base)?;
        status_counts.increment(evaluation.status)?;
        let hashes = if let Some(certificates) = evaluation.certificates {
            retained_count += 1;
            let candidate_hash = sha256(&certificates.candidate);
            let factorization_hash = sha256(&certificates.factorization);
            append_store_entry(
                &mut certificate_entries,
                index,
                &certificates.candidate,
                &certificates.factorization,
            )?;
            Some((candidate_hash, factorization_hash))
        } else {
            None
        };
        disposition_records.extend_from_slice(&encode_disposition_record(
            evaluation.status,
            index,
            hashes,
        )?);
    }

    if candidate_records.len() != 16_400 {
        return Err("REF-0 regeneration did not produce 16400 candidate bytes".to_owned());
    }
    if status_counts.total() != REF0_COUNT
        || status_counts.nonprimitive_duplicate != 256
        || status_counts.invalid != 0
        || status_counts.unresolved != 0
    {
        return Err("REF-0 regenerated status conservation failed".to_owned());
    }

    let mut certificate_store =
        Vec::with_capacity(CERTIFICATE_STORE_DOMAIN.len() + 8 + certificate_entries.len());
    certificate_store.extend_from_slice(CERTIFICATE_STORE_DOMAIN);
    certificate_store.extend_from_slice(&retained_count.to_be_bytes());
    certificate_store.extend_from_slice(&certificate_entries);

    let candidate_shard_sha256 =
        candidate_shard_digest(REF0_SHELL_ID, 0, REF0_COUNT, &candidate_records)?;
    let candidate_shell_sha256 = candidate_shell_digest(
        REF0_SHELL_ID,
        &[ShardDigest {
            shard_index: 0,
            first_candidate_index: 0,
            record_count: REF0_COUNT,
            sha256: candidate_shard_sha256,
        }],
    )?;
    let disposition_shard_sha256 =
        disposition_shard_digest(REF0_SHELL_ID, 0, 0, REF0_COUNT, &disposition_records)?;
    let disposition_shell_sha256 = disposition_shell_digest(
        REF0_SHELL_ID,
        &[ShardDigest {
            shard_index: 0,
            first_candidate_index: 0,
            record_count: REF0_COUNT,
            sha256: disposition_shard_sha256,
        }],
    )?;

    Ok(RegeneratedReference {
        candidate_records,
        disposition_records,
        certificate_store,
        candidate_shard_sha256,
        candidate_shell_sha256,
        disposition_shard_sha256,
        disposition_shell_sha256,
        status_counts,
        retained_count,
    })
}

fn check_artifact(descriptor: &ArtifactRecord, expected_path: &str, bytes: &[u8]) -> Result<()> {
    if descriptor.path != expected_path
        || descriptor.byte_length != bytes.len() as u64
        || descriptor.sha256 != hex::encode(sha256(bytes))
    {
        return Err(format!("artifact descriptor mismatch for {expected_path}"));
    }
    Ok(())
}

fn read_regular(source: &Path, relative: &str) -> Result<Vec<u8>> {
    let path = source.join(relative);
    let metadata = fs::symlink_metadata(&path)
        .map_err(|error| format!("metadata {}: {error}", path.display()))?;
    if !metadata.file_type().is_file() || metadata.file_type().is_symlink() {
        return Err(format!(
            "{} is not a regular non-symlink file",
            path.display()
        ));
    }
    fs::read(&path).map_err(|error| format!("read {}: {error}", path.display()))
}

pub fn verify_reference_box(
    source: &Path,
    bytes: Vec<u8>,
    protocol_commit: &str,
    source_commit: &str,
    factor_base: &VerifiedFactorBase,
) -> Result<VerifiedReference> {
    let candidate_bytes = read_regular(source, "reference-box/candidate-records.bin")?;
    let disposition_bytes = read_regular(source, "reference-box/disposition.bin")?;
    let certificate_bytes = read_regular(source, "reference-box/certificates.bin")?;
    verify_reference_box_bytes(
        bytes,
        &candidate_bytes,
        &disposition_bytes,
        &certificate_bytes,
        protocol_commit,
        source_commit,
        factor_base,
    )
}

pub fn verify_reference_box_bytes(
    bytes: Vec<u8>,
    candidate_bytes: &[u8],
    disposition_bytes: &[u8],
    certificate_bytes: &[u8],
    protocol_commit: &str,
    source_commit: &str,
    factor_base: &VerifiedFactorBase,
) -> Result<VerifiedReference> {
    let value = parse_canonical_json(&bytes)?;
    let manifest: ReferenceBox =
        serde_json::from_value(value).map_err(|error| format!("reference-box schema: {error}"))?;
    if manifest.schema != REFERENCE_SCHEMA
        || manifest.experiment_id != "EXP-SCURVE-1a8daf"
        || manifest.protocol_version != 2
        || manifest.protocol_commit != protocol_commit
        || manifest.source_commit != source_commit
        || manifest.curve_uid != CURVE_UID
        || manifest.factor_base_sha256 != hex::encode(factor_base.sha256)
        || manifest.shell_id != REF0_SHELL_ID
        || manifest.bounds
            != (Bounds {
                v_min: 1,
                v_max: 2,
                x_min: 0,
                x_max_inclusive: 1024,
            })
        || manifest.candidate_count != REF0_COUNT
        || manifest.shard_count != 1
        || manifest.candidate_shards.len() != 1
        || manifest.disposition_shards.len() != 1
    {
        return Err("reference-box envelope differs from frozen REF-0".to_owned());
    }

    let regenerated = regenerate(factor_base)?;
    if candidate_bytes != regenerated.candidate_records {
        return Err("REF-0 candidate records differ from exact regeneration".to_owned());
    }
    if disposition_bytes != regenerated.disposition_records {
        return Err("REF-0 disposition records differ from exact regeneration".to_owned());
    }
    if certificate_bytes != regenerated.certificate_store {
        return Err("REF-0 certificate store differs from exact replay".to_owned());
    }
    check_artifact(
        &manifest.candidate_records,
        "reference-box/candidate-records.bin",
        candidate_bytes,
    )?;
    check_artifact(
        &manifest.disposition_records,
        "reference-box/disposition.bin",
        disposition_bytes,
    )?;
    if manifest.certificates.path != "reference-box/certificates.bin"
        || manifest.certificates.byte_length != certificate_bytes.len() as u64
        || manifest.certificates.sha256 != hex::encode(sha256(certificate_bytes))
        || manifest.certificates.record_count != regenerated.retained_count
    {
        return Err("reference certificate descriptor mismatch".to_owned());
    }

    let candidate_shard = &manifest.candidate_shards[0];
    if candidate_shard.shard_index != 0
        || candidate_shard.first_candidate_index != 0
        || candidate_shard.record_count != REF0_COUNT
        || candidate_shard.records_sha256 != hex::encode(sha256(candidate_bytes))
        || candidate_shard.candidate_shard_sha256 != hex::encode(regenerated.candidate_shard_sha256)
        || manifest.candidate_shell_sha256 != hex::encode(regenerated.candidate_shell_sha256)
    {
        return Err("reference candidate shard/root descriptor mismatch".to_owned());
    }
    let disposition_shard = &manifest.disposition_shards[0];
    if disposition_shard.shard_index != 0
        || disposition_shard.first_candidate_index != 0
        || disposition_shard.record_count != REF0_COUNT
        || disposition_shard.byte_offset != 0
        || disposition_shard.byte_length != disposition_bytes.len() as u64
        || disposition_shard.records_sha256 != hex::encode(sha256(disposition_bytes))
        || disposition_shard.disposition_shard_sha256
            != hex::encode(regenerated.disposition_shard_sha256)
        || manifest.disposition_shell_sha256 != hex::encode(regenerated.disposition_shell_sha256)
    {
        return Err("reference disposition shard/root descriptor mismatch".to_owned());
    }
    if manifest.status_counts != regenerated.status_counts
        || manifest.status_counts.total() != REF0_COUNT
        || !manifest.complete
    {
        return Err("reference status counters or completeness mismatch".to_owned());
    }
    Ok(VerifiedReference {
        manifest,
        regenerated,
        manifest_bytes: bytes,
    })
}

pub fn parse_certificate_store(bytes: &[u8]) -> Result<BTreeMap<u64, (Vec<u8>, Vec<u8>)>> {
    let mut offset = 0usize;
    if !bytes.starts_with(CERTIFICATE_STORE_DOMAIN) {
        return Err("certificate store domain mismatch".to_owned());
    }
    offset += CERTIFICATE_STORE_DOMAIN.len();
    let count = take_u64(bytes, &mut offset)?;
    let mut entries = BTreeMap::new();
    let mut previous = None;
    for _ in 0..count {
        let index = take_u64(bytes, &mut offset)?;
        if previous.is_some_and(|old| index <= old) {
            return Err("certificate store indices are not strictly increasing".to_owned());
        }
        let candidate_length = take_u32(bytes, &mut offset)? as usize;
        let candidate = take(bytes, &mut offset, candidate_length)?.to_vec();
        let factorization_length = take_u32(bytes, &mut offset)? as usize;
        let factorization = take(bytes, &mut offset, factorization_length)?.to_vec();
        entries.insert(index, (candidate, factorization));
        previous = Some(index);
    }
    if offset != bytes.len() || entries.len() as u64 != count {
        return Err("certificate store has a count mismatch or trailing bytes".to_owned());
    }
    Ok(entries)
}

fn take<'a>(bytes: &'a [u8], offset: &mut usize, length: usize) -> Result<&'a [u8]> {
    let end = offset
        .checked_add(length)
        .ok_or_else(|| "binary length overflow".to_owned())?;
    let result = bytes
        .get(*offset..end)
        .ok_or_else(|| "truncated binary artifact".to_owned())?;
    *offset = end;
    Ok(result)
}

fn take_u32(bytes: &[u8], offset: &mut usize) -> Result<u32> {
    let raw: [u8; 4] = take(bytes, offset, 4)?
        .try_into()
        .map_err(|_| "internal u32 slice error".to_owned())?;
    Ok(u32::from_be_bytes(raw))
}

fn take_u64(bytes: &[u8], offset: &mut usize) -> Result<u64> {
    let raw: [u8; 8] = take(bytes, offset, 8)?
        .try_into()
        .map_err(|_| "internal u64 slice error".to_owned())?;
    Ok(u64::from_be_bytes(raw))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::factor_base::{regenerate_entries, FactorBaseManifest};
    use crate::params::{ALGEBRAIC_BOUND, D_DEC, MAPPABLE_BOUND, P_DEC, T_DEC};

    #[test]
    fn exact_ref0_index_mapping() {
        assert_eq!(pair_at(0).unwrap(), (1, 1));
        assert_eq!(pair_at(511).unwrap(), (1, 1023));
        assert_eq!(pair_at(512).unwrap(), (2, 0));
        assert_eq!(pair_at(1024).unwrap(), (2, 1024));
        assert!(pair_at(1025).is_err());
    }

    #[test]
    fn empty_certificate_store_is_canonical() {
        let mut bytes = CERTIFICATE_STORE_DOMAIN.to_vec();
        bytes.extend_from_slice(&0u64.to_be_bytes());
        assert!(parse_certificate_store(&bytes).unwrap().is_empty());
        bytes.push(0);
        assert!(parse_certificate_store(&bytes).is_err());
    }

    #[test]
    fn full_reference_regeneration_conserves_and_replays_certificates() {
        let manifest = FactorBaseManifest {
            schema: crate::factor_base::FACTOR_BASE_SCHEMA.to_owned(),
            curve_uid: CURVE_UID.to_owned(),
            p: P_DEC.to_owned(),
            t: T_DEC.to_owned(),
            discriminant: D_DEC.to_owned(),
            algebraic_bound: ALGEBRAIC_BOUND,
            mappable_bound: MAPPABLE_BOUND,
            source_commit: "0".repeat(40),
            entries: regenerate_entries().unwrap(),
            entry_count: 0,
        };
        let mut manifest = manifest;
        manifest.entry_count = manifest.entries.len() as u64;
        let factor_base = VerifiedFactorBase {
            manifest,
            canonical_bytes: Vec::new(),
            sha256: [0x5a; 32],
        };
        let reference = regenerate(&factor_base).unwrap();
        assert_eq!(reference.candidate_records.len(), 16_400);
        assert_eq!(reference.status_counts.total(), REF0_COUNT);
        assert_eq!(reference.status_counts.nonprimitive_duplicate, 256);
        assert_eq!(reference.status_counts.invalid, 0);
        assert_eq!(reference.status_counts.unresolved, 0);
        let certificates = parse_certificate_store(&reference.certificate_store).unwrap();
        assert_eq!(certificates.len() as u64, reference.retained_count);

        let (v, x) = pair_at(512).unwrap();
        assert_eq!(
            evaluate_candidate(512, v, x, &factor_base).unwrap().status,
            4
        );
        for index in (513..=1023).step_by(2) {
            let (v, x) = pair_at(index).unwrap();
            assert_eq!(
                evaluate_candidate(index, v, x, &factor_base)
                    .unwrap()
                    .status,
                0
            );
        }
    }

    fn manifest_for(
        reference: &RegeneratedReference,
        factor_base: &VerifiedFactorBase,
    ) -> ReferenceBox {
        ReferenceBox {
            schema: REFERENCE_SCHEMA.to_owned(),
            experiment_id: "EXP-SCURVE-1a8daf".to_owned(),
            protocol_version: 2,
            protocol_commit: "a".repeat(40),
            source_commit: "b".repeat(40),
            curve_uid: CURVE_UID.to_owned(),
            factor_base_sha256: hex::encode(factor_base.sha256),
            shell_id: REF0_SHELL_ID.to_owned(),
            bounds: Bounds {
                v_min: 1,
                v_max: 2,
                x_min: 0,
                x_max_inclusive: 1024,
            },
            candidate_count: REF0_COUNT,
            shard_count: 1,
            candidate_records: ArtifactRecord {
                path: "reference-box/candidate-records.bin".to_owned(),
                byte_length: reference.candidate_records.len() as u64,
                sha256: hex::encode(sha256(&reference.candidate_records)),
            },
            candidate_shards: vec![CandidateShardRecord {
                shard_index: 0,
                first_candidate_index: 0,
                record_count: REF0_COUNT,
                records_sha256: hex::encode(sha256(&reference.candidate_records)),
                candidate_shard_sha256: hex::encode(reference.candidate_shard_sha256),
            }],
            candidate_shell_sha256: hex::encode(reference.candidate_shell_sha256),
            disposition_records: ArtifactRecord {
                path: "reference-box/disposition.bin".to_owned(),
                byte_length: reference.disposition_records.len() as u64,
                sha256: hex::encode(sha256(&reference.disposition_records)),
            },
            disposition_shards: vec![DispositionShardRecord {
                shard_index: 0,
                first_candidate_index: 0,
                record_count: REF0_COUNT,
                byte_offset: 0,
                byte_length: reference.disposition_records.len() as u64,
                records_sha256: hex::encode(sha256(&reference.disposition_records)),
                disposition_shard_sha256: hex::encode(reference.disposition_shard_sha256),
            }],
            disposition_shell_sha256: hex::encode(reference.disposition_shell_sha256),
            certificates: CertificateRecord {
                path: "reference-box/certificates.bin".to_owned(),
                byte_length: reference.certificate_store.len() as u64,
                sha256: hex::encode(sha256(&reference.certificate_store)),
                record_count: reference.retained_count,
            },
            status_counts: reference.status_counts.clone(),
            complete: true,
        }
    }

    fn test_factor_base() -> VerifiedFactorBase {
        let entries = regenerate_entries().unwrap();
        VerifiedFactorBase {
            manifest: FactorBaseManifest {
                schema: crate::factor_base::FACTOR_BASE_SCHEMA.to_owned(),
                curve_uid: CURVE_UID.to_owned(),
                p: P_DEC.to_owned(),
                t: T_DEC.to_owned(),
                discriminant: D_DEC.to_owned(),
                algebraic_bound: ALGEBRAIC_BOUND,
                mappable_bound: MAPPABLE_BOUND,
                source_commit: "b".repeat(40),
                entry_count: entries.len() as u64,
                entries,
            },
            canonical_bytes: Vec::new(),
            sha256: [0x5a; 32],
        }
    }

    #[test]
    fn semantic_stream_mutations_fail_even_with_rehashed_descriptors() {
        let factor_base = test_factor_base();
        let reference = regenerate(&factor_base).unwrap();
        let manifest = manifest_for(&reference, &factor_base);
        let bytes =
            crate::encoding::canonical_json(&serde_json::to_value(&manifest).unwrap()).unwrap();
        verify_reference_box_bytes(
            bytes,
            &reference.candidate_records,
            &reference.disposition_records,
            &reference.certificate_store,
            &"a".repeat(40),
            &"b".repeat(40),
            &factor_base,
        )
        .unwrap();

        let mut candidate = reference.candidate_records.clone();
        candidate[15] ^= 1;
        let mut mutated = manifest.clone();
        mutated.candidate_records.sha256 = hex::encode(sha256(&candidate));
        mutated.candidate_shards[0].records_sha256 = hex::encode(sha256(&candidate));
        let shard = candidate_shard_digest(REF0_SHELL_ID, 0, REF0_COUNT, &candidate).unwrap();
        mutated.candidate_shards[0].candidate_shard_sha256 = hex::encode(shard);
        mutated.candidate_shell_sha256 = hex::encode(
            candidate_shell_digest(
                REF0_SHELL_ID,
                &[ShardDigest {
                    shard_index: 0,
                    first_candidate_index: 0,
                    record_count: REF0_COUNT,
                    sha256: shard,
                }],
            )
            .unwrap(),
        );
        let mutated_json =
            crate::encoding::canonical_json(&serde_json::to_value(&mutated).unwrap()).unwrap();
        assert!(verify_reference_box_bytes(
            mutated_json,
            &candidate,
            &reference.disposition_records,
            &reference.certificate_store,
            &"a".repeat(40),
            &"b".repeat(40),
            &factor_base,
        )
        .is_err());

        let mut disposition = reference.disposition_records.clone();
        disposition[0] ^= 1;
        let mut mutated = manifest.clone();
        mutated.disposition_records.sha256 = hex::encode(sha256(&disposition));
        mutated.disposition_shards[0].records_sha256 = hex::encode(sha256(&disposition));
        let shard =
            disposition_shard_digest(REF0_SHELL_ID, 0, 0, REF0_COUNT, &disposition).unwrap();
        mutated.disposition_shards[0].disposition_shard_sha256 = hex::encode(shard);
        mutated.disposition_shell_sha256 = hex::encode(
            disposition_shell_digest(
                REF0_SHELL_ID,
                &[ShardDigest {
                    shard_index: 0,
                    first_candidate_index: 0,
                    record_count: REF0_COUNT,
                    sha256: shard,
                }],
            )
            .unwrap(),
        );
        let mutated_json =
            crate::encoding::canonical_json(&serde_json::to_value(&mutated).unwrap()).unwrap();
        assert!(verify_reference_box_bytes(
            mutated_json,
            &reference.candidate_records,
            &disposition,
            &reference.certificate_store,
            &"a".repeat(40),
            &"b".repeat(40),
            &factor_base,
        )
        .is_err());

        let mut certificates = reference.certificate_store.clone();
        let last = certificates.len() - 1;
        certificates[last] ^= 1;
        let mut mutated = manifest;
        mutated.certificates.sha256 = hex::encode(sha256(&certificates));
        let mutated_json =
            crate::encoding::canonical_json(&serde_json::to_value(&mutated).unwrap()).unwrap();
        assert!(verify_reference_box_bytes(
            mutated_json,
            &reference.candidate_records,
            &reference.disposition_records,
            &certificates,
            &"a".repeat(40),
            &"b".repeat(40),
            &factor_base,
        )
        .is_err());
    }
}
