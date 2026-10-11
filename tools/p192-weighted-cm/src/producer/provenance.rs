use super::Result;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BuildProvenance {
    pub protocol_commit: &'static str,
    pub source_commit: &'static str,
    pub source_dirty: bool,
    pub git_metadata_present: bool,
    pub build_profile: &'static str,
}

fn full_lower_hex(value: &str) -> bool {
    value.len() == 40
        && value.bytes().any(|byte| byte != b'0')
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

pub fn embedded() -> Result<BuildProvenance> {
    let protocol_commit = option_env!("P192_WCM_PROTOCOL_COMMIT")
        .ok_or_else(|| "binary has no embedded protocol commit".to_owned())?;
    let source_commit = option_env!("P192_WCM_SOURCE_COMMIT")
        .ok_or_else(|| "binary has no embedded source commit".to_owned())?;
    if !full_lower_hex(protocol_commit) {
        return Err("embedded protocol commit is not full lower-case 40-hex".to_owned());
    }
    if !full_lower_hex(source_commit) {
        return Err("embedded source commit is not full lower-case 40-hex".to_owned());
    }
    let source_dirty = match option_env!("P192_WCM_SOURCE_DIRTY") {
        Some("true") => true,
        Some("false") => false,
        _ => return Err("binary has no valid embedded dirty-state flag".to_owned()),
    };
    let git_metadata_present = match option_env!("P192_WCM_GIT_METADATA_PRESENT") {
        Some("true") => true,
        Some("false") => false,
        _ => return Err("binary has no valid embedded Git-metadata flag".to_owned()),
    };
    let build_profile = option_env!("P192_WCM_BUILD_PROFILE")
        .ok_or_else(|| "binary has no embedded Cargo build profile".to_owned())?;
    Ok(BuildProvenance {
        protocol_commit,
        source_commit,
        source_dirty,
        git_metadata_present,
        build_profile,
    })
}

pub fn clean() -> Result<BuildProvenance> {
    let provenance = embedded()?;
    if provenance.source_dirty {
        return Err("preflight refuses a binary built from a dirty source tree".to_owned());
    }
    if !provenance.git_metadata_present {
        return Err("preflight refuses a binary built without Git metadata".to_owned());
    }
    if provenance.build_profile != "release" {
        return Err("preflight requires an admission-eligible release build".to_owned());
    }
    Ok(provenance)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_script_embeds_strict_provenance() {
        let provenance = embedded().unwrap();
        assert!(full_lower_hex(provenance.protocol_commit));
        assert!(full_lower_hex(provenance.source_commit));
        assert!(!provenance.build_profile.is_empty());
        assert!(provenance.git_metadata_present);
        assert!(!full_lower_hex("0000000000000000000000000000000000000000"));
        assert_eq!(
            clean().is_ok(),
            !provenance.source_dirty
                && provenance.git_metadata_present
                && provenance.build_profile == "release"
        );
    }
}
