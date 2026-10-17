//! Native verification that the fixture-only rebuild keeps runtime code,
//! allocated data and loader layout unchanged. Build IDs are metadata.
use std::{collections::BTreeMap, fs};

fn word(bytes: &[u8], offset: usize, length: usize) -> u64 {
    let mut value = [0u8; 8];
    value[..length].copy_from_slice(&bytes[offset..offset + length]);
    u64::from_le_bytes(value)
}

#[derive(Debug, PartialEq, Eq)]
struct Section {
    kind: u64,
    flags: u64,
    address: u64,
    offset: u64,
    size: u64,
    alignment: u64,
}

fn sections(bytes: &[u8]) -> BTreeMap<String, Section> {
    assert_eq!(&bytes[..6], b"\x7fELF\x02\x01", "ELF64/LE required");
    let start = word(bytes, 40, 8) as usize;
    let size = word(bytes, 58, 2) as usize;
    let count = word(bytes, 60, 2) as usize;
    let strings = start + word(bytes, 62, 2) as usize * size;
    let strings = word(bytes, strings + 24, 8) as usize;
    let mut result = BTreeMap::new();
    for i in 0..count {
        let header = start + i * size;
        let flags = word(bytes, header + 8, 8);
        if flags & 2 == 0 {
            continue;
        }
        let name = strings + word(bytes, header, 4) as usize;
        let end = bytes[name..].iter().position(|&b| b == 0).unwrap() + name;
        let name = String::from_utf8(bytes[name..end].to_vec()).unwrap();
        let section = Section {
            kind: word(bytes, header + 4, 4),
            flags,
            address: word(bytes, header + 16, 8),
            offset: word(bytes, header + 24, 8),
            size: word(bytes, header + 32, 8),
            alignment: word(bytes, header + 48, 8),
        };
        assert!(result.insert(name, section).is_none());
    }
    result
}

fn program_headers(bytes: &[u8]) -> &[u8] {
    let offset = word(bytes, 32, 8) as usize;
    let length = word(bytes, 54, 2) as usize * word(bytes, 56, 2) as usize;
    &bytes[offset..offset + length]
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 3, "MEASURED_ELF REBUILT_ELF");
    let a = fs::read(&args[1]).unwrap();
    let b = fs::read(&args[2]).unwrap();
    let sa = sections(&a);
    let sb = sections(&b);
    assert_eq!(sa, sb, "allocated section layouts differ");
    assert_eq!(
        program_headers(&a),
        program_headers(&b),
        "loader layout differs"
    );
    for (name, section) in &sa {
        if section.kind == 8 || name == ".note.gnu.build-id" {
            continue;
        }
        let start = section.offset as usize;
        let end = start + section.size as usize;
        assert_eq!(
            &a[start..end],
            &b[start..end],
            "runtime section differs: {name}"
        );
    }
    println!(
        "{{\"format\":\"ELF64/LE\",\"allocated_sections\":{},\"program_headers_identical\":true,\"allocated_layouts_identical\":true,\"runtime_code_and_data_identical_except_build_id\":true,\"ignored_content_section\":\".note.gnu.build-id\"}}",
        sa.len()
    );
}
