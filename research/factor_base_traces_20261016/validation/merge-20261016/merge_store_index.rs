use crypto_lib::store::index::Index;
use std::path::Path;

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 4);
    let mut merged = Index::load(Path::new(&args[1]))?;
    let other = Index::load(Path::new(&args[2]))?;
    let shared = other
        .entries
        .keys()
        .filter(|k| merged.entries.contains_key(*k))
        .count();
    for (name, entry) in other.entries {
        if let Some(old) = merged.entries.get(&name) {
            if old != &entry {
                return Err(format!(
                    "Conflicting existing entry {name}; refusing to choose a side"
                ));
            }
        } else {
            merged.put(&name, entry, false)?;
        }
    }
    merged.save(Path::new(&args[3]))?;
    let reread = Index::load(Path::new(&args[3]))?;
    assert_eq!(merged, reread);
    println!("Preserved {} entries from both sides ({} shared); native Index::load/put/save round trip passed", merged.entries.len(), shared);
    Ok(())
}
