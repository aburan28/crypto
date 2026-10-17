//! Native disk-backed, resumable prime-field isogeny frontier and paired runs.
//! Linux: dynamically links the system's libsqlite3.so.0; no Python or FUSE.
#[cfg(target_os = "linux")]
#[path = "support/disk_sqlite.rs"]
mod disk_sqlite;
#[cfg(target_os = "linux")]
#[path = "support/isogeny_disk.rs"]
mod isogeny_disk;
#[cfg(target_os = "linux")]
fn main() {
    isogeny_disk::main();
}
#[cfg(not(target_os = "linux"))]
fn main() {
    panic!("This checkpoint example currently supports Linux with libsqlite3.so.0");
}
