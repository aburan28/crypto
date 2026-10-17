//! Generic process launch only; scientific budgets remain in the native tool.
use std::{
    env,
    fs::OpenOptions,
    os::unix::process::CommandExt,
    process::{Command, Stdio},
};
unsafe extern "C" {
    fn setsid() -> i32;
}
fn main() -> std::io::Result<()> {
    let args: Vec<_> = env::args_os().collect();
    if args.len() < 4 {
        eprintln!("usage: detach LOG COMMAND [ARG ...]");
        std::process::exit(64);
    }
    let log = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[1])?;
    let mut command = Command::new(&args[2]);
    command
        .args(&args[3..])
        .stdin(Stdio::null())
        .stdout(log.try_clone()?)
        .stderr(log);
    unsafe {
        command.pre_exec(|| {
            if setsid() == -1 {
                Err(std::io::Error::last_os_error())
            } else {
                Ok(())
            }
        });
    }
    let child = command.spawn()?;
    println!("{}", child.id());
    Ok(())
}
