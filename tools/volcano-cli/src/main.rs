fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let (value, code) = match volcano_cli::run(&args) {
        Ok(result) => result,
        Err(error) => (serde_json::json!({"status":"ERROR","error":error}), 2),
    };
    println!("{}", serde_json::to_string_pretty(&value).unwrap());
    std::process::exit(code);
}
