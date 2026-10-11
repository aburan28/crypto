//! Render the EXP7 all-pattern extension chart from its frozen JSON.
#[path = "../research/prime_fourier_flatness_exp7_20261007/render_all_patterns.rs"]
mod render;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    render::run()
}
