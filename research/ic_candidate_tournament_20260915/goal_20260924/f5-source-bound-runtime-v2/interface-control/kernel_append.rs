
/// Correctness-only interface probe. Never calls execute_job or curve/PDP code.
pub fn inspect_job_schema(input: &str) -> Result<serde_json::Value, String> {
    let job: Job = serde_json::from_str(input).map_err(|error| error.to_string())?;
    Ok(json!({
        "mode": job.mode,
        "degree": job.degree,
        "curve_a": job.curve_a,
        "public_targets": job.public_targets,
        "target_seeds": job.target_seeds,
        "algorithm_seed": job.algorithm_seed,
        "factor_base": job.factor_base,
        "config": job.config,
        "exclusive_phases": job.exclusive_phases
    }))
}
