// Shared tables, read from the analyzer output.
#let jobs = (("k600", "IC K600"), ("k3300", "IC K3300"))
#let arm-names = (
  ic: "IC (frozen producer)",
  rho_plain_fold: "rho, empty table (fold, Tier 1 worker)",
  bl_a_time_matched: "rho + table, (a) 650 s, seed …601",
  bl_a2_wall_equal: "rho + table, (a′) = IC wall, seed …605",
  bl_b_memory_matched: "rho + table, (b) = IC max RSS",
)
#let f(x, d: 2) = str(calc.round(x, digits: d))
#let ci(p, d: 2) = "[" + f(p.at(0), d: d) + ", " + f(p.at(1), d: d) + "]"
#let analysis(dir, job) = json(dir + "/" + job + "/analysis.json")
#let cpu(dir, job) = {
  let lines = read(dir + "/" + job + "/host/lscpu.txt").split("\n")
  let m = lines.filter(l => l.starts-with("Model name"))
  if m.len() > 0 { m.at(0).split(":").at(1).trim() } else { "?" }
}

#let run-section(run-label, dir, root: "") = [
  #for (job, job-label) in jobs [
    #let a = analysis(root + dir, job)
    == #run-label — #job-label job (#cpu(root + dir, job))
    All recovered scalars replayed with the generic library curve and agreeing
    across arms and with the hosted answers: *#a.all_scalars_verified_and_agreeing*
    (#a.targets targets).

    #table(
      columns: (auto, auto, auto, auto, auto),
      align: (left, right, left, right, left),
      table.header[*arm*][*median ms*][*95% CI*][*mean ms*][*95% CI*],
      ..for (key, label) in arm-names {
        if key in a.arms {
          let s = a.arms.at(key)
          (label, f(s.median_ms, d: 3), ci(s.median_ms_ci95, d: 3), f(s.mean_ms, d: 3), ci(s.mean_ms_ci95, d: 3))
        }
      }
    )

    #table(
      columns: (auto, auto, auto, auto, auto, auto),
      align: (left, right, left, right, left, right),
      table.header[*paired ratio*][*median of ratios*][*95% CI*][*ratio of means*][*95% CI*][*num. slower*],
      ..for (key, p) in a.paired {
        (p.numerator + " / " + p.denominator, f(p.median_of_paired_ratios), ci(p.median_of_paired_ratios_ci95),
         f(p.ratio_of_means), ci(p.ratio_of_means_ci95), str(p.targets_where_numerator_slower) + "/16")
      }
    )

    #let ic = a.setup.ic
    #table(
      columns: (auto, auto, auto, auto, auto),
      align: (left, right, right, right, left),
      table.header[*precomputation*][*wall s*][*CPU s*][*memory*][*notes*],
      ..("IC", f(ic.target_independent_ms / 1000, d: 1),
       if ic.time_v_user_seconds != none { ic.time_v_user_seconds } else { "—" },
       f(ic.time_v_max_rss_bytes / 1e9, d: 2) + " GB max RSS",
       "rank threads " + str(ic.threads)),
      ..for key in ("bl_a_time_matched", "bl_a2_wall_equal", "bl_b_memory_matched") {
        if key in a.setup {
          let s = a.setup.at(key)
          (arm-names.at(key), f(s.setup_ms / 1000, d: 1), f(s.setup_cpu_seconds, d: 0),
           f(s.table_bytes / 1e9, d: 3) + " GB table",
           str(s.table_distinct_entries) + " entries, " + f(s.table_walk_rate_steps_per_second / 1e6, d: 2) + " M steps/s, stop: " + s.table_stopped_by)
        }
      }
    )
  ]
]

