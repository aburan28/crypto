walk(if type == "object" then
  (if has("ns_per_op") then del(.gae, .ns_per_op) else . end)
  | del(.wall_ns, .cpu_at_start, .cpu_at_end, .phases_ns, .framework_total_gae_before_solver_removal, .solve_wall_ns, .cpu_ns, .run_delay_ns, .timeslices, .placement, .hw_solve)
else . end)
