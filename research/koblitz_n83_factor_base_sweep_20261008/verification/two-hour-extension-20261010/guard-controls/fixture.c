#include <stdio.h>
#include <unistd.h>
#ifndef CONTROL_SLEEP
#define CONTROL_SLEEP 0
#endif
int main(int argc, char **argv) {
  if (argc != 6) return 2;
  sleep(CONTROL_SLEEP);
  FILE *f = fopen(argv[5], "w");
  if (!f) return 3;
  fputs("{\"schema\":\"n83.factored-s4-capacity-worker/v1\",\"status\":\"PASS_model_construction_only\",\"orbit_columns\":64,\"memory_cgroup_limit_bytes\":4294967296,\"memory_cgroup_swap_limit_bytes\":0,\"sat_variables\":1,\"sat_clauses\":1,\"solver_search_executed\":false,\"model_lifting_executed\":false,\"total_index_calculus_runtime_ms\":null,\"selected_best_total_runtime\":null}\n", f);
  return fclose(f) ? 4 : 0;
}
