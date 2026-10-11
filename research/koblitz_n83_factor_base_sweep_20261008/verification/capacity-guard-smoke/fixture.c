#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char **argv) {
    if (argc != 6 || strcmp(argv[1], "primary-sat-build") != 0) return 2;
    if (atoi(argv[3]) == 256) { sleep(3); return 0; }
    if (atoi(argv[3]) == 600) {
        char *memory = malloc(256 * 1024 * 1024);
        if (!memory) return 5;
        for (size_t index = 0; index < 256 * 1024 * 1024; index += 4096) memory[index] = 1;
        sleep(3);
        return 0;
    }
    FILE *source = fopen("/sys/fs/cgroup/memory.max", "r");
    unsigned long long limit = 0;
    if (!source || fscanf(source, "%llu", &limit) != 1) return 3;
    fclose(source);
    FILE *output = fopen(argv[5], "w");
    if (!output) return 4;
    fprintf(output,
        "{\"schema\":\"n83.factored-s4-capacity-worker/v1\","
        "\"status\":\"PASS_model_construction_only\","
        "\"orbit_columns\":%d,\"memory_cgroup_limit_bytes\":%llu,"
        "\"memory_cgroup_swap_limit_bytes\":0,"
        "\"solver_search_executed\":false,\"model_lifting_executed\":false,"
        "\"total_index_calculus_runtime_ms\":null,"
        "\"selected_best_total_runtime\":null}\n", atoi(argv[3]), limit);
    fclose(output);
    return 0;
}
