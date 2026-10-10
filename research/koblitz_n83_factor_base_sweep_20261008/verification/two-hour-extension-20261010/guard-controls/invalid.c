#include <stdio.h>
int main(int argc, char **argv) {
  if (argc != 6) return 2;
  FILE *f = fopen(argv[5], "w");
  if (!f) return 3;
  fputs("{\"schema\":\"wrong-schema\",\"status\":\"PASS_model_construction_only\"}\n", f);
  return fclose(f) ? 4 : 0;
}
