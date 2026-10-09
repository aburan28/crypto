#include <stdio.h>

int main(int argc, char **argv) {
    if (argc != 6) return 2;
    FILE *output = fopen(argv[5], "w");
    if (!output) return 3;
    fputs("[]\n", output);
    fclose(output);
    return 0;
}
