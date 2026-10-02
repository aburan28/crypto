/* isolab-launch: the in-container process launcher and timer.
 *
 * Bind-mounted into every container (static, so it runs in any image), it
 * either holds the container open (--hold) or runs the user's command and
 * measures exactly that process tree: CLOCK_MONOTONIC wall around fork/exec
 * to exit, and wait4(2) rusage for the child and everything it waited for.
 * The record is written as JSON to --record so the worker reads a timing
 * taken inside the box, free of container start-up and exec overhead.
 *
 *   isolab-launch --hold
 *   isolab-launch --record FILE [--stdin FILE] -- argv...
 */
#define _GNU_SOURCE
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

static volatile sig_atomic_t got_term = 0;
static pid_t child = 0;

static void on_term(int sig) {
    got_term = sig;
    if (child > 0) kill(child, SIGTERM);
}

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec / 1e9;
}

int main(int argc, char **argv) {
    const char *record = NULL, *stdin_file = NULL;
    int i = 1;
    for (; i < argc; i++) {
        if (!strcmp(argv[i], "--hold")) {
            signal(SIGTERM, on_term);
            signal(SIGINT, on_term);
            while (!got_term) pause();
            return 0;
        } else if (!strcmp(argv[i], "--record") && i + 1 < argc) {
            record = argv[++i];
        } else if (!strcmp(argv[i], "--stdin") && i + 1 < argc) {
            stdin_file = argv[++i];
        } else if (!strcmp(argv[i], "--")) {
            i++;
            break;
        } else {
            fprintf(stderr, "isolab-launch: unknown argument %s\n", argv[i]);
            return 125;
        }
    }
    if (i >= argc) {
        fprintf(stderr, "isolab-launch: no command\n");
        return 125;
    }
    signal(SIGTERM, on_term);
    signal(SIGINT, on_term);
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    child = fork();
    if (child < 0) {
        perror("fork");
        return 125;
    }
    if (child == 0) {
        if (stdin_file) {
            int fd = open(stdin_file, O_RDONLY);
            if (fd >= 0) { dup2(fd, 0); close(fd); }
        } else {
            int fd = open("/dev/null", O_RDONLY);
            if (fd >= 0) { dup2(fd, 0); close(fd); }
        }
        execvp(argv[i], &argv[i]);
        fprintf(stderr, "isolab-launch: cannot exec %s: %s\n", argv[i], strerror(errno));
        _exit(127);
    }
    int status = 0;
    struct rusage ru;
    for (;;) {
        pid_t w = wait4(child, &status, 0, &ru);
        if (w == child) break;
        if (w < 0 && errno == EINTR) continue;
        perror("wait4");
        return 125;
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) / 1e9;
    int exit_code = -1, sig = 0;
    if (WIFEXITED(status)) exit_code = WEXITSTATUS(status);
    else if (WIFSIGNALED(status)) sig = WTERMSIG(status);
    if (record) {
        FILE *f = fopen(record, "w");
        if (f) {
            fprintf(f, "{\"wall_s\": %.9f, \"exit_code\": %d, \"signal\": %d, "
                       "\"user_s\": %.6f, \"sys_s\": %.6f, \"max_rss_kb\": %ld, "
                       "\"minor_faults\": %ld, \"major_faults\": %ld, "
                       "\"voluntary_switches\": %ld, \"involuntary_switches\": %ld, "
                       "\"start_monotonic\": %.9f, \"terminated_by_launcher\": %d}\n",
                    wall, exit_code, sig,
                    ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1e6,
                    ru.ru_stime.tv_sec + ru.ru_stime.tv_usec / 1e6,
                    ru.ru_maxrss, ru.ru_minflt, ru.ru_majflt, ru.ru_nvcsw, ru.ru_nivcsw,
                    ts0.tv_sec + ts0.tv_nsec / 1e9, (int)got_term);
            fclose(f);
        }
    }
    (void)now;
    if (sig) return 128 + sig;
    return exit_code;
}
