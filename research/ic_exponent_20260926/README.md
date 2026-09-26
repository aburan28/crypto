# Ledger §20: the Koblitz thread's exponent against batch rho

The question: at `k = 32` targets, with every phase priced on current
`main`, how does the Koblitz collection thread's `m = 3` method compare
with batch rho at the same `k`, and how does the ratio scale with `r`?
The write-up is §20 of
`research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`, and the
page panel is `#koblitz-exponent-20260926` in
`docs/index-calculus-scoreboard.html`.

## Result

- **No crossing.**
  - Across nine sizes from `r = 2^18` to `2^47.2`, the thread sits
    4.86× [4.21, 5.55] batch rho at its best (`K_0/GF(2^41)`, `2^39`).
  - It sits 24.8× at `2^18` and 13.6× at `2^47.2`.
- **The page's `r^{1/6}` law, carried from §19.5's 6.38×, is falsified at
  the small end.** It predicted 0.8–1.2× below `2^25`, and the thread
  measures 12.7–24.8× there.
- **The frozen phase model's minimum is not confirmed.** It predicted
  2.85× at `2^28`; the measured minimum is 4.86×, eleven bits higher.
- **The top end's exponent is not decided.** The declared fit gives
  `β = 0.137 [−0.10, 0.37]`, which contains both `1/6` and the model's
  `0.141`.
- **The cause at the small end** is the descent, whose cost per target
  exceeds batch rho's below about `2^37` and does not amortise over `k`.
  On top of it come the workflow's constructions in big-integer
  arithmetic, 7–70% of `S`.
- **The frozen headline re-priced (Control 2).** It reproduces every
  frozen count and reads 12.53× batch rho, against §19.5's 6.38×.

Class: **accounting**.

## Files

| file | what |
|:--|:--|
| `PROTOCOL.md` | the declared protocol (v1), and Amendment 1, before any declared run |
| `predict.py`, `prediction.json` | the two predictions, from frozen constants, before any run |
| `make_params.py` | the recipe rules: one parameter file per size, column count, descent and seed set |
| `run.py` | the harness: manifest, sweep, measurement, controls; resumable, never overwrites |
| `control.py` | Control 1: a price report against `ic workflow`'s report and relation files |
| `analyse.py`, `analysis.json` | means, intervals, the fit and the graded targets |
| `render_rows.py` | the note's and the page's table rows, printed from `analysis.json` |
| `host.json` | the host manifest (AGENTS.md §10) |
| `runs/<curve>/sweep/` | the seed set `W` sweep: every `(columns, descent)` report and `chosen.json` |
| `runs/<curve>/measure/` | `M1`–`M4`: price reports (a `-double` rerun where the spread rule fired), workflow reports, control verdicts, `uptime` before and after |
| `runs/controls/` | the frozen headline and the thread's own `n = 41` and `n = 53` recipes |
| `runs/run.log` | the harness's log |

## Reproducing

    cargo build --release --bin ic
    cd research/ic_exponent_20260926
    python3 run.py all                 # RAYON_NUM_THREADS=1 under taskset -c 2, resumable
    python3 analyse.py > analysis.json
    python3 render_rows.py md

Counts are deterministic and reproduce on any host. The prices are this
host's: one x86-64 cloud container, Intel Xeon at 2.10 GHz with AVX-512,
rustc 1.94.1.
