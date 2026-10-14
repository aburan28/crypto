# Paper figures, 2026-10-14

SVG and 150-dpi PNG copies of the fourteen figures in the ECC2K-130 Blackwell
paper (`paper/ecc2k130/`, Draft 8), plus `figures.pdf` with all of them, one per
page. The figure sources are `paper/ecc2k130/figures/*.tex`; regenerate this
directory with `make -C paper/ecc2k130 figures-export`. Every value plotted is a
recorded measurement or a quantity derived from one in the paper, with its
evidence record cited there; nothing here is a new measurement.

| file | what it shows | paper |
|---|---|---|
| `fig-history` | ECC2K-130 rates 2009–2012 against today, log scale (GTX 295 bars per two-die card) | §2 |
| `fig-ledger` | the sixteen paired ledger steps on one RTX PRO 6000, marked by kind of change | §8 |
| `fig-gains` | where the speedup came from: 229x over 2010 split into chip and per-SM-clock; the ledger's 17.2x split by kind | §8 |
| `fig-kernel` | kernel anatomy: grid and waves, SM, the 16 walks of a thread, where state lives | §8 |
| `fig-step` | one step's data flow: forward pass, one inversion, reverse pass | §8 |
| `fig-clmad` | the packed 131-bit product on CLMAD and its logic-pipe reduction | §5 |
| `fig-tiling` | memory layout: untiled vs 256-worker tiles, the compact 17-byte record | §8 |
| `fig-roofline` | pipe utilisation of the 20.08 table-walk build | §8 |
| `fig-cycles` | cycle rules v2 and v3 of the table walk | §8 |
| `fig-gpus` | eleven GPUs: sigma-walk throughput and iterations per dollar | §8 |
| `fig-cpu` | CPU rates per core by instruction set | §9 |
| `fig-isa` | timeline of the instructions this arithmetic uses, GPU and CPU | §6 |
| `fig-bounds` | measured rates against derived ceilings per platform | §15 |
| `fig-solvetime` | expected days to an ECC2K-130 collision against fleet size, per device, with the cutoff-32 detection delay and the 2010 GTX 295 | §16 |

`solvetime.py` is the arithmetic behind §16's table and figure (expected
work, rho's spread, the cutoff-32 detection delay, list-price cost, and the
NIST expected work); `python3 solvetime.py` prints every value used.
