# Ring-polymer quantum annealing of silicon

RPQA applied to a periodic silicon cell, with forces from a NEP machine-learning
potential through an ASE socket client. The method is described in

> A. Fiorentino and N. Marzari, *Quantum annealing for materials*,
> [arXiv:2606.03405](https://arxiv.org/abs/2606.03405) (2026)

```
./run.sh
```

About 19 minutes with 8 clients. The driver needs **pyNEP**
(<https://github.com/bigd4/PyNEP>) and ASE; if they are not in the `python` on
your PATH, point `PYTHON` at an interpreter that has them:

```
PYTHON=~/miniforge3/envs/nep_env/bin/python ./run.sh
```

## What it does

The ring polymer is propagated at 50 K while `lambdaqkin` is ramped down. Every
`pinning_interval` steps the `rpqa` smotion relaxes every bead, pins the replica
sitting in the deepest minimum, and hands control back to the dynamics, which
carries on around the frozen replica. The delocalization stage is not separate:
it is the first `start_step` steps, during which the smotion does nothing.

The target is the **Fd-3m diamond structure**, which for these 32 atoms in this
cell relaxes to -148.0318 eV, i.e. -4.625994 eV/atom, with coordination exactly
4 and a nearest-neighbour distance of 2.3714 A.

| | |
|---|---|
| system | 32 Si, periodic, cell from the structure file (2x2x1 conventional, a = 5.4765 A) |
| beads | 32 |
| temperature | 50 K |
| `lambdaqkin` | 1200 -> 2, `sqrtscale` |
| timestep / `tau` | 4 fs / 40 fs |
| `nmts` | 20 |
| delocalization | 1250 steps (5 ps) |
| `pinning_interval` | 1250 steps (5 ps), 8 events |
| total | 12500 steps = 50 ps |
| relaxation | `cg_rp`, 0.01 eV/A, 1 meV, 0.05 A, capped at 1000 steps |

## Choosing the `lambdaqkin` range

Nuclear quantum effects scale as `lambdaqkin * hbar^2 / M`, and silicon is
heavy: a value that delocalizes a light nucleus does nothing here. The useful
handle is the spread of the ring polymer about its centroid,
`r_g = hbar sqrt(lambdaqkin * beta / 12 m)`:

| `lambdaqkin` | r_g per component | in 3D | |
|---|---|---|---|
| 1200 | 1.86 A | 3.22 A | wider than the 2.37 A bond: fully delocalized |
| 300 | 0.93 A | 1.61 A | ~1 A, not small |
| 30 | 0.29 A | 0.51 A | still not small |
| **2** | **0.076 A** | **0.13 A** | << 1 A: essentially localized |

Start where the ring is wider than a bond, so replicas can sit in different
basins; end where the delocalization is far below a bond length, so they collapse
onto one. Ending at 30 rather than 2 leaves the replicas spread over ~4 eV and
they never agree on a minimum — the run then does not anneal, it just drifts.

## What a run looks like

`reference/` holds the full output of the shipped configuration (seed 31415).
The annealing shows up clearly in the pinning log:

```
   ev  lambda  pinned   best(eV)   vs Fd-3m/atom  spread(eV)  replicas at best
    1     929      14  -140.3122      +0.2412        4.35         1/32
    2     693      28  -140.5245      +0.2346        4.52         1/32
    3     492      25  -140.5890      +0.2326        4.34         1/32
    4     325      13  -140.9260      +0.2221        5.96         1/32
    5     192      11  -148.0316      +0.0000       10.83         1/32
    6      94      15  -148.0317       0.0000       10.77        12/32
    7      31      30  -148.0317       0.0000        0.0003      32/32
    8       2      30  -148.0317       0.0000        0.0003      32/32
    9       2      30  -148.0317       0.0000        0.0003      32/32
```

Four events grind at +0.23 eV/atom; then at `lambdaqkin` = 192 a single replica
finds the crystal, the pin moves to it, and over the next two events the whole
ring collapses onto it — spread 10.8 eV to 3e-4 eV, 1 replica to all 32. Each of
events 1 to 5 is a genuine improvement, which is what re-pinning is for: the
pinned replica is not a commitment, and a replica that later quenches deeper
takes over.

**This does not happen every time.** With a different seed and everything else
identical, the run annealed steadily (-139.47 to -145.48 eV, recovering 0.16 of
the 0.24 eV/atom it started above Fd-3m) and the replicas still collapsed onto a
single minimum, but that minimum was amorphous rather than the crystal. At these
settings, finding Fd-3m in 50 ps is roughly a coin flip. Production runs use
25000-50000 steps; if you need the crystal reliably rather than as a
demonstration, lengthen the schedule.

## The starting structures

`structures/nbeads_32/optimal_beads_*.xyz` hold four independent ring-polymer
configurations, 32 frames each, one per bead. `input.xml` uses
`optimal_beads_0.xyz`; change that one line to try the others.

**None of them contains the crystal.** Quenching all 128 replicas individually
at fixed cell:

| file | best replica | worst | best vs Fd-3m |
|---|---|---|---|
| `optimal_beads_0.xyz` | -141.32 eV | -136.54 | **+0.210 eV/atom** |
| `optimal_beads_1.xyz` | -140.09 | -135.58 | +0.248 |
| `optimal_beads_2.xyz` | -139.77 | -135.47 | +0.258 |
| `optimal_beads_3.xyz` | -139.77 | -136.20 | +0.258 |

So the annealing has to find a basin no starting replica is in. The smaller
`structures/N16/` set is a much easier problem — there the replicas already
quench to essentially the answer — and is kept for quick tests.

## Why `cg_rp` and not `lbfgs`

At high `lambdaqkin` the replicas are spread over several eV and reach their
minima at very different rates. `cg_rp` gives each bead its own convergence test
and step length; every other optimizer shares one test across all beads, so the
relaxation waits for the worst replica.

The difference is not marginal. With `lbfgs` this system ran 2000 iterations
with the energy frozen and the maximum force stuck at 1.7e-2 Ha/bohr, never
converging. With `cg_rp` the relaxations take 24 to 550 steps.

## Number of pinning events

Halving the events to four, at 2500-step intervals, with the same seed, total
steps and lambda schedule, ended 0.17 eV/atom above Fd-3m instead of reaching
it — and two of its four events produced no improvement at all. Its first
relaxation also hit the 1000-step cap, the only cap hit in any of these runs,
because at 2500-step intervals the first event comes later with the ring more
delocalized. Eight events cost about 7% more wall time than four. Measured once,
on one seed, so treat it as a hint rather than a result.

## Reading the output

`rpqa.rpqa_pin` — one line per pinning event: the replica pinned, its energy, the
spread across replicas, relaxation steps used, and the current `lambdaqkin`.

`rpqa.rpqa_inherent.out` — with `<inherent_data>True</inherent_data>`, the
quenched energy of every bead at every event, one column per bead, in eV, with
the structures in `rpqa.rpqa_inherent.pos_*.xyz`. The count of replicas away from
the pinned basin is what says whether the search is still exploring.

Note that i-PI writes trajectories in **atomic units**, with the cell on a
`# CELL(abcABC):` comment rather than an extxyz `Lattice=`. `ase.io.read` takes
neither: it reads bohr as Angstrom and finds no cell. Use i-PI's own reader for
these files. `reference/found_fd3m.xyz` is written out properly for ASE.

## Parallelism

Clients evaluate beads concurrently, but the gain saturates on a system this
small:

| clients | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| speedup | 1.00x | 1.17x | 1.90x | 2.51x | 2.59x |

A NEP evaluation on a few tens of atoms takes milliseconds, so i-PI's own
per-step cost sets a floor that no number of clients removes. Eight workers get
essentially all of the available speedup.

## Requirements

`<normal_modes propagator='bab'>` is **required**: freezing a bead is not a
diagonal constraint in the ring-polymer normal-mode basis, so the exact and
Cayley propagators refuse it. `bab` also cannot be combined with a barostat, so
RPQA runs are NVT.
