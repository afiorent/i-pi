# Ring-polymer quantum annealing of silicon

RPQA applied to a periodic silicon cell, with forces from a NEP machine-learning
potential through an ASE socket client. The method is described in

> A. Fiorentino and N. Marzari, *Quantum annealing for materials*,
> [arXiv:2606.03405](https://arxiv.org/abs/2606.03405) (2026)

```
./run.sh
```

About 13 minutes with 8 clients. The driver needs **pyNEP**
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

| | |
|---|---|
| system | 16 Si, periodic, cell from the structure file |
| beads | 32 |
| temperature | 50 K |
| `lambdaqkin` | 1200 → 300, `sqrtscale` |
| timestep / `tau` | 2 fs / 20 fs |
| `nmts` | 20 |
| delocalization | 2500 steps (5 ps) |
| `pinning_interval` | 2500 steps (5 ps), 4 events |
| relaxation | `cg_rp`, 0.01 eV/A, 1 meV, 0.05 A, capped at 1000 steps |

## Two things that differ from the LJ examples

**`lambdaqkin` starts very high.** Nuclear quantum effects scale as
`lambdaqkin * hbar^2 / M`, and silicon is heavy: a value that delocalizes a
light nucleus does nothing here. 1200 is the right order for Si; see the paper
for how to choose the range.

**The anneal is partial, on purpose.** Production runs decrement lambda slowly —
in the reference setup `sqrt(lambdaqkin)` drops by about 0.9 per 2500-step
interval, so reaching `lambdaqkin` = 1, the physical value, takes on the order of
forty intervals and 10^5 steps. This example covers four intervals and stops at
300. It demonstrates the machinery, not a converged anneal.

That distinction matters, because the ramp rate is not a free parameter. Forcing
`lambdaqkin` from 1200 to 1 in 200 steps, as a first attempt at this example did,
stiffens the ring-polymer springs by a factor of 35 almost instantly and does
enormous work on the system: the temperature spiked to 5000 K. At the rate used
here the run holds 50.3 K on average against a 50 K target.

## Why `cg_rp` and not `lbfgs`

At high `lambdaqkin` the replicas are spread over several eV and reach their
minima at very different rates. `cg_rp` gives each bead its own convergence test
and step length; every other optimizer shares one convergence test across all
beads, so the whole relaxation waits for the worst replica.

The difference is not marginal here. With `lbfgs` this system ran 2000
iterations with the energy frozen and the maximum force stuck at 1.7e-2 Ha/bohr,
never converging. With `cg_rp` the four relaxations took 194, 173, 519 and 165
steps.

## Reading the output

`rpqa.rpqa_pin` — one line per pinning event: the replica pinned, its energy,
the spread across replicas, the relaxation steps used, and the current
`lambdaqkin`. A typical run:

```
     step  sys  bead   potential/eV        spread/eV   relaxsteps      lambdaqkin
     2500    0    13  -7.40146356e+01   6.96430744e+00      194     9.18645003e+02
     5000    0    21  -7.40147185e+01   6.81873494e+00      173     6.74910003e+02
     7500    0     4  -7.40147200e+01   6.06113430e+00      519     4.68675003e+02
    10000    0     4  -7.40147219e+01   7.43356298e+00      165     3.00000000e+02
```

The pinned replica changes — 13, then 21, then 4 — which is the point: pinning
is not a commitment, and a replica that later quenches deeper takes over.

`rpqa.rpqa_inherent.out` — with `<inherent_data>True</inherent_data>`, the
quenched energy of *every* bead at every event, one column per bead, in eV, with
the structures in the companion `rpqa.rpqa_inherent.pos_*.xyz`. This is what
tells you whether the search is still exploring: replicas outside the pinned
basin are what can still find something better.

## Parallelism

Clients evaluate beads concurrently, but the gain saturates quickly on a system
this small:

| clients | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| speedup | 1.00x | 1.17x | 1.90x | 2.51x | 2.59x |

A NEP evaluation on 16 atoms takes a few ms, so i-PI's own per-step cost — socket
traffic, the normal-mode transforms, the dependency graph — sets a floor around
50 ms/step that no number of clients removes. Eight workers get essentially all
of the available speedup. A larger cell or a more expensive potential shifts the
balance and scales further.

## Requirements

`<normal_modes propagator='bab'>` is **required**: freezing a bead is not a
diagonal constraint in the ring-polymer normal-mode basis, so the exact and
Cayley propagators refuse it. Note also that `bab` cannot be combined with a
barostat, so RPQA runs are NVT.
