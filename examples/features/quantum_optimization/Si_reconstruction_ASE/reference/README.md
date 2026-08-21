# Reference output

Output of the run that `../input.xml` produces, kept so the example's claims can
be checked without re-running it. Seed 31415, 8 clients, 18.6 min wall.

| file | |
|---|---|
| `input.xml` | the exact configuration that produced this |
| `rpqa.rpqa_pin` | one line per pinning event |
| `rpqa.rpqa_inherent.out` | quenched energy of all 32 replicas at each event, eV |
| `rpqa.out` | step, time, conserved, temperature, kinetic_cv, potential, volume, lambdaqkin |
| `rpqa.beadpot` | per-replica potential along the trajectory |
| `inherent_pinned_bead30.xyz` | inherent structures of the finally-pinned replica, one frame per event |
| `found_fd3m.xyz` | the crystal it ends on, in Angstrom with the cell, ready for ASE |

## What the run does

Nine pinning events. The first four sit at +0.23 eV/atom, then at
`lambdaqkin` = 192 one replica finds the diamond structure, the pin moves to it,
and the ring collapses onto it over the next two events:

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

`found_fd3m.xyz` has coordination exactly 4.000 and a nearest-neighbour distance
of 2.3714 +- 0.0000 A, matching ideal diamond; its energy is -148.031745 eV
against -148.031823 eV for the relaxed Fd-3m reference in the same cell, a
difference of 2.4 ueV/atom.

Mean temperature over the run is 50.1 K against the 50 K target, and no
relaxation hit `max_relax_steps`.

## Note on reading the xyz files

i-PI writes trajectories in **atomic units**, with the cell on a
`# CELL(abcABC):` comment line rather than an extxyz `Lattice=`. ASE's reader
takes neither: `ase.io.read` on `inherent_pinned_bead30.xyz` returns positions
interpreted as Angstrom (they are bohr) and no cell at all. Use i-PI's own
reader, or `found_fd3m.xyz`, which is written out properly for ASE.

## Caveat

This is one run, one starting structure, one seed. Event 5 is a stochastic
escape; how reliably it happens has not been measured here.
