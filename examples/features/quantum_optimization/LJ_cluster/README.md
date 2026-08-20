# Ring-polymer quantum annealing of an LJ$_{25}$ cluster

Finds the global minimum of a 25-atom Lennard-Jones cluster in a single i-PI
run, driven by the `rpqa` super-motion.

```
./run.sh
```

Two clients, about 12 minutes. The answer to look for is **-10.2373 eV**, the
LJ$_{25}$ global minimum (-102.372663 $\epsilon$ in the
[Cambridge Cluster Database](https://doye.chem.ox.ac.uk/jon/structures/LJ.html),
with $\epsilon$ = 0.00367493 Ha = 0.1 eV).

## What it does

The ring polymer is propagated at 50 K while `lambdaqkin` -- an effective
$\hbar^2$ that scales the quantum delocalization -- is ramped from 19 down to
0.01. Every `pinning_interval` steps the smotion:

1. relaxes every bead to its nearest minimum (an *inherent structure*),
2. pins the replica sitting in the deepest one, freezing it for good,
3. puts the other beads back where the dynamics left them and hands control
   back, so the ring polymer keeps exploring around the frozen replica.

**The delocalization stage is not separate.** It is the first `start_step`
steps, during which the smotion does nothing at all. Its length is the single
most important parameter here -- see below.

## Parameters

| | |
|---|---|
| beads / temperature | 16 / 50 K |
| `lambdaqkin` | 19 → 0.01, `sqrtscale` |
| timestep / thermostat `tau` | 1.0 fs / 20 fs |
| `start_step` (delocalization) | 5000 steps = 5.0 ps |
| `pinning_interval` | 5000 steps = 5.0 ps, 5 pinning events |
| `total_steps` | 30000 |
| ramp `total_steps` | **25000** — the last pinning event, not the end of the run |
| relaxation | `lbfgs` |

Two things are easy to get wrong:

- `<normal_modes propagator='bab'/>` is **required**. Freezing a bead is not a
  diagonal constraint in the ring-polymer normal-mode basis, so the exact and
  Cayley propagators refuse it.
- The ramp's `<total_steps>` is the step at which the anneal must be
  *finished*. Setting it to the length of the run leaves `lambdaqkin` at ~0.002
  when the run stops, with the ring still delocalized.

## Reading the output

`rpqa.rpqa_pin` — one line per pinning event: which bead was pinned, its
energy, the spread across beads, how many optimizer steps the relaxation took,
and the current `lambdaqkin`.

`rpqa.rpqa_inherent.out` — with `<inherent_data>True</inherent_data>`, the
quenched energy of *every* bead at every event, one column per bead (eV). The
companion `rpqa.rpqa_inherent.pos_*.xyz` hold the corresponding structures, one
frame per event.

That table is what makes a failure diagnosable. A run that misses the global
minimum does so because no replica ever quenched into the right basin -- which
is visible here, and not in the energies alone.

## How long does the delocalization need to be?

Long enough for at least one replica to reach the global-minimum basin. **One
is enough**: in testing, a run with exactly one bead in that basin at the first
pinning event pinned it and had all sixteen there by the second event. This
example is less marginal than that -- it typically gets several -- but the
margin is what the delocalization length buys, so it is worth understanding.

What matters is physical *time*, not step count:

| delocalization | timestep | time | outcome |
|---|---|---|---|
| 10000 steps | 0.5 fs | 5.0 ps | global minimum |
| **5000 steps** | **1.0 fs** | **5.0 ps** | global minimum (this example, half the cost) |
| 5000 steps | 0.5 fs | 2.5 ps | misses it on harder structures |

Halving the steps while doubling the timestep and `tau` keeps the result and
halves the cost. The springs stay very well resolved: they are integrated at
`dt/nmts/2` = 0.05 fs with the default `nmts`=10, so $\omega_{max}\Delta t
\approx 0.1$ even at the stiff end of the anneal.

Extra pinning events do not rescue a delocalization that was too short -- they
only re-relax around the replica already pinned.

## The starting structure

`init.xyz` is a quenched-random start whose own local minimum is
**-10.0801 eV**, 0.157 eV above the global one, so the annealing is doing the
work rather than the initial guess. Across ten such structures the starting
minima span -9.56 to -10.09 eV, and how deep the start is does not predict
whether the run succeeds.

## Why it is one run

RPQA was originally driven from a shell script that restarted i-PI for every
stage -- relax, then a Python script to pick the replica to pin, then the
annealing interval -- so a five-iteration run meant eleven i-PI startups and a
directory full of per-stage fragments to stitch back together.

Doing the same thing from inside the engine costs the same in force
evaluations, since it runs the same dynamics and the same relaxations. What it
saves is the repeated startup, and it keeps momenta across the pinning events
instead of restarting each stage from rest. One trajectory, one RESTART that
resumes anywhere.
