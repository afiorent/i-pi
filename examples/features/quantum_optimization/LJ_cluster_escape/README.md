# Recovering from a bad first pin

A companion to [../LJ_cluster](../LJ_cluster), showing what happens when the
delocalization is too short: the first replica gets pinned in the wrong basin,
and the annealing recovers from it anyway.

The method is described in

> A. Fiorentino and N. Marzari, *Quantum annealing for materials*,
> [arXiv:2606.03405](https://arxiv.org/abs/2606.03405) (2026)

```
./run.sh
```

Two clients, about 15 minutes.

## The point

Pinning is not a commitment. When the smotion relaxes the beads at a pinning
event it pins whichever replica is deepest *at that moment*; if a different
replica later quenches into a deeper basin, the pin moves there. So a poor
first choice can be undone.

This run makes that visible by setting `start_step` to 200 steps -- 0.2 ps,
against the 5 ps of `../LJ_cluster`. The ring polymer has barely spread by the
first pinning event, so it pins a mediocre minimum:

| event | step | pinned bead | energy (eV) | spread (eV) | beads in the global basin |
|---|---|---|---|---|---|
| 1 |   200 |  2 | **-10.102270** | 0.474 | 0 / 16 |
| 2 |  5200 | 14 | **-10.237267** | 0.304 | 5 / 16 |
| 3 | 10200 | 14 | -10.237267 | 0.0005 | 16 / 16 |
| 4 | 15200 | 14 | -10.237267 | 0.300 | 14 / 16 |
| 5 | 20200 | 14 | -10.237267 | 0.00003 | 16 / 16 |
| 6 | 25200 | 14 | -10.237267 | 0.0003 | 16 / 16 |
| 7 | 30200 | 14 | -10.237267 | 0.001 | 14 / 16 |
| 8 | 35200 | 14 | -10.237267 | 0.00007 | 16 / 16 |

Bead 2 is pinned at event 1, 0.135 eV above the global minimum, with no replica
anywhere near the right basin. During the next interval five of them find it,
and the pin moves to bead 14. The spread is still 0.304 eV at that point: the
ring is delocalized and exploring, which is precisely what makes the recovery
possible.

What follows is *not* a clean collapse. The spread keeps jumping back up --
0.300 eV at event 4, when two replicas have wandered back out of the pinned
basin, and again at event 7 -- because `lambdaqkin` is still appreciable
(6.3 and 0.50 at those events). Only at the end, fully quenched, does the ring
sit still. Those excursions are the search still working: a replica outside the
pinned basin is a replica that could find something better.

The escape has to happen early, while `lambdaqkin` is large. Extra pinning
events cannot supply it later: once the anneal has taken the ring down to small
`lambdaqkin` it stops leaving the pinned basin at all, and further events only
re-relax the same minimum. `escape.rpqa_inherent.out` is where to look -- the
count of replicas outside the pinned basin, rather than the pinned energy
itself, is what says whether the search is still alive.

## These settings are not for production

They are chosen to make the recovery visible in a short run, not to converge
reliably. For real work use a **longer delocalization and longer
trajectories**, as in `../LJ_cluster`:

| | this example | ../LJ_cluster |
|---|---|---|
| delocalization | 200 steps (0.2 ps) | 5000 steps (5.0 ps) |
| interval | 5000 steps | 5000 steps |
| pinning events | 8 | 5 |

Recovery is a safety net, not a substitute for delocalizing properly. Over five
runs at these settings -- the five starting structures below, each with the seed
matching its index -- four pinned the global minimum immediately and only one
had to recover, so even 0.2 ps is usually enough for this system.

Recovery is not guaranteed, though. In a separate run at 0.5 ps, one trajectory
pinned -10.1867 eV at the first event and was still there after eight, with the
replicas having settled into that basin by event 3 and never leaving it again.
Whether a badly-started run recovers depends on a better basin turning up while
`lambdaqkin` is still large enough for the ring to reach it.

## The starting structures

Five different starting configurations ship here, `init_0.xyz` to
`init_4.xyz`. `input.xml` uses `init_4.xyz`, the one that produces the escape;
change that one line to try the others.

**None of them contains the global minimum.** Quenching each directly, without
any annealing, gives:

| file | quenched energy (eV) | above the global minimum |
|---|---|---|
| `init_0.xyz` | -10.080072 | +0.157 |
| `init_1.xyz` |  -9.942079 | +0.295 |
| `init_2.xyz` |  -9.746572 | +0.491 |
| `init_3.xyz` |  -10.085294 | +0.152 |
| `init_4.xyz` |  -9.780971 | +0.456 |

The global minimum is -10.2373 eV (-102.372663 $\epsilon$ in the
[Cambridge Cluster Database](https://doye.chem.ox.ac.uk/jon/structures/LJ.html),
with $\epsilon$ = 0.1 eV). So in every case the annealing has to find a basin
the starting structure is not in -- the result is not an artefact of a lucky
initial guess. How deep the starting quench is does not predict how the run
goes: `init_3.xyz` starts closest and `init_2.xyz` furthest, and neither is the
one that has to recover.

## Reproducibility

Each `input.xml` carries an explicit `<prng><seed>`, and the trajectory is
bit-reproducible with the two clients `run.sh` starts: which client evaluates
which bead does not matter, because the force is deterministic. Running this
example twice gives identical output, so the table above is what you should
see.
