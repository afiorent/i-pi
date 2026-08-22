# Quantum annealing

Finding minima by annealing the quantum nuclear density: the system is first
delocalized over a large effective $\hbar$, so that it samples several basins at
once, and then quenched. Depending on where the anneal ends, the same machinery
is either a global optimizer or a route to nuclear quantum effects directly.

The method is described in

> A. Fiorentino and N. Marzari, *Quantum annealing for materials*,
> [arXiv:2606.03405](https://arxiv.org/abs/2606.03405) (2026)

## [RPQA/](RPQA) — replica-pinned quantum annealing

The ring polymer is annealed while `lambdaqkin`, a prefactor on the quantum
kinetic energy equivalent to scaling $\hbar^2$, is ramped down. Periodically the
`rpqa` super-motion relaxes every replica, pins the one sitting in the deepest
minimum, and lets the dynamics carry on around it.

- **[LJ_cluster/](RPQA/LJ_cluster)** — a 25-atom Lennard-Jones cluster, and the
  place to start: settings that converge reliably, in about 12 minutes.
- **[LJ_cluster_escape/](RPQA/LJ_cluster_escape)** — the same system with a
  deliberately too-short delocalization, showing that a replica pinned in the
  wrong basin can still be displaced later by a better one. Explains what the
  pinning does; not settings to copy.
- **[Si_reconstruction_ASE/](RPQA/Si_reconstruction_ASE)** — 32 silicon atoms in
  a periodic cell with a NEP machine-learning potential driven through ASE.
  Condensed phase, and the anneal ends at the physical `lambdaqkin` rather than
  quenching to a classical minimum.

## [QA/](QA) — quantum annealing without pinning

Reserved for the asymmetric double-well example.

## Reference data

Reference structures and energies for LJ clusters, in LJ units, are in the
Cambridge Cluster Database: https://doye.chem.ox.ac.uk/jon/structures/LJ.html
