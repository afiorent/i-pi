# Quantum optimization

Finding global minima by annealing a ring polymer: the path integral is
delocalized over a large effective $\hbar$, then quenched, so that the replicas
explore several basins at once and settle into a deep one.

- **[LJ_cluster/](LJ_cluster)** — ring-polymer quantum annealing (RPQA) of a
  25-atom Lennard-Jones cluster, as a single i-PI run driven by the `rpqa`
  super-motion. Start here.
- **[legacy_LJ_cluster/](legacy_LJ_cluster)** — the earlier scripted versions of
  the same ideas, where i-PI is restarted for each stage and a Python script
  picks the replica to pin in between: `QA/` (quantum annealing) and `RPQA/`
  (replica-pinned quantum annealing). Kept for reference.

Reference structures and energies for LJ clusters, in LJ units, are in the
Cambridge Cluster Database: https://doye.chem.ox.ac.uk/jon/structures/LJ.html
