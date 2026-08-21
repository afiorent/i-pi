# Quantum optimization

Finding global minima by annealing a ring polymer: the path integral is
delocalized over a large effective $\hbar$, then quenched, so that the replicas
explore several basins at once and settle into a deep one.

- **[LJ_cluster/](LJ_cluster)** — ring-polymer quantum annealing (RPQA) of a
  25-atom Lennard-Jones cluster, as a single i-PI run driven by the `rpqa`
  super-motion. Start here: it uses settings that converge reliably.
- **[LJ_cluster_escape/](LJ_cluster_escape)** — the same system with a
  deliberately too-short delocalization, showing that a replica pinned in the
  wrong basin can still be displaced later by a better one. Useful for
  understanding what the pinning does; not settings to copy.

Reference structures and energies for LJ clusters, in LJ units, are in the
Cambridge Cluster Database: https://doye.chem.ox.ac.uk/jon/structures/LJ.html
