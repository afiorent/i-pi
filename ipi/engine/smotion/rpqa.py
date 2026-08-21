"""Ring-polymer quantum annealing, driven as a super-motion.

Periodically interrupts the annealing dynamics to relax every bead, pin the
replica sitting in the best minimum, and let the dynamics carry on around it.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import time

import numpy as np

from ipi.engine.smotion import Smotion
from ipi.engine.motion import Dynamics, GeopMotion
from ipi.engine.barostats import Barostat
from ipi.engine.thermostats import ThermoSVR
from ipi.engine.normalmodes import active_beads_mask
from ipi.utils import io
from ipi.utils.depend import dstrip
from ipi.utils.messages import verbosity, info, warning
from ipi.utils.units import unit_to_user

__all__ = ["RPQA"]


def _find_dynamics(motion):
    """Digs the Dynamics motion out of a (possibly nested) MultiMotion."""

    if isinstance(motion, Dynamics):
        return motion
    for m in getattr(motion, "mlist", []):
        found = _find_dynamics(m)
        if found is not None:
            return found
    return None


class RPQA(Smotion):
    """Ring-polymer quantum annealing.

    Runs the RPQA loop that used to be a shell script driving i-PI four times
    per iteration: relax every bead, pick the one in the deepest minimum, move
    that one replica to its minimum while the rest of the ring polymer stays
    where the dynamics left it, freeze it, and hand control back to the
    annealing dynamics for another interval.

    The dynamics itself is *not* run from here -- it is the system's own
    <motion>, stepped by the main loop as usual, and typically a <multi> of
    <dynamics> and <qkin_ramp>. This class only owns the relaxation and the
    pinning, so the annealing schedule stays where a reader expects it.

    Attributes:
        pinning_interval: How many steps of dynamics between two pinning events.
        start_step: First step at which a pinning event happens; the run before
            it is the equilibration ("delocalization") stage.
        max_relax_steps: Safety cap on the optimizer iterations per event.
        pinfile: Where the pinning history is logged.
        inherent_data: Whether to save the inherent structures and their energies.
        geop: One owned GeopMotion per system, used for the relaxations.
        pinned_bead: The currently pinned replica of each system, -1 if none.
    """

    def __init__(
        self,
        pinning_interval=1000,
        start_step=-1,
        max_relax_steps=5000,
        pinfile="rpqa_pin",
        inherent_data=False,
        optimizer=None,
        pinned_bead=None,
    ):
        """Initialises RPQA.

        Args:
           pinning_interval: Steps of dynamics between pinning events.
           start_step: Step of the first pinning event. Negative means "one
              pinning_interval", i.e. equilibrate for as long as one interval.
           max_relax_steps: Maximum optimizer iterations per pinning event.
           pinfile: Name of the pinning log.
           inherent_data: If True, save the relaxed configuration of every bead
              at each pinning event, along with its potential energy.
           optimizer: A dict of GeopMotion options, as produced by InputGeop.
           pinned_bead: Restart state -- the pinned replica of each system.
        """

        super(RPQA, self).__init__()

        self.pinning_interval = int(pinning_interval)
        if self.pinning_interval < 1:
            raise ValueError("RPQA pinning_interval must be a positive number of steps")

        if start_step < 0:
            self.start_step = self.pinning_interval
        else:
            self.start_step = int(start_step)

        self.max_relax_steps = int(max_relax_steps)
        self.pinfile = pinfile
        self.inherent_data = bool(inherent_data)

        if optimizer is None:
            optimizer = {}
        self.optimizer_options = dict(optimizer)
        # The relaxation must never terminate the simulation when it converges,
        # and its convergence must not be routed through motion.finished: the
        # main loop skips smotion.step() entirely on the step where any motion
        # reports finished, so we would be shut out of our own workflow.
        self.optimizer_options["exit_on_convergence"] = False

        if pinned_bead is None:
            self.pinned_bead = np.zeros(0, int)
        else:
            self.pinned_bead = np.asarray(pinned_bead, int).copy()

        self.mode = "rpqa"

    def bind(self, syslist, prng, omaker):
        super(RPQA, self).bind(syslist, prng, omaker)

        if len(self.pinned_bead) == 0:
            self.pinned_bead = -np.ones(len(self.syslist), int)
        elif len(self.pinned_bead) != len(self.syslist):
            raise ValueError(
                "Number of pinned beads does not match the number of systems"
            )

        self.dynamics = []
        self.geop = []
        for s in self.syslist:
            if s.beads.nbeads == 1:
                raise ValueError(
                    "RPQA needs a ring polymer with more than one bead; "
                    "there is nothing to pin in a classical simulation."
                )
            # Checked here rather than at the first pinning event, which could
            # be thousands of steps into a run.
            if s.nm.propagator != "bab":
                raise ValueError(
                    "RPQA freezes a bead, which only the Cartesian "
                    "free-ring-polymer propagator supports. Set "
                    "<normal_modes propagator='bab'>."
                )

            dyn = _find_dynamics(s.motion)
            if dyn is None:
                raise ValueError(
                    "RPQA drives the annealing through the system's own "
                    "dynamics, but no <motion mode='dynamics'> was found."
                )
            # The number of frozen degrees of freedom changes at the first
            # pinning event, but fixdof was handed to the thermostat and the
            # barostat once, at bind, when nothing was pinned yet. Langevin
            # ignores fixdof, so the common case is safe; the consumers that do
            # not are refused rather than left quietly wrong.
            if type(dyn.barostat) is not Barostat:
                raise ValueError(
                    "RPQA does not support a barostat: the barostat's dof count "
                    "is fixed at bind time and cannot follow the pinning."
                )
            if isinstance(dyn.thermostat, ThermoSVR):
                raise ValueError(
                    "RPQA does not support the SVR thermostat: its ndof is "
                    "fixed at bind time and cannot follow the pinning. Use "
                    "<thermostat mode='langevin'>."
                )

            self.dynamics.append(dyn)

            geop = GeopMotion(
                fixcom=dyn.fixcom,
                fixatoms_dof=dyn.fixatoms_dof,
                **self.optimizer_options,
            )
            # Binding a second motion to the system is safe: the optimizer
            # works on its own clones of the beads and forces, and only writes
            # back through beads.q / transfer_forces.
            geop.bind(s.ensemble, s.beads, s.nm, s.cell, s.forces, self.prng, omaker)
            self.geop.append(geop)

        # restores the constraint of a run resumed from a checkpoint
        for isys, k in enumerate(self.pinned_bead):
            if k >= 0:
                self._pin(isys, k)

        self.pf = self.output_maker.get_output(self.pinfile)
        if self.output_maker.f_start:  # a fresh run, not a resumed one
            self.pf.write(
                "#     step  sys  bead   potential/eV        spread/eV"
                "   relaxsteps      lambdaqkin\n"
            )
            self.pf.force_flush()

        self._bind_inherent()

    def _bind_inherent(self):
        """Opens the inherent-structure outputs, if they were asked for.

        One xyz per bead, appended once per pinning event, plus a table of the
        relaxed energies. get_output picks "w" or "a" from the restart state,
        so a resumed run keeps adding to the existing files.
        """

        self.ifile = None
        self.ipos = []
        if not self.inherent_data:
            return

        # bead index is zero-padded the same way trajectory outputs do it
        tagged = len(self.syslist) > 1
        for isys, s in enumerate(self.syslist):
            digits = int(1 + np.floor(np.log(s.beads.nbeads) / np.log(10)))
            tag = ("_s%d" % isys) if tagged else ""
            self.ipos.append(
                [
                    self.output_maker.get_output(
                        "rpqa_inherent%s.pos_%0*d.xyz" % (tag, digits, b)
                    )
                    for b in range(s.beads.nbeads)
                ]
            )

        self.ifile = self.output_maker.get_output("rpqa_inherent.out")
        if self.output_maker.f_start:
            self.ifile.write(
                "# RPQA inherent structures: the relaxed configuration of every bead at\n"
                "# each pinning event. Positions are in the companion .pos_*.xyz files,\n"
                "# one frame per event.\n"
                "#     step  sys  pinned   potential/eV of bead 0, 1, ... in order\n"
            )
            self.ifile.force_flush()

    def _write_inherent(self, isys, step, pots, k):
        """Saves the relaxed structures and energies of one pinning event.

        Must be called while the beads still hold the relaxed positions, i.e.
        before the pre-relax configuration is put back.
        """

        s = self.syslist[isys]
        pots_ev = [unit_to_user("energy", "electronvolt", p) for p in pots]

        for b in range(s.beads.nbeads):
            io.print_file(
                "xyz",
                s.beads[b],
                s.cell,
                self.ipos[isys][b],
                # the trailing space matters: print_file appends the key and
                # units to whatever it is given
                title=(
                    "RPQA inherent  Step:  %10d  Bead:   %5d  Potential: %15.8e eV%s "
                    % (step, b, pots_ev[b], "  PINNED" if b == k else "")
                ),
                key="positions",
                dimension="length",
            )
            self.ipos[isys][b].force_flush()

        self.ifile.write(
            "% 10d % 5d % 7d" % (step, isys, k)
            + "".join(" %15.8e" % p for p in pots_ev)
            + "\n"
        )
        self.ifile.force_flush()

    def _pin(self, isys, k):
        """Freezes bead k of system isys, releasing whichever was frozen before."""

        s = self.syslist[isys]
        dyn = self.dynamics[isys]

        dyn.fixbeads = np.asarray([k], int)
        # The thermostat's fixdof was computed once at bind from get_fixdof(),
        # so the *number* of frozen degrees of freedom must not change here --
        # only which bead they belong to. Pinning several beads at once would
        # need that count refreshed too.
        assert len(dyn.fixbeads) == 1, "RPQA pins exactly one bead at a time"

        mask = active_beads_mask(
            s.beads.nbeads, s.beads.natoms, dyn.fixatoms_dof, dyn.fixbeads
        )
        s.nm.activebeads_mask = mask
        dyn.integrator.fixbeads = dyn.fixbeads
        dyn.integrator.activebeads_mask = mask

        # Properties reads fixbeads off system.motion, which for the usual
        # <motion mode='multi'> input is the MultiMotion and not the Dynamics
        # we just pinned. Without this the reported temperature is low by
        # (frozen dof)/(total dof), since get_temp cannot compensate for a
        # constraint it cannot see.
        if s.motion is not dyn:
            s.motion.fixbeads = dyn.fixbeads

        self.pinned_bead[isys] = k

    def _relax(self, isys):
        """Minimizes every bead, returning the number of iterations used."""

        geop = self.geop[isys]
        geop.reset()

        for i in range(self.max_relax_steps):
            # The counter is local and starts at zero because bfgs and lbfgs
            # key their first search direction on step == 0; a global step
            # would silently skip that initialisation.
            geop.step(i)
            if geop.optimizer.converged:
                return i + 1

        warning(
            "@RPQA: relaxation hit max_relax_steps=%d without converging"
            % self.max_relax_steps,
            verbosity.low,
        )
        return self.max_relax_steps

    def step(self, step=None):
        """Relaxes, picks the best replica and pins it."""

        if step is None or step < self.start_step:
            return
        if (step - self.start_step) % self.pinning_interval != 0:
            return

        t_start = time.time()

        for isys, s in enumerate(self.syslist):
            # The relaxation moves the beads, so the pre-relax configuration
            # has to be kept: everything except the winning replica is put back.
            q0 = dstrip(s.beads.q).copy()
            econs0 = s.ensemble.econs

            nrelax = self._relax(isys)

            pots = dstrip(s.forces.pots).copy()
            k = int(np.argmin(pots))

            if self.inherent_data:
                # while beads.q is still the relaxed configuration
                self._write_inherent(isys, step, pots, k)

            # Only the best bead keeps its minimum. The others go back to where
            # the dynamics left them, so the ring polymer stays delocalised;
            # the relaxation serves to find and place the replica to pin.
            q_new = q0
            q_new[k] = dstrip(s.beads.q)[k]
            s.beads.q[:] = q_new

            self._pin(isys, k)

            # The pinning is not a physical move, so its energy goes to the
            # ensemble's "external" reservoir and the conserved quantity stays
            # continuous across the event.
            s.ensemble.eens += econs0 - s.ensemble.econs

            info(
                "@RPQA: step %d, system %d: pinned bead %d at %f, "
                "spread %e, %d relaxation steps"
                % (step, isys, k, pots[k], pots.max() - pots.min(), nrelax),
                verbosity.medium,
            )

            self.pf.write(
                "% 10d % 5d % 5d %15.8e %15.8e % 8d %15.8e\n"
                % (
                    step,
                    isys,
                    k,
                    unit_to_user("energy", "electronvolt", pots[k]),
                    unit_to_user("energy", "electronvolt", pots.max() - pots.min()),
                    nrelax,
                    s.ensemble.lambdaqkin,
                )
            )
        self.pf.force_flush()

        info(
            "# RPQA pinning evaluated in %f sec." % (time.time() - t_start),
            verbosity.debug,
        )
