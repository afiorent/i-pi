"""Deals with creating the RPQA smotion class.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.

Classes:
   InputRPQA: Deals with creating the RPQA object from a file, and
      writing the checkpoints.
"""

import numpy as np
from ipi.utils.inputvalue import *
from ipi.inputs.motion.geop import InputGeop

__all__ = ["InputRPQA"]


class InputRPQA(InputDictionary):
    """Ring-polymer quantum annealing options.

    Contains the relaxation options and the schedule on which replicas are
    pinned. The annealing dynamics itself is not configured here: it is the
    system's own <motion>, typically a <multi> of <dynamics> and <qkin_ramp>.
    """

    fields = {
        "pinning_interval": (
            InputValue,
            {
                "dtype": int,
                "default": 1000,
                "help": "How many steps of dynamics to run between two pinning events.",
            },
        ),
        "start_step": (
            InputValue,
            {
                "dtype": int,
                "default": -1,
                "help": "Step at which the first pinning event happens; the run up to it "
                "is the equilibration (delocalization) stage. A negative value means one "
                "pinning_interval.",
            },
        ),
        "max_relax_steps": (
            InputValue,
            {
                "dtype": int,
                "default": 5000,
                "help": "Safety cap on the number of optimizer iterations spent on a "
                "single relaxation. Reaching it produces a warning.",
            },
        ),
        "pinfile": (
            InputValue,
            {
                "dtype": str,
                "default": "rpqa_pin",
                "help": "File to keep track of which replica is pinned when.",
            },
        ),
        "pinned_bead": (
            InputArray,
            {
                "dtype": int,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The replica currently pinned in each system, -1 if none. "
                "Restart state, normally not set by hand.",
            },
        ),
        "optimizer": (
            InputGeop,
            {
                "default": {},
                "help": "Options for the relaxation performed at each pinning event. "
                "exit_on_convergence is forced off: convergence hands control back to "
                "RPQA rather than ending the simulation.",
            },
        ),
    }

    default_help = "Ring-polymer quantum annealing"
    default_label = "RPQA"

    def store(self, rpqa):
        if rpqa == {}:
            return
        self.pinning_interval.store(rpqa.pinning_interval)
        self.start_step.store(rpqa.start_step)
        self.max_relax_steps.store(rpqa.max_relax_steps)
        self.pinfile.store(rpqa.pinfile)
        self.pinned_bead.store(rpqa.pinned_bead)
        # the owned optimizers are per system and identically configured, so
        # the first one carries the options; before bind there is none yet
        if len(getattr(rpqa, "geop", [])) > 0:
            self.optimizer.store(rpqa.geop[0])

    def fetch(self):
        rv = super(InputRPQA, self).fetch()
        return rv
