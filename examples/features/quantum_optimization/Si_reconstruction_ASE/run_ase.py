"""ASE socket client driving a NEP potential for the Si RPQA example.

Only the species and the cell of the atoms object matter here: i-PI sends the
positions of each bead over the socket and this just evaluates the potential.

Needs `pynep` and `ase`; see README.md.
"""

import os

from ase.calculators.socketio import SocketClient
from ase.io import read
from pynep.calculate import NEP

HERE = os.path.dirname(os.path.abspath(__file__))

# paths are resolved relative to this file so the driver can be started from
# anywhere, and so the example does not carry one machine's scratch paths
POTENTIAL = os.path.join(HERE, "potential", "Si_Fan_GAP.txt")
STRUCTURE = os.path.join(HERE, "structures", "N16", "nbeads_32", "optimal_beads_0.xyz")

# the structure file holds one frame per bead; any of them will do to set up
# the calculator, since i-PI overwrites the positions on every call
atoms = read(STRUCTURE, index=0)
atoms.calc = NEP(POTENTIAL)

# must match <address> in input.xml
client = SocketClient(unixsocket="rpqa_si")
client.run(atoms)
