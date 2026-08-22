#!/bin/bash
# RPQA on amorphous silicon: one i-PI run driven by the 'rpqa' smotion, with
# forces from a NEP potential through ASE socket clients.
#
# The driver needs `pynep` and `ase` in the environment. If they are not in the
# python on your PATH, point PYTHON at one that has them, e.g.
#   PYTHON=~/miniforge3/envs/nep_env/bin/python ./run.sh

source ../../../../../env.sh

# Clients evaluate beads concurrently, but the gain saturates: measured on this
# example, 1/2/4/8/16 clients give 1.00/1.17/1.90/2.51/2.59x. A NEP evaluation
# on 16 atoms is only a few ms, so i-PI's own per-step cost is the floor and
# more than 8 workers buys nothing. A larger cell or a costlier potential
# would scale further.
nworker=${NWORKER:-8}
PYTHON=${PYTHON:-python}
IPI_ADDRESS=rpqa_si

if ! $PYTHON -c "import pynep, ase" 2>/dev/null; then
    echo "ERROR: the driver needs pynep and ase, and '$PYTHON' has neither." >&2
    echo "       Install pyNEP (https://github.com/bigd4/PyNEP) or set PYTHON" >&2
    echo "       to an interpreter that has it, e.g. a conda/mamba env." >&2
    exit 1
fi

rm -f /tmp/ipi_${IPI_ADDRESS}

echo "Start $(date)"
i-pi input.xml &> log.ipi &

# give the server a few seconds to open the socket; for heavier setups
# this may need increasing
sleep 5

for w in $(seq 1 $nworker); do
    $PYTHON run_ase.py &> log.driver.$w &
done

wait
echo "End $(date)"

echo
echo "Which replica was pinned when:"
cat rpqa.rpqa_pin
echo
echo "rpqa.rpqa_inherent.out holds the quenched energy of every bead at each"
echo "pinning event, one column per bead, in eV."
