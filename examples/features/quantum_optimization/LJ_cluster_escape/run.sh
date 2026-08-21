#!/bin/bash
# RPQA on a 25-atom LJ cluster: one i-PI run, driven by the 'rpqa' smotion.
#
# A deliberately under-delocalized RPQA run; see README.md.
# The first pinning event lands on a poor minimum and the run recovers from it.

source ../../../../env.sh

nworker=2
IPI_ADDRESS=rpqa_escape
sigma=5.270446
epsilon=0.00367493
cutoff=36.893122     # 7*sigma
IPI_DRIVER="i-pi-driver -m lj -o $sigma,$epsilon,$cutoff -u -a $IPI_ADDRESS"

rm -f /tmp/ipi_${IPI_ADDRESS}

echo "Start $(date)"
i-pi input.xml &> log.ipi &

# give the server a few seconds to open the socket; for heavier setups
# this may need increasing
sleep 5

for w in $(seq 1 $nworker); do
    $IPI_DRIVER &> log.driver.$w &
done

wait
echo "End $(date)"

echo
echo "Which replica was pinned when (watch the energy jump at event 2):"
cat escape.rpqa_pin
echo
echo "escape.rpqa_inherent.out holds the quenched energy of every bead at each"
echo "pinning event; the global minimum is -10.2373 eV."
