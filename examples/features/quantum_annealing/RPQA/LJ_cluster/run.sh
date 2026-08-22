#!/bin/bash
# RPQA on a 25-atom LJ cluster: one i-PI run, driven by the 'rpqa' smotion.
#
# i-PI is started once and the relax/pick/pin logic runs inside the engine, so
# there is a single continuous trajectory and a single RESTART covering the
# whole run.

source ../../../../../env.sh

nworker=2
IPI_ADDRESS=rpqa_lj
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
echo "Which replica was pinned when:"
cat rpqa.rpqa_pin
echo
echo "The last column of rpqa.rpqa_inherent.out holds the quenched energy of"
echo "every bead at each pinning event; the global minimum is -10.2373 eV."
