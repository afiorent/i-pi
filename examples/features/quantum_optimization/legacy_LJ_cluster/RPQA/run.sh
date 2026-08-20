source ~/useful_repos/fork_ipi/env.sh

echo 'Start' $(date)
nbeads=16
Niteration=5
Nsplit=$((Niteration -1))
nworker=2
RP_interval=2000
IPI_ADDRESS=localhost
sigma=5.270446
epsilon=0.00367493
cutoff=36.893122 ## 7*sigma 
IPI_DRIVER="i-pi-driver -m lj -o $sigma,$epsilon,$cutoff -u -a $IPI_ADDRESS" #sigma,epsilon,cutoff
rm /tmp/ipi_${IPI_ADDRESS}

## We assume i-pi (and the driver code) are in the path, otherwise
## you have to set this environment variable
IPI_PATH= 


## Driver command

#export PATH=$PATH:${IPI_PATH}/bin



#### Equilibration
## Launches i-PI
i-pi delo.xml &> log.ipi-delo &

## Gives a few seconds to allow the server to open the Unix socket
## For *very* complicated simulations you may need to increase a bit
sleep 5; 

## Launches the driver code
for nbead in `seq 1 $nworker`; do
    $IPI_DRIVER &> log.driver.$nbead &
done

## Waits for all jobs to be finished
wait

echo "Done equilibration"

for i in $(seq 0 $Nsplit)
do

echo "Local minimization step minimization"
## Launches i-PI
cp relax_template.xml min_$i.xml
sed -i -e "s/'min'/'min$i'/g" min_$i.xml
sed -i -e "s/nbeads='32'/nbeads='$nbeads'/g" min_$i.xml
IPI_INPUT=min_${i}.xml
# increase relaxation accuracy for last step
if [[ "$i" == "$Nsplit" ]]; then
    # Run sed substitution
    sed -i -e 's/<force> 1e-3/<force> 1e-5/g' "$IPI_INPUT"
    sed -i -e 's/<position> 1e-3/<position> 1e-5/g' "$IPI_INPUT"
    echo "Reducing relaxation threshold"
fi
rm /tmp/ipi_${IPI_ADDRESS}
## Launches i-PI
i-pi ${IPI_INPUT} &> log.min_$i &

## Gives a few seconds to allow the server to open the Unix socket
## For *very* complicated simulations you may need to increase a bit
sleep 5; 

## Launches the driver code
for nbead in `seq 1 $nworker`; do
    $IPI_DRIVER &> log.driver.$nbead &
done

## Waits for all jobs to be finished
wait

echo "Done minimization"
echo 'update RPQA input'
python3 update_RPQA.py --iteration $i --hbar2-0 16 --hbar2-1 0.0001 --N-iteration ${Niteration} --simu-steps $RP_interval --nbeads $nbeads >> update.out 
## Input file
IPI_INPUT=input_${i}.xml
rm /tmp/ipi_${IPI_ADDRESS}
## Launches i-PI
i-pi ${IPI_INPUT} &> log.ipi_$i &

## Gives a few seconds to allow the server to open the Unix socket
## For *very* complicated simulations you may need to increase a bit
sleep 5; 

## Launches the driver code
for nbead in `seq 1 $nworker`; do
    $IPI_DRIVER &> log.driver.$nbead &
done

## Waits for all jobs to be finished
wait
echo "Done QA part $i"


done

echo 'joining traj and cleaning'
python3 join_traj_and_output.py --niterations $Niteration --nbeads $nbeads --out-json complete_simulation.json
rm min*xyz
echo 'End' $(date)
