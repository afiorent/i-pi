source ~/useful_repos/fork_ipi/env.sh

echo 'Start' $(date)
nworker=2
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

## Input file
IPI_INPUT=input.xml
rm /tmp/ipi_${IPI_ADDRESS}
## Launches i-PI
i-pi ${IPI_INPUT} &> log.ipi &

## Gives a few seconds to allow the server to open the Unix socket
## For *very* complicated simulations you may need to increase a bit
sleep 5; 

## Launches the driver code
for nbead in `seq 1 $nworker`; do
    $IPI_DRIVER &> log.driver.$nbead &
done

## Waits for all jobs to be finished
wait
echo "Done QA part"

## Input file
IPI_INPUT=relax.xml
rm /tmp/ipi_${IPI_ADDRESS}
## Launches i-PI
i-pi ${IPI_INPUT} &> log.min &

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




echo 'End' $(date)
