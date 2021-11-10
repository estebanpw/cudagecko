#!/bin/bash

if [ $# != 6 ]; then
   echo " ==== ERROR === "
   echo ""
   echo "   usage:  $0 <query> <name> <splits> <ram> <device id> <min len>"
   echo ""
   exit -1
fi

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
QUERY=$1
NAME=$2
SPLITS=$3
RAM=$4
DEV=$5
LEN=$6

# Tasks that run simultaneously go here
echo "Launching tasks"

for((i=0; i<$SPLITS; i++)); do

        $BINDIR/gpu_cuda_workflow -query $QUERY -ref ${NAME}${i}.fasta -dev $DEV -len $LEN -ram $RAM &

done


for job in `jobs -p`
do
    wait $job
done

echo "Completed all MPS executions"
