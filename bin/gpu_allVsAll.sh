#!/usr/bin/env bash



DIR=$1
EXT=$2
LEN=$3
DEV=$4

array=()
x=0


BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


if [ $# != 4 ]; then
    echo "***ERROR*** Use: $0 <genomesDirectory> <extension> <minLen> <device>"
    exit -1
fi

for elem in $(ls -d $DIR/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
    array[$x]=$elem
    x=`expr $x + 1`
done

total=$(echo "$x" | awk '{print ($0 * ($0-1)/2) }')
curr=0

for ((i=0 ; i < ${#array[@]} ; i++))
do
    for ((j=i ; j < ${#array[@]} ; j++))
    do
        if [ $i != $j ]; then
                seqX=${array[$i]}
                seqY=${array[$j]}
                echo "[GPU $DEV]: ($curr/$total) ${seqX}-${seqY}"
                { time $BINDIR/gpu_cuda_workflow -query $DIR/${seqX}.$EXT -ref $DIR/${seqY}.$EXT -dev $DEV -len $LEN >/dev/null 2> errors; } 2>&1 | grep real | awk '{print $2}'
                cat errors

                curr=`expr $curr + 1`
            
        fi
    done
done

