#!/usr/bin/env bash
DIR=$1
DIR2=$2
LEN=$3
DEV=$4


if [ $# != 4 ]; then
    echo "***ERROR*** Use: $0 <genomesDirectory1> <genomesDirectory2> <minLen> <device>"
    exit -1
fi

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



array=()
x=0
array2=()

for elem in $(ls -d $DIR/* | awk -F "/" '{print $NF}')
do
    array[$x]=$elem
    x=`expr $x + 1`
done

y=0

for elem in $(ls -d $DIR2/* | awk -F "/" '{print $NF}')
do
    array2[$y]=$elem
    y=`expr $y + 1`
done

total=`expr $x \* $y`
curr=0

for ((i=0 ; i < ${#array[@]} ; i++))
do
    for ((j=0 ; j < ${#array2[@]} ; j++))
    do
                seqX=${array[$i]}
                seqY=${array2[$j]}
                echo "[GPU $DEV] ($curr/$total): ${seqX}-${seqY}"
                { time $BINDIR/gpu_cuda_workflow -query $DIR/${seqX} -ref $DIR2/${seqY} -dev $DEV -len $LEN > /dev/null 2> errors; } 2>&1 | grep "real" | awk '{print $2}'
                cat errors
                curr=`expr $curr + 1`
    done
done
