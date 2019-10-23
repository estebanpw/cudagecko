#!/usr/bin/env bash
DIR=$1
L=$2
S=$3
WL=$4
EXT=$5

array=()
x=0

if [ $# != 5 ]; then
	echo "***ERROR*** Use: $0 genomesDirectory L(200) S(40) K(8) fastaFilesExtension"
	exit -1
fi

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for elem in $(ls -d $DIR/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	array[$x]=$elem
	x=`expr $x + 1`
	#echo "X: $elem"
done

for ((i=0 ; i < ${#array[@]} ; i++))
do
	for ((j=i ; j < ${#array[@]} ; j++))
	do
		if [ $i != $j ]; then
				seqX=${array[$i]}
				seqY=${array[$j]}
				echo "----------${seqX}-${seqY}-----------"
			if [[ ! -f frags/${seqX}-${seqY}.frags ]];	then
				
				echo "${BINDIR}/workflow.sh $DIR/${seqX}.$EXT $DIR/${seqY}.$EXT $L $S $WL 1"
				${BINDIR}/workflow.sh $DIR/${seqX}.$EXT $DIR/${seqY}.$EXT $L $S $WL 1
			
			fi
		fi
	done
done
