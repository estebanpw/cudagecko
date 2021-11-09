#!/bin/bash

if [ $# != 2 ]; then
   echo " ==== ERROR === "
   echo ""
   echo "   usage:  $0 <fasta> <splits>"
   echo ""
   exit -1
fi

FASTA=$1
SPLITS=$2

NAME=$(basename "$FASTA")
EXTENSION="${NAME##*.}"
NAME="${NAME%.*}"

# First step, count bytes of sequence
grep -v ">" $FASTA > $FASTA.nohead
TBYTES=$(wc -c $FASTA.nohead)
echo "Sequence has $TBYTES bytes"
SPLITBYTES=$(awk -v var1="$TBYTES" -v var2="$SPLITS" 'BEGIN { print  int( var1 / var2 ) }')

# Second step, grep name of sequence
SEQHEADER=$(grep ">" $FASTA)
echo "Splitting $SEQHEADER into $2 splits of $SPLITBYTES bytes each"

# Third step: separate into sequences
START=1
MINUSONE=`expr $SPLITS - 1`
for((i=0; i<$SPLITS; i++)); do
        echo "$SEQHEADER split $i" > $NAME.split$i.fasta
		if [ "$i" -eq "$MINUSONE"  ]; then
			# Last split
			tail -c +$START $FASTA.nohead >> $NAME.split$i.fasta
		else
	        tail -c +$START $FASTA.nohead | head -c $SPLITBYTES >> $NAME.split$i.fasta
		fi
        START=`expr $START + $SPLITBYTES`

done

rm $FASTA.nohead
