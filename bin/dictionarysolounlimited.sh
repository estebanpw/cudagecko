#!/bin/bash 

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# != 1 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName.fasta"
   echo ""
   exit -1
fi

seqName=$(basename "$1")
extension="${seqName##*.}"
seqName="${seqName%.*}"

# find words and order
echo "${BINDIR}/words $1 ${seqName}.words.unsort"
/usr/bin/time -f "%e" ${BINDIR}/words $1 ${seqName}.words.unsort
echo "${BINDIR}/sortWords 10000000 1 ${seqName}.words.unsort ${seqName}.words.sort"
/usr/bin/time -f "%e" ${BINDIR}/sortWords 10000000 1 ${seqName}.words.unsort ${seqName}.words.sort

# Create hash table in disk
echo "${BINDIR}/w2hd ${seqName}.words.sort ${seqName}"
/usr/bin/time -f "%e" ${BINDIR}/w2hd ${seqName}.words.sort ${seqName}

