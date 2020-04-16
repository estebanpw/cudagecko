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
time ${BINDIR}/words $1 ${seqName}.words.unsort
echo "${BINDIR}/sortWords 10000000 1 ${seqName}.words.unsort ${seqName}.words.sort"
time ${BINDIR}/sortWords 10000000 1 ${seqName}.words.unsort ${seqName}.words.sort

# Create hash table in disk
echo "${BINDIR}/w2hd ${seqName}.words.sort ${seqName}"
time ${BINDIR}/w2hd ${seqName}.words.sort ${seqName}

