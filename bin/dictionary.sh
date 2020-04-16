#!/bin/bash 

FL=1000   # frequency limit
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
${BINDIR}/words $1 ${seqName}.words.unsort
echo "${BINDIR}/sortWords 10000000 32 ${seqName}.words.unsort ${seqName}.words.sort"
${BINDIR}/sortWords 10000000 32 ${seqName}.words.unsort ${seqName}.words.sort

# Create hash table in disk
echo "${BINDIR}/w2hd ${seqName}.words.sort ${seqName}"
${BINDIR}/w2hd ${seqName}.words.sort ${seqName}

