#!/bin/bash

FL=1000   # frequency limit

if [ $# != 7 ]; then
   echo " ==== ERROR ...."
   echo ""
   echo "   usage:  $0 seqXName seqYName lenght similarity WL evalue GPUdevice"
   echo ""
   exit -1
fi


# {

dirNameX=$(readlink -f $1 | xargs dirname)
seqXName=$(basename "$1")
extensionX="${seqXName##*.}"
seqXName="${seqXName%.*}"

dirNameY=$(readlink -f $2 | xargs dirname)
seqYName=$(basename "$2")
extensionY="${seqYName##*.}"
seqYName="${seqYName%.*}"

#seqXName=`basename $1 .fasta`
#seqYName=`basename $2 .fasta`

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

length=${3}
similarity=${4}
WL=${5} # wordSize
evalue=${6}
gpu=${7}
strand=${8}
fixedL=1

mkdir intermediateFiles

mkdir intermediateFiles/${seqXName}-${seqYName}
mkdir results

# Copiamos los fastas
ln -s ${dirNameX}/${seqXName}.${extensionX} intermediateFiles/${seqXName}-${seqYName}
ln -s ${dirNameY}/${seqYName}.${extensionY} intermediateFiles/${seqXName}-${seqYName}

cd intermediateFiles/${seqXName}-${seqYName}

###############


echo "${BINDIR}/reverseComplement ${seqYName}.${extensionX} ${seqYName}-revercomp.${extensionY}"
${BINDIR}/reverseComplement ${seqYName}.${extensionY} ${seqYName}-revercomp.${extensionY}

# Run forward

time ${BINDIR}/gpu_workflow -q ${seqXName}.${extensionX} -r ${seqYName}.${extensionY} -frags ${seqXName}-${seqYName}-sf.frags -power 24 -dev $gpu -l $length -s $similarity -kmer $WL

# Run reverse

time ${BINDIR}/gpu_workflow -q ${seqXName}.${extensionX} -r ${seqYName}-revercomp.${extensionY} -frags ${seqXName}-${seqYName}-revercomp-sr.frags -power 24 -dev $gpu -l $length -s $similarity --r -kmer $WL


echo "${BINDIR}/combineFrags ${seqXName}-${seqYName}-sf.frags ${seqXName}-${seqYName}-revercomp-sr.frags ${seqXName}-${seqYName}.frags"
${BINDIR}/combineFrags ${seqXName}-${seqYName}-sf.frags ${seqXName}-${seqYName}-revercomp-sr.frags ${seqXName}-${seqYName}.frags


# Get Info from frags 
echo "${BINDIR}/gpu_convert_frags_csv ${seqXName}-${seqYName}.frags ${length} ${similarity} > ${seqXName}-${seqYName}.csv.tmp"
${BINDIR}/gpu_convert_frags_csv ${seqXName}-${seqYName}.frags ${length} ${similarity} ${seqXName} ${seqYName} > ${seqXName}-${seqYName}.csv

	
if [[ -L "../../${seqXName}.fasta" ]]
then
	rm ../../${seqXName}.fasta
fi

if [[ -L "../../${seqYName}.fasta" ]]
then
	rm ../../${seqYName}.fasta
fi

#Movemos los frags y los info
mv ${seqXName}-${seqYName}.frags ../../results
mv ${seqXName}-${seqYName}.frags.INF ../../results
#mv ${seqXName}-${seqYName}.old.frags ../../results
mv ${seqXName}-${seqYName}.csv ../../results

#echo "Borrando ${seqXName}-${seqYName}"
cd ..
#rm -rf ${seqXName}-${seqYName}

#} &> /dev/null
