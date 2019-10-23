#!/bin/bash 

FL=1000   # frequency limit

if [ $# != 8 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName seqYName lenght similarity WL fixedL strand evalue"
   echo ""
   exit -1
fi

seqXName=$(basename "$1")
extensionX="${seqXName##*.}"
seqXName="${seqXName%.*}"

seqYName=$(basename "$2")
extensionY="${seqYName##*.}"
seqYName="${seqYName%.*}"

#seqXName=`basename $1 .fasta`
#seqYName=`basename $2 .fasta`

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

length=${3}
similarity=${4}
WL=${5} # wordSize
fixedL=${6}
strand=${7}
evalue=${8}
distance=$((4*${WL}))

if [[ ! -f ../hits/${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered ]]; then
	echo "(time ${BINDIR}/hits ${seqXName} ${seqYName} ${seqXName}-${seqYName}-K${WL}.hits ${FL} ${WL}) &> log-hits-$strand"
	(time ${BINDIR}/hits ${seqXName} ${seqYName} ${seqXName}-${seqYName}-K${WL}.hits ${FL} ${WL}) &> log-hits-$strand
	
	#echo "${BINDIR}/sortHits 10000000 32 ${seqXName}-${seqYName}-K${WL}.hits ${seqXName}-${seqYName}-K${WL}.hits.sorted"
	#${BINDIR}/sortHits 10000000 32 ${seqXName}-${seqYName}-K${WL}.hits ${seqXName}-${seqYName}-K${WL}.hits.sorted
	
	#echo "(time ${BINDIR}/gpu_sort_hits -f ${seqXName}-${seqYName}-K${WL}.hits -wpi 1 -power 23 -dev 3) &> log-gpu-sorthits-$strand"
	#(time ${BINDIR}/gpu_sort_hits -f ${seqXName}-${seqYName}-K${WL}.hits -wpi 1 -power 23 -dev 3) &> log-gpu-sorthits-$strand

	
	#echo "${BINDIR}/filterHits ${seqXName}-${seqYName}-K${WL}.hits.sorted ${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered ${WL}"
	#(time ${BINDIR}/filterHits ${seqXName}-${seqYName}-K${WL}.hits.sorted ${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered ${WL}) &> log-filterhits-$strand

	#mv ${seqXName}-${seqYName}-K${WL}.hits.sorted ../hits/${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered
	#mv ${seqXName}-${seqYName}-K${WL}.hits.sorted.nseqs ../hits/${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered.nseqs
fi

#ln -s ../hits/${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered .
#ln -s ../hits/${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered.nseqs .

#echo "${BINDIR}/FragHits $1 $2 ${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered ${seqXName}-${seqYName}-s${strand}.frags ${length} ${similarity} ${WL} ${fixedL} ${strand}"
#${BINDIR}/FragHits $1 $2 ${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered ${seqXName}-${seqYName}-s${strand}.frags ${length} ${similarity} ${WL} ${fixedL} ${strand}

if [[ "$strand" == "f" ]]; then

	echo "(time ${BINDIR}/gpu_sort_hits_plus_frags -q $1 -r $2 -hits ${seqXName}-${seqYName}-K${WL}.hits -out ${seqXName}-${seqYName}-s${strand}.frags -power 23 -wpi 1 -kmer 32 -dev 3 -l ${length} -s ${similarity}) &> log-gpu-frags-$strand"
	(time ${BINDIR}/gpu_sort_hits_plus_frags -q $1 -r $2 -hits ${seqXName}-${seqYName}-K${WL}.hits -out ${seqXName}-${seqYName}-s${strand}.frags -power 23 -wpi 1 -kmer 32 -dev 3 -l ${length} -s ${similarity}) &> log-gpu-frags-$strand

	#echo "(time ${BINDIR}/gpu_frags -q $1 -r $2 -hits ${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered -l ${length} -s ${similarity} -e ${evalue} -dev 3) &> log-gpu-frags-$strand"
	#(time ${BINDIR}/gpu_frags -q $1 -r $2 -hits ${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered -l ${length} -s ${similarity} -e ${evalue} -dev 3) &> log-gpu-frags-$strand

	#mv ${seqXName}.${extensionX}-${seqYName}.${extensionY}-f.frags ${seqXName}-${seqYName}-s${strand}.frags
	#mv ${seqXName}.${extensionX}-${seqYName}.${extensionY}-f.frags.INF ${seqXName}-${seqYName}-s${strand}.frags.INF
fi

if [[ "$strand" == "r" ]]; then

	echo "(time ${BINDIR}/gpu_sort_hits_plus_frags -q $1 -r $2 -hits ${seqXName}-${seqYName}-K${WL}.hits -out ${seqXName}-${seqYName}-s${strand}.frags -power 23 -wpi 1 -kmer 32 -dev 3 -l ${length} -s ${similarity}) &> log-gpu-frags-$strand"
	(time ${BINDIR}/gpu_sort_hits_plus_frags -q $1 -r $2 -hits ${seqXName}-${seqYName}-K${WL}.hits -out ${seqXName}-${seqYName}-s${strand}.frags -power 23 -wpi 1 -kmer 32 -dev 3 -l ${length} -s ${similarity}) &> log-gpu-frags-$strand
	#echo "(time ${BINDIR}/gpu_frags -q $1 -r $2 -hits ${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered -l ${length} -s ${similarity} -e ${evalue} --r -dev 3) &> log-gpu-frags-$strand"
	#(time ${BINDIR}/gpu_frags -q $1 -r $2 -hits ${seqXName}-${seqYName}-K${WL}.hits.sorted.filtered -l ${length} -s ${similarity} -e ${evalue} --r -dev 3) &> log-gpu-frags-$strand

	mv ${seqXName}.${extensionX}-${seqYName}.${extensionY}-r.frags ${seqXName}-${seqYName}-s${strand}.frags
fi

echo "--------------------DONE------------------"

