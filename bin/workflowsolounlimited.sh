#!/bin/bash

if [ $# != 6 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName seqYName lenght similarity WL fixedL"
   echo ""
   exit -1
fi


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
fixedL=${6}
strand=${7}

mkdir intermediateFiles

mkdir intermediateFiles/${seqXName}-${seqYName}
mkdir results
mkdir intermediateFiles/dictionaries
mkdir intermediateFiles/hits

# Copiamos los fastas
ln -s ${dirNameX}/${seqXName}.${extensionX} intermediateFiles/${seqXName}-${seqYName}
ln -s ${dirNameY}/${seqYName}.${extensionY} intermediateFiles/${seqXName}-${seqYName}

cd intermediateFiles/${seqXName}-${seqYName}

###############


echo "${BINDIR}/reverseComplement ${seqYName}.${extensionX} ${seqYName}-revercomp.${extensionY}"
${BINDIR}/reverseComplement ${seqYName}.${extensionX} ${seqYName}-revercomp.${extensionY}

if [[ ! -f ../dictionaries/${seqXName}.d2hP ]];	then
	echo "${BINDIR}/dictionarysolounlimited.sh ${seqXName}.${extensionX}"
	${BINDIR}/dictionarysolounlimited.sh ${seqXName}.${extensionX}
fi
		
if [[ ! -f ../dictionaries/${seqYName}.d2hP ]];	then
	echo "${BINDIR}/dictionarysolounlimited.sh ${seqYName}.${extensionY} "
	${BINDIR}/dictionarysolounlimited.sh ${seqYName}.${extensionY}
fi
		
if [[ ! -f ../dictionaries/${seqYName}-revercomp.d2hP ]];	then
	echo "${BINDIR}/dictionarysolounlimited.sh ${seqYName}-revercomp.${extensionY}"
	${BINDIR}/dictionarysolounlimited.sh ${seqYName}-revercomp.${extensionY} 
fi		




mv ${seqXName}.d2hP ../dictionaries/
mv ${seqXName}.d2hW ../dictionaries/
mv ${seqYName}.d2hP ../dictionaries/
mv ${seqYName}.d2hW ../dictionaries/
mv ${seqYName}-revercomp.d2hP ../dictionaries/
mv ${seqYName}-revercomp.d2hW ../dictionaries/
		
# Hacemos enlace simbolico
ln -s ../dictionaries/${seqXName}.d2hP .
ln -s ../dictionaries/${seqXName}.d2hW .

ln -s ../dictionaries/${seqYName}.d2hP .
ln -s ../dictionaries/${seqYName}.d2hW .

ln -s ../dictionaries/${seqYName}-revercomp.d2hP .
ln -s ../dictionaries/${seqYName}-revercomp.d2hW .

echo "${BINDIR}/comparisonsolounlimited.sh ${seqXName}.${extensionX} ${seqYName}.${extensionY} ${length} ${similarity} ${WL} ${fixedL} f"
${BINDIR}/comparisonsolounlimited.sh ${seqXName}.${extensionX} ${seqYName}.${extensionY} ${length} ${similarity} ${WL} ${fixedL} f

echo "${BINDIR}/comparisonsolounlimited.sh ${seqXName}.${extensionX} ${seqYName}-revercomp.${extensionY} ${length} ${similarity} ${WL} ${fixedL} r"
${BINDIR}/comparisonsolounlimited.sh ${seqXName}.${extensionX} ${seqYName}-revercomp.${extensionY} ${length} ${similarity} ${WL} ${fixedL} r


#echo "rm ${seqYName}-revercomp.${extensionY}"
#rm ${seqYName}-revercomp.${extensionY}

echo "${BINDIR}/combineFrags ${seqXName}-${seqYName}-sf.frags ${seqXName}-${seqYName}-revercomp-sr.frags ${seqXName}-${seqYName}.frags"
${BINDIR}/combineFrags ${seqXName}-${seqYName}-sf.frags ${seqXName}-${seqYName}-revercomp-sr.frags ${seqXName}-${seqYName}.frags

#echo "${BINDIR}/newFragToBalazsVersion ${seqXName}-${seqYName}.frags ${seqXName}-${seqYName}.old.frags"
#${BINDIR}/newFragToBalazsVersion ${seqXName}-${seqYName}.frags ${seqXName}-${seqYName}.old.frags

#echo "${BINDIR}/af2pngrev ${seqXName}-${seqYName}.frags ${seqXName}-${seqYName}.png ${seqXName} ${seqYName}"
#${BINDIR}/af2pngrev ${seqXName}-${seqYName}.frags ${seqXName}-${seqYName}.png ${seqXName} ${seqYName}

#Borramos todo menos los frags y los diccionarios

# Get Info from frags 
echo "${BINDIR}/getInfo ${seqXName}-${seqYName}.frags > ${seqXName}-${seqYName}.csv"
${BINDIR}/getInfo ${seqXName}-${seqYName}.frags > ${seqXName}-${seqYName}.csv.tmp
cat ${seqXName}-${seqYName}.frags.INF ${seqXName}-${seqYName}.csv.tmp > ${seqXName}-${seqYName}.csv
rm -rf ${seqXName}-${seqYName}.csv.tmp
	
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
mv ${seqXName}-${seqYName}.frags.MAT ../../results
#mv ${seqXName}-${seqYName}.old.frags ../../results
mv ${seqXName}-${seqYName}.csv ../../results

#echo "Borrando ${seqXName}-${seqYName}"
cd ..
#rm -rf ${seqXName}-${seqYName}

#} &> /dev/null
