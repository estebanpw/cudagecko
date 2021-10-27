#!/usr/bin/env bash

set +e # disables exit from error 

SEQX=$1
SEQY=$2
LEN=$3
DEV=$4

if [ $# -lt 4 ]; then
	echo "***ERROR*** Use: $0 <query fasta> <ref fasta> <min len{32,64,92,etc}> <device{0,1,..n}>"
	exit -1
fi

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


myfactor="0.15"

$BINDIR/gpu_cuda_workflow -query $SEQX -ref $SEQY -dev $DEV -len $LEN -factor $myfactor 2>error

run=1


S1="Reached maximum limit of hits"
S2="Not enough memory in pool"
S3="Error generating hits"

if grep -qF "$S1" error || grep -qF "$S2" error || grep -qF "$S3" error ; then

	myfactor="0.12"
	while grep -qF "$S1" error || grep -qF "$S2" error || grep -qF "$S3" error ; do
		
		cat error
		echo "[WARNING] Restarting execution with factor=$myfactor"
		
		rm error
		
		$BINDIR/gpu_cuda_workflow -query $SEQX -ref $SEQY -factor $myfactor -len $LEN -dev $DEV 2>error
		
		outtext="[SUCCESS] Completed with factor $myfactor"
		
		run=`expr $run + 1`
		myfactor=$(echo $myfactor | awk '{print $1/2}')
		
		
	done

	echo "$outtext"

elif [ $(wc -c error | awk '{print $1}') -gt 1 ]; then

	echo "[ERROR] Encountered the following error:"
	cat error

else

	outtext="[SUCCESS] Completed with factor $myfactor"

fi

rm error
			











