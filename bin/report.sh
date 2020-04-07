#!/bin/bash

if [ $# -lt 3 ]; then
    echo "***ERROR*** Use: $0 <fastax> <fastay> <minlen>"
    exit -1
fi


FASTAX=$1
FASTAY=$2
MINLEN=$3

nameX=$(basename $FASTAX)
pureNameX="${nameX%.*}"
nameY=$(basename $FASTAY)
pureNameY="${nameY%.*}"

BINGPU=/home/estebanpw/blueberry/openCL/cudagecko/ggecko/bin
BINBLA=/home/estebanpw/blueberry/openCL/Gblastn/gblastn/c++/GCC480-ReleaseMT64/bin
BINGEC=/home/estebanpw/blueberry/openCL/cudagecko/ggecko/bin
BINQUA=/home/estebanpw/blueberry/openCL/scripts
BINBLASTTOGECKO=/home/estebanpw/blueberry/openCL/scripts
BINCOV=/home/estebanpw/blueberry/openCL/coverage-flattener
RES=/home/estebanpw/blueberry/res/reports/$nameX-$nameY

mkdir -p $RES


REPORT=$RES/report-$nameX-$nameY.log

echo "#REPORT" > $REPORT
echo "#" >> $REPORT



# Starting report

echo "# GECKO 1 CORE UNLIMITED HITS"
echo "# GECKO 1 CORE UNLIMITED HITS" >> $REPORT
echo "/usr/bin/time -v $BINGEC/workflowsolounlimited.sh $FASTAX $FASTAY $MINLEN 80 32 1" >> $REPORT

(/usr/bin/time -v $BINGEC/workflowsolounlimited.sh $FASTAX $FASTAY $MINLEN 80 32 1) &>> $REPORT &

pid=$!

while true;
do

    count=$(ps -A| grep $pid |wc -l)

    if [[ $count -eq 0 ]] #if process is already terminated, then there can be two cases, the process executed and stop successfully or it is terminated abnormally
    then
        break
    else
        MBsFree=$(df -BM | grep "home" | awk '{print $4}' | sed 's/M//g')
        echo "$pid has still $MBsFree MBs to use"

        if [[ $MBsFree -lt 10000 ]]
        then
            echo "NO MORE SPACE, ABORTING!"
            kill $pid
            rm -rf intermediateFiles results
            exit
        fi
    fi

    sleep 5

done

cp results/$pureNameX-$pureNameY.csv $RES
rm -rf intermediateFiles results

echo "# RUNNING QUALITY FOR GECKO"
echo "# RUNNING QUALITY FOR GECKO" >> $REPORT
echo "$BINQUA/qualityrunner.sh $RES/$pureNameX-$pureNameY.csv $FASTAX $FASTAY ungapped 0 > $RES/$pureNameX-$pureNameY.gecko.quality " >> $REPORT
$BINQUA/qualityrunner.sh $RES/$pureNameX-$pureNameY.csv $FASTAX $FASTAY ungapped 0 > $RES/$pureNameX-$pureNameY.gecko.quality 


echo "# GBLASTN NO DUST NO MASK"
echo "# GBLASTN NO DUST NO MASK" >> $REPORT
echo "# MAKEBLASTDB IS INCLUDED" >> $REPORT
echo "time $BINBLA/makeblastdb -dbtype nucl -in $FASTAY" >> $REPORT

dbtime=$({ time $BINBLA/makeblastdb -dbtype nucl -in $FASTAY >/dev/null 2> errors; } 2>&1 | grep real | awk '{print $2}')

echo "# TIME MAKEBLASTDB" >> $REPORT
echo "$dbtime" >> $REPORT


echo "/usr/bin/time -v $BINBLA/blastn -task megablast -dust no -soft_masking false -ungapped -perc_identity 80 -num_alignments 1000000 -word_size 32 -query $FASTAX -db $FASTAY -outfmt "7 qlen slen qstart qend sstart send length nident pident" -out $RES/$nameX-$nameY.gblastn" >> $REPORT

(/usr/bin/time -v $BINBLA/blastn -task megablast -dust no -soft_masking false -ungapped -perc_identity 80 -num_alignments 1000000 -word_size 32 -query $FASTAX -db $FASTAY -outfmt "7 qlen slen qstart qend sstart send length nident pident" -out $RES/$nameX-$nameY.gblastn) &>> $REPORT

grep -v "#" $RES/$nameX-$nameY.gblastn | awk -v len="$MINLEN" '{if($7 >= len) print $0}' > temp
mv temp $RES/$nameX-$nameY.gblastn


echo "# RUNNING BLAST TO GECKO CONVERSOR" >> $REPORT

$BINBLASTTOGECKO/blastntocsv.sh $RES/$nameX-$nameY.gblastn

echo "# RUNNING QUALITY FOR GBLASTN"
echo "# RUNNING QUALITY FOR GBLASTN" >> $REPORT

echo "$BINQUA/qualityrunner.sh $RES/$nameX-$nameY.gblastn.csv $FASTAX $FASTAY ungapped 1 > $RES/$nameX-$nameY.gblastn.quality " >> $REPORT
$BINQUA/qualityrunner.sh $RES/$nameX-$nameY.gblastn.csv $FASTAX $FASTAY ungapped 1 > $RES/$nameX-$nameY.gblastn.quality 


echo "# GPUGECKO "
echo "# GPUGECKO " >> $REPORT
echo "/usr/bin/time -v $BINGPU/gpu_cuda_workflow -query $FASTAX -ref $FASTAY -dev 3 -len $MINLEN" >> $REPORT
(/usr/bin/time -v $BINGPU/gpu_cuda_workflow -query $FASTAX -ref $FASTAY -dev 3 -len $MINLEN) &>> $REPORT

mv $nameX-$nameY.csv $RES/$nameX-$nameY.gpugecko.csv

echo "# RUNNING QUALITY FOR GPUGECKO"
echo "# RUNNING QUALITY FOR GPUGECKO" >> $REPORT
echo "$BINQUA/qualityrunner.sh $RES/$nameX-$nameY.gpugecko.csv $FASTAX $FASTAY ungapped 1 > $RES/$nameX-$nameY.gpugecko.quality " >> $REPORT

$BINQUA/qualityrunner.sh $RES/$nameX-$nameY.gpugecko.csv $FASTAX $FASTAY ungapped 1 > $RES/$nameX-$nameY.gpugecko.quality 

echo "# WAITING FOR JOBS" >> $REPORT

for job in `jobs -p`
do
    #echo $job
    wait $job
done

echo "# RUNNING COVERAGE ANALYSIS" >> $REPORT

$BINCOV/basic-coverage $RES/$pureNameX-$pureNameY.csv > $RES/$pureNameX-$pureNameY.gecko.cov
$BINCOV/basic-coverage $RES/$nameX-$nameY.gblastn.csv > $RES/$nameX-$nameY.gblastn.cov
$BINCOV/basic-coverage $RES/$nameX-$nameY.gpugecko.csv > $RES/$nameX-$nameY.gpugecko.cov

covGECKO=$(cat $RES/$pureNameX-$pureNameY.gecko.cov | grep "TBases")
covBLAST=$(cat $RES/$nameX-$nameY.gblastn.cov | grep "TBases")
covGPUGK=$(cat $RES/$nameX-$nameY.gpugecko.cov | grep "TBases")


qualityGPUGECKO=$(grep "@" $RES/$nameX-$nameY.gpugecko.quality | awk '{print $4}' | sed 's/(//g' | sed 's/)//g' | sed 's/%//g'  | awk 'BEGIN{min=100; max=0; avg=0;}{ if($0 > max)max=$0; if($0<min)min=$0; avg += $0; } END{print "min", min, "max", max, "avg", avg/NR}')
qualityGECKO=$(grep "@" $RES/$pureNameX-$pureNameY.gecko.quality | awk '{print $4}' | sed 's/(//g' | sed 's/)//g' | sed 's/%//g'  | awk 'BEGIN{min=100; max=0; avg=0;}{ if($0 > max)max=$0; if($0<min)min=$0; avg += $0; } END{print "min", min, "max", max, "avg", avg/NR}')
qualityGBLASTN=$(grep "@" $RES/$nameX-$nameY.gblastn.quality | awk '{print $4}' | sed 's/(//g' | sed 's/)//g' | sed 's/%//g'  | awk 'BEGIN{min=100; max=0; avg=0;}{ if($0 > max)max=$0; if($0<min)min=$0; avg += $0; } END{print "min", min, "max", max, "avg", avg/NR}')

echo "Total coverage GECKO:    $covGECKO"
echo "Quality GECKO:           $qualityGECKO"
echo "Quality GECKO:           $qualityGECKO" >> $REPORT
echo "Total coverage GBLASTN:  $covBLAST"
echo "Quality GBLASTN:         $qualityGBLASTN"
echo "Quality GBLASTN:         $qualityGBLASTN" >> $REPORT
echo "Total coverage GPUGECKO: $covGPUGK"
echo "Quality GPUGECKO:        $qualityGPUGECKO"
echo "Quality GPUGECKO:        $qualityGPUGECKO" >> $REPORT

echo "Total coverage GECKO:    $covGECKO" >> $REPORT
echo "Total coverage GBLASTN:  $covBLAST" >> $REPORT
echo "Total coverage GPUGECKO: $covGPUGK" >> $REPORT


echo "# REPORT COMPLETED" >> $REPORT

echo "FINISHED!"
