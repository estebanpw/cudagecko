#!/bin/bash

if [ $# -lt 3 ]; then
    echo "***ERROR*** Use: $0 <csv> <fastax> <fastay>"
    exit -1
fi

CSV=$1
FASTAX=$2
FASTAY=$3

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

REVCOM=$BINDIR
EASYX=$(basename $FASTAX)
EASYY=$(basename $FASTAY)
SUBFASTAS="subfastas-$EASYX-$EASYY"

(mkdir $SUBFASTAS) &> /dev/null

grep -v ">" $FASTAX | tr -d '\n' > $SUBFASTAS/x.fasta
grep -v ">" $FASTAY | tr -d '\n' > $SUBFASTAS/y.fasta

$REVCOM/reverseComplement $FASTAY $SUBFASTAS/yrev.temp
grep -v ">" $SUBFASTAS/yrev.temp | tr -d '\n' > $SUBFASTAS/yrev.fasta

FILTERX=./$SUBFASTAS/x.fasta
FILTERY=./$SUBFASTAS/y.fasta

lenX=$(wc -c $FILTERX | awk 'BEGIN{FS=" "}{print $1}')
lenY=$(wc -c $FILTERY | awk 'BEGIN{FS=" "}{print $1}')

tAlignments=0

while IFS= read -r line; do

    coords=$(echo $line | awk 'BEGIN{FS=","}{print $2+1,$4+1,$3+1,$5+1,$6}')
    arraycoords=(${coords})

    xstart=${arraycoords[0]}
    xend=${arraycoords[1]}
    ystart=${arraycoords[2]}
    yend=${arraycoords[3]}
    strand=${arraycoords[4]}




    if [[ $strand == "f" ]]; then

        len=`expr $xend - $xstart`
        echo "> $2 [$xstart, $xend]"> $SUBFASTAS/extractX-$xstart-$xend.fasta
        echo "> $3 [$ystart, $yend]"> $SUBFASTAS/extractY-$ystart-$yend.fasta
        tail -c +$xstart $SUBFASTAS/x.fasta | head -c $len >> $SUBFASTAS/extractX-$xstart-$xend.fasta
        tail -c +$ystart $SUBFASTAS/y.fasta | head -c $len >> $SUBFASTAS/extractY-$ystart-$yend.fasta

        echo "# Align. $tAlignments; query: subfastas/extractX-$xstart-$xend.fasta ref: subfastas/extractY-$ystart-$yend.fasta"
        awk 'NR==FNR && FNR==2{l1=$0} NR!=FNR && FNR==2{l2=$0} 
                END{good=0; split(l1, s1, ""); split(l2, s2, ""); print l1; 
                for(i=1; i<=length(l1); i++){ if(s1[i]==s2[i]) {good++; printf("|");} else printf(" ");  }; printf("\n");
                print l2; print "@ Identity:",good "/" length(l1), "("100*good/length(l1)"%, forward)";  }' $SUBFASTAS/extractX-$xstart-$xend.fasta $SUBFASTAS/extractY-$ystart-$yend.fasta
        rm $SUBFASTAS/extractX-$xstart-$xend.fasta $SUBFASTAS/extractY-$ystart-$yend.fasta

    else

        len=`expr $xend - $xstart`
        # Undo the y change
        
        ystart=`expr $lenY - $ystart`
        ystart=`expr $ystart + 1`

        yend=`expr $lenY - $yend`
        yend=`expr $yend + 1`
    

        echo "> $2 [$xstart, $xend]"> $SUBFASTAS/extractX-$xstart-$xend.fasta
        echo "> $3 [$ystart, $yend]"> $SUBFASTAS/extractYrev-$ystart-$yend.fasta
        tail -c +$xstart $SUBFASTAS/x.fasta | head -c $len >> $SUBFASTAS/extractX-$xstart-$xend.fasta
        tail -c +$ystart $SUBFASTAS/yrev.fasta | head -c $len >> $SUBFASTAS/extractYrev-$ystart-$yend.fasta

        echo "# Align. $tAlignments; query: subfastas/extractX-$xstart-$xend.fasta ref: subfastas/extractYrev-$ystart-$yend.fasta"
        #echo "Aligning $xstart-$xend with $ystart-$yend"
        awk 'NR==FNR && FNR==2{l1=$0} NR!=FNR && FNR==2{l2=$0}
                END{good=0; split(l1, s1, ""); split(l2, s2, ""); print l1;
                for(i=1; i<=length(l1); i++){ if(s1[i]==s2[i]) {good++; printf("|");} else printf(" ");  }; printf("\n");
                print l2; print "@ Identity:",good "/" length(l1), "("100*good/length(l1)"%, reverse)";  }' $SUBFASTAS/extractX-$xstart-$xend.fasta $SUBFASTAS/extractYrev-$ystart-$yend.fasta
        
        rm $SUBFASTAS/extractX-$xstart-$xend.fasta $SUBFASTAS/extractYrev-$ystart-$yend.fasta

    fi

    tAlignments=`expr $tAlignments + 1`
    echo " "
   

done < <(tail -n +18 $CSV)

rm $SUBFASTAS/*

#coords=$(tail -n +18 $CSV | head -1 |  awk 'BEGIN{FS=","}{print $2,$4,$3,$5}')





