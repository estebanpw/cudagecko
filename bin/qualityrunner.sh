#!/bin/bash

if [ $# -lt 3 ]; then
    echo "***ERROR*** Use: $0 <csv> <fastax> <fastay> [1 for printing ungapped alignment]"
    exit -1
fi

CSV=$1
FASTAX=$2
FASTAY=$3
PRINTALIGN=$4

mkdir subfastas

while IFS= read -r line; do

    coords=$(echo $line | awk 'BEGIN{FS=","}{print $2,$4,$3,$5,$6}')
    arraycoords=(${coords})

    xstart=${arraycoords[0]}
    xend=${arraycoords[1]}
    ystart=${arraycoords[2]}
    yend=${arraycoords[3]}
    strand=${arraycoords[4]}




    if [[ $strand == "f" ]]; then

        len=`expr $xend - $xstart`
        echo "> $2 [$xstart, $xend]"> subfastas/extractX-$xstart-$xend.fasta
        echo "> $3 [$ystart, $yend]"> subfastas/extractY-$ystart-$yend.fasta
        grep -v ">" $FASTAX | tr -d '\n' | tail -c +$xstart | head -c $len >> subfastas/extractX-$xstart-$xend.fasta
        grep -v ">" $FASTAY | tr -d '\n' | tail -c +$ystart | head -c $len >> subfastas/extractY-$ystart-$yend.fasta

        #echo "Aligning $xstart-$xend with $ystart-$yend"
        echo "needle -asequence subfastas/extractX-$xstart-$xend.fasta -bsequence subfastas/extractY-$ystart-$yend.fasta -outfile subfastas/align-$xstart-$ystart.align"
        needle -asequence subfastas/extractX-$xstart-$xend.fasta -bsequence subfastas/extractY-$ystart-$yend.fasta -outfile subfastas/align-$xstart-$ystart.align -gapopen 10 -gapextend 0.5 &> /dev/null
        echo "Gapped alignment"
        cat subfastas/align-$xstart-$ystart.align | grep "Identity\|Similarity\|Gaps"
        echo "Ungapped alignment"
        if [[ $PRINTALIGN == "1" ]]; then
            awk 'NR==FNR && FNR==2{l1=$0} NR!=FNR && FNR==2{l2=$0} 
            END{good=0; split(l1, s1, ""); split(l2, s2, ""); print l1; 
            for(i=0; i<length(l1); i++){ if(s1[i]==s2[i]) {good++; printf("|");} else printf(" ");  }; printf("\n");
            print l2; print "@ Identity:",good "/" length(l1), "("100*good/length(l1)"%)";  }' subfastas/extractX-$xstart-$xend.fasta subfastas/extractY-$ystart-$yend.fasta
        else

            awk 'NR==FNR && FNR==2{l1=$0} NR!=FNR && FNR==2{l2=$0} 
            END{good=0; split(l1, s1, ""); split(l2, s2, "");  
            for(i=0; i<length(l1); i++){ if(s1[i]==s2[i])good++; }; 
            print "@ Identity:",good "/" length(l1), "("100*good/length(l1)"%)"  }' subfastas/extractX-$xstart-$xend.fasta subfastas/extractY-$ystart-$yend.fasta
        fi

    fi

done < <(tail -n +18 $CSV)

#coords=$(tail -n +18 $CSV | head -1 |  awk 'BEGIN{FS=","}{print $2,$4,$3,$5}')




