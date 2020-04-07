#!/bin/bash

if [[ $# -ne 1 ]]; then

    echo "Use: $0 <input.blast>"
    echo "Blast must be used with time blastn -task megablast -ungapped -use_gpu true -gpu_id 0 -perc_identity 80 -num_alignments 1000000 -word_size 32 -outfmt "\""7 qlen slen qstart qend sstart send length nident pident"\"" -query ... -db ... -outfmt 7 -out ..."
    exit

fi

lenX=$(grep -v "#" $1 | awk '{print $1}' | head -1)
lenY=$(grep -v "#" $1 | awk '{print $2}' | head -1)

echo "All by-Identity Ungapped Fragments (Hits based approach)
[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>
SeqX filename : undef
SeqY filename : undef
SeqX name : undef
SeqY name : undef
SeqX length : $lenX
SeqY length : $lenY
Min.fragment.length : undef
Min.Identity : undef
Tot Hits (seeds) : undef
Tot Hits (seeds) used: undef
Total fragments : undef
========================================================
Total CSB: 0
========================================================
Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY" > $1.csv

grep -v "#" $1 | awk -v leny="$lenY" 'BEGIN{OFS=","}{ if($5 < $6) print "Frag",$3,$5,$4,$6,"f",0,$8,$8*$9/100,$8*$9/100,$9,$9,0,0; else { newystart=$6-2; newyend=$5-2; print "Frag",$3,newyend,$4,newystart,"r",0,$8,$8*$9/100,$8*$9/100,$9,$9,0,0;}}' >> $1.csv
