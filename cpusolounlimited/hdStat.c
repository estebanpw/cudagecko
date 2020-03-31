/* leehd read and displays the hash table from disk 
   Syntax: leehd prefixNameOUT

    prefixNameOUT.h2dW  : index of words-Pos-Ocurrences
    prefixNameOUT.h2dP  : positions
    both must be available

    Any char as third argument means "Verbose mode"
    Feb.2012: computes word frequencies	

                              ortrelles@uma.es / Dic.2011
    ---------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <inttypes.h>

#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"

#define PEQ 1001

int main(int ac, char** av){

	char fname[1024], *W;
	W=(char *)malloc(33*sizeof(char));
	FILE *f1, *f2, *f3;
	hashentry he;
	uint64_t i=0;
        location spos;
	uint64_t nW=0,maxF=0, aveF=0;
	int flagV=0;
	int64_t freq[PEQ];

	if(ac<2)terror("USE: leehd  prefixNameOUT [v=verbose]\n");
	if (ac==3) flagV=1;
	for (i=0;i<PEQ;i++) freq[i]=0;
	sprintf(fname,"%s.d2hW",av[1]); // Words file (first level of hash table)
	if ((f1 = fopen(fname,"rb"))==NULL) terror("opening prefix.h2dW file");
	sprintf(fname,"%s.d2hP",av[1]); // Positions file
	if ((f2 = fopen(fname,"rb"))==NULL) terror("opening prefix.h2dP file");

	sprintf(fname,"%s.freq",av[1]); // output
	if ((f3 = fopen(fname,"wt"))==NULL) terror("opening prefix.freq OUT file");

	// kick-off
	if(fread(&he,sizeof(hashentry),1,f1)!=1)
		terror("Empty dictionary");
       
        while(!feof(f1)){

             if (flagV) {showWord(&he.w, W);fprintf(stdout, "%.32s", W);}
             if (flagV) fprintf(stdout,"  : num=%-7" PRIu64 ":",he.num);
		if (he.num>=PEQ) {
			fprintf(f3, "%" PRIu64 "\t", he.num);
			showWord(&he.w, W);
			fprintf(f3, "%.32s", W);
			fprintf(f3, "%" PRIu64 "\n", he.num);
		}
		else freq[he.num]++;
	     	nW++; 
		if (he.num>maxF) maxF=he.num;
		aveF+=he.num;

             fseek(f2,0, he.pos);
	     if (flagV) {

             for (i=0;i<he.num;i++){
	       if(fread(&spos,sizeof(location),1,f2)!=1)
			terror("Error reading the word occurrences");
               fprintf(stdout,"(%" PRIu64 ",%" PRIu64 ") ",spos.pos,spos.seq);
             }
             fprintf(stdout,"\n");
	    }
	     if(fread(&he,sizeof(hashentry),1,f1)!=1)
		if(ferror(f1))
			terror("Error reading a dictionary entry");
        }
	free(W);

	fclose(f1);
	fclose(f2);
	// store PEQ freqs--------
	fprintf(f3,"freqs of words that appear\nTimes\tnWords\n");
	for (i=0;i<PEQ;i++)
		if (freq[i]) fprintf(f3,"%" PRId64 "\t%" PRId64 "\n",i,freq[i]);

	fclose(f3);
	fprintf(stdout,"Num.Words=%" PRIu64 " MaxFreq=%" PRIu64 " TotRepeat=%" PRIu64 " AveragFreq=%f\n",nW,maxF,aveF, (float)aveF/(float)nW);

	exit(0);
}



