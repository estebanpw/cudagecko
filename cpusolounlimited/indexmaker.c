
/********************* 
File		indexmaker.c
Author		Bitlab
Description	Computes a series of statistics and stores them in natural order and sorted order
	
USAGE		<sequences.fasta>	The fasta sequences
			<output.index>		The output name for the two indexes
 * ------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"
#include "comparisonFunctions.h"
#include "commonFunctions.h"


void copyRead(struct rIndex2 *To, struct rIndex2 *From);
void QsortRead(struct rIndex2 *array,uint64_t x, uint64_t y);

int main(int ac, char **av) {

	FILE   *f, *fOut;
	struct rIndex2 R;//, *RR;
	//uint64_t *positions;
	uint64_t  c;
	uint64_t tLen=0, tLmask=0, tN=0, nReads=0;
	uint64_t maxRlen=0, maxRlenMasked=0, maxNonACGT=0;
	uint64_t N=0, i;
	//char tmp[MAXLID];

	if (ac!=3) terror("USE: indexmaker <sequences.fasta> <output.index>");

	
	if ((f=fopen(av[1],"rt"))==NULL) terror("Could not open sequences fasta file");

	if ((fOut=fopen(av[2],"wb"))==NULL) terror("Could not open output file");

    R.Lac=0;
    R.rLen=0;
    c=fgetc(f);
    while (!feof(f)) {
       	while(c!='>') c=(char)fgetc(f); // only first time
		R.pos=ftell(f) - 1;
		i=0;
        c=getc(f);
        if(i==0 && c== ' ') c=getc(f); //Handle sequences with a ' ' after the '>'
        while(i< MAXLID && c!='\n' && c!=' ') {
         	R.id[i++] = c;
          	c=getc(f);
        }
        R.id[i]=0x00;

        while(c!='\n') c=getc(f);
		
		i=0;
		R.rLen=0;
		R.rLmasked=0;
		R.nonACGT=0;
       
		while(c!='>' && !feof(f)) {
			if (c=='N' || c=='n') R.nonACGT++; 
			else if (c>96 && c<123) R.rLmasked++;
		   	c=toupper(c);
		  	if (c>64 && c<91) R.rLen++;
		  	if (c>47 && c<58) R.rLen++; //Handling color space reads
			
			c=(char)fgetc(f);
		}
		R.rNumber = N;
		if (R.rLen > maxRlen) maxRlen=R.rLen;
		if (R.rLmasked > maxRlenMasked) maxRlenMasked = R.rLmasked;
		if (R.nonACGT  > maxNonACGT) maxNonACGT = R.nonACGT;
		
		fwrite(&R,sizeof(struct rIndex2),1,fOut); 

		N++;

		R.Lac +=R.rLen;
		tLen  +=R.rLen;
		tN    +=R.nonACGT;
		tLmask+=R.rLmasked;

		nReads++;

	}

	fclose(f);
	fclose(fOut);

	fprintf(stdout,"....done\n");
	return 0;
}

void copyRead(struct rIndex2 *To, struct rIndex2 *From){

    strcpy(To->id, From->id);
    To->rNumber = From->rNumber;
    To->rLen    = From->rLen;
    To->rLmasked= From->rLmasked;
    To->nonACGT = From->nonACGT;
    To->pos     = From->pos;      //Start position of read
    To->Lac     = From->Lac;      // accumulated length
}


void QsortRead(struct rIndex2 *array,uint64_t x, uint64_t y) {

	struct rIndex2 pivote, aux;
	uint64_t x1, y1;

	copyRead(&pivote, &array[(x+y)/2]);
	x1 = x;
	y1 = y;

 	do{
		while (strcmp(pivote.id, array[x1].id)>0) x1++;
		while (strcmp(pivote.id, array[y1].id)<0) y1--;
    	if (x1 < y1) { 
			copyRead(&aux,&array[x1]);
			copyRead(&array[x1], &array[y1]);
			copyRead(&array[y1], &aux);
			x1++;
			y1--;
    	}else if (x1 == y1) x1++;
   } while (x1 <=y1);

  if (x<y1) QsortRead(array,x,y1);
  if (x1<y) QsortRead(array,x1,y);
 }
