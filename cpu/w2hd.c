/* creates a hash table in disk from a set of ordered words
   Syntax: w2hd  wordsSort.In  prefixNameOUT

    wordsSort is a bin file with Word-Pos-Seq
    prefixNameOUT.h2dW  : index of words-Pos-Ocurrences
    prefixNameOUT.h2dP  : positions(Pos+seq)

    Feb.2011: add a new parameter: PrefixSize

    PrefixSize: defines the word-prefix size to be used to identify when two
		words are the "same"

                              ortrelles@uma.es / Dic.2011
    ---------------------------------------------------------*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"

int main(int ac, char** av){

	char fname[1024];
	uint64_t nWords=0;
	FILE* fw, *fOut1, *fOut2;

	wentry we;
	hashentry he;
	location loc;

	if(ac!=3)terror("USE: w2hd  wordsSort.In  prefixNameOUT\n");

	if ((fw = fopen(av[1],"rb"))==NULL) terror("opening IN file");

	sprintf(fname,"%s.d2hW",av[2]);
	if ((fOut1 = fopen(fname,"wb"))==NULL) terror("opening prefix.d2hW file");
	sprintf(fname,"%s.d2hP",av[2]);
	if ((fOut2 = fopen(fname,"wb"))==NULL) terror("opening prefix.d2hP file");

	if(fread(&we,sizeof(wentry),1,fw)!=1){
		terror("empty words file");
	}
	memcpy(&he.w.b[0],&we.w.b[0],8);
	he.pos=0;
	he.num=0;
   
	while(!feof(fw)){
		  loc.pos=we.pos;
		  loc.seq=we.seq;
		  if (wordcmp(&he.w.b[0],&we.w.b[0],32)!=0) {
			 fwrite(&he,sizeof(hashentry),1,fOut1);
			 memcpy(&he.w.b[0],&we.w.b[0],8);
			 he.pos=ftell(fOut2);
			 he.num=0;
		  }

		  fwrite(&loc,sizeof(location),1,fOut2);
		  he.num++;
		  nWords++;

		  if(fread(&we,sizeof(wentry),1,fw)!=sizeof(wentry)){
			if(ferror(fw))terror("error reading words file");
		  }
	}

	fwrite(&he,sizeof(hashentry),1,fOut1);

	fprintf(stdout,"\nw2hd: %s tot words=%" PRIu64 "\n",av[1],nWords);

	fclose(fOut1);
	fclose(fOut2);
	fclose(fw);

	return 0;
}



