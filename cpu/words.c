/* words.c
 It is recommended to mask the Low-Complexity regions before using this program

 This program generates a set of 32-mers for the given input sequence.

 Usage: "./words seq.IN words.OUT
 where seq.IN is a plain-text sequence
 words.OUT is a binary file of "wentry" structures with the 32-mers
 -----------------------------------------------------
 oscart@uma.es
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"

#define WORD_SIZE 32
#define BYTES_IN_WORD 8

void shift_word(word * w){
	int i;
	for(i=0;i<BYTES_IN_WORD-1;i++){
		w->b[i]<<=2;
		w->b[i]|=(w->b[i+1]>>6);
	}
	w->b[BYTES_IN_WORD-1]<<=2;
}

void main_FILE(char * inFile, char * outFile){
	FILE *f;
	FILE *f2;
	char c;
	char *seq = NULL;
	uint64_t i = 0, r = 0;

	if ((f=fopen(inFile,"rt"))==NULL){
	perror("opening sequence file");
	}
	if ((f2=fopen(outFile,"wb"))==NULL) {
		terror("opening OUT sequence Words file");
	}
	if ((seq = calloc(READBUF, sizeof(char))) == NULL) {
		terror("not enough memory for read buffer");
	}

	//To force the read
	i = READBUF + 1;

	c = buffered_fgetc(seq, &i, &r, f);
	while (c != '\n') {
		c = buffered_fgetc(seq, &i, &r, f);
	}

	wentry WE;
	WE.seq=0;
	unsigned long index=0;
	unsigned long inEntry=0;
	unsigned long NW=0;
	unsigned long Tot=0;
	unsigned long NoACGT=0;
	unsigned long NoC=0;

	c = buffered_fgetc(seq, &i, &r, f);
	while (!feof(f) || (feof(f) && i < r)){
		if (!isupper(toupper(c))){
			if(c=='>'){
				c = buffered_fgetc(seq, &i, &r, f);
				while (c != '\n')
					c = buffered_fgetc(seq, &i, &r, f);
				WE.seq++;
				inEntry=0;
				index++;
			}
			NoC++;
			c = buffered_fgetc(seq, &i, &r, f);
			continue;
		}
		shift_word(&WE.w);
		switch (c) {
			case 'A': inEntry++; break;
			case 'C':
				WE.w.b[BYTES_IN_WORD-1]|=1;
				inEntry++;
				break;
			case 'G':
				WE.w.b[BYTES_IN_WORD-1]|=2;
				inEntry++;
				break;
			case 'T':
				WE.w.b[BYTES_IN_WORD-1]|=3;
				inEntry++;
				break;
			default :
				inEntry=0; NoACGT++; break;
		}
		index++;
		Tot++;
		if(inEntry>=(unsigned long)WORD_SIZE){
			WE.pos=index-WORD_SIZE;
			NW++;
			fwrite(&WE,sizeof(wentry),1,f2);
		}
		c = buffered_fgetc(seq, &i, &r, f);
	}

	fclose(f);

}

int main(int ac, char** av){
	if(ac!=3){
		terror("USE: words seqFile.IN words.OUT");
	}
	main_FILE(av[1], av[2]);
	return 0;
}

