/* hits : determine HITS between two sequences based on their 
 dictionaries from disk
 Syntax: hits prefixNameX prefixNameY Outfile

 prefixNameX & Y :refers to *.d2hW : index of words-Pos-Ocurrences
 and *.d2hP : words positions
 Outfile is in the form of [Diag][posX][posY]

 Jan.2012: define an I-O buffer to reduce disk activity
 Feb.2012: - define a function to compare words instead of chars
 - use static buffer for positions
 - New parameter: FreqThre (frequency threshold) to avoid high
 frequency words (word.freq> FreqThr are skipped)

 Feb.6   : new parameter (same meaning as in w2hd): PrefixSize

 May 22  : for long repetitions some of them will be skipped (the step is
 computed as NREP / MaxREP

 ortrelles@uma.es / Dic.2011
 ---------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "comparisonFunctions.h"

#define MAXBUF 10000000
#define MaxREP 10000

int main(int ac, char** av) {

	char fname[1024];
	int hitsInBuf = 0;
	int i, j;
	int comp;
	int firstMatch = 0, endMatch = 0;
	uint64_t freqThr;
	uint64_t nHits = 0, wordMatches = 0;
	int wSize;
	int stepX, stepY;
	long offsetWordBefore, offsetWordAfter, offsetLocationBefore, offsetLocationAfter;
	FILE *fX1, *fX2, *fY1, *fY2, *fOut;

	location *posX = NULL, *posY = NULL;
	hashentry heX, heY;
	hit *hBuf;

	if (ac != 6)
		terror("USE: hits prefixNameX prefixNameY  Outfile FreqThre PrefixSize");
	freqThr = (uint64_t) atoi(av[4]); // IGNORED
	wSize = atoi(av[5]);

	// I-O buffer
	if ((hBuf = (hit*) calloc(sizeof(hit), MAXBUF)) == NULL)
		terror("HITS: memory for I-O buffer");
	// word positions buffers
	if ((posX = (location*) calloc(MAXBUF, sizeof(location))) == NULL)
		terror("memory for posX buffer.. using MAXBUF=10MB");
	if ((posY = (location*) calloc(MAXBUF, sizeof(location))) == NULL)
		terror("memory for posY buffer.. using MAXBUF=10MB");

	// Sequence X files
	sprintf(fname, "%s.d2hW", av[1]);
	if ((fX1 = fopen(fname, "rb")) == NULL)
		terror("opening seqX.d2hW file");
	sprintf(fname, "%s.d2hP", av[1]);
	if ((fX2 = fopen(fname, "rb")) == NULL)
		terror("opening seqX.d2hP file");

	// Sequence Y files
	sprintf(fname, "%s.d2hW", av[2]);
	if ((fY1 = fopen(fname, "rb")) == NULL)
		terror("opening seqY.d2hW file");
	sprintf(fname, "%s.d2hP", av[2]);
	if ((fY2 = fopen(fname, "rb")) == NULL)
		terror("opening seqY.d2hP file");

	// OUT file
	if ((fOut = fopen(av[3], "wb")) == NULL)
		terror("opening OUT file");

	// kick-off
	if (readHashEntry(&heX, fX1) == -1)
		terror("no hash (1)");
	if (readHashEntry(&heY, fY1) == -1)
		terror("no hash (2)");

	while (!feof(fX1) && !feof(fY1)) {

		comp = wordcmp(&heX.w.b[0], &heY.w.b[0], wSize);
		if (comp < 0) {
			readHashEntry(&heX, fX1);
			//Save position of first missmatch after matches and rewind
			if(firstMatch){
				offsetWordAfter = ftell(fY1) - sizeof(hashentry);
				offsetLocationAfter = ftell(fY2);
				fseek(fY1, offsetWordBefore, SEEK_SET);
				readHashEntry(&heY, fY1);
				fseek(fY2, offsetLocationBefore, SEEK_SET);
				firstMatch = 0;
				endMatch = 1;
			}
			continue;
		}
		if (comp > 0) {
			//No more matches, go to the next word
			if(endMatch){
				fseek(fY1, offsetWordAfter, SEEK_SET);
				fseek(fY2, offsetLocationAfter, SEEK_SET);
				endMatch = 0;
			}
			readHashEntry(&heY, fY1);
			continue;
		}

		wordMatches++;

		// Load word position for seqX
		if(!firstMatch)loadWordOcurrences(heX, &posX, &fX2);

		// Saving the offset of the first match
		if(wSize < 32 && !firstMatch){
			offsetWordBefore = ftell(fY1) - sizeof(hashentry);
			offsetLocationBefore = ftell(fY2);
			firstMatch = 1;
		}

		// Load word position for seqY
		loadWordOcurrences(heY, &posY, &fY2);

		

		stepX = 1;
		stepY = 1;

		for (i = 0; i < heX.num; i += stepX)
			for (j = 0; j < heY.num; j += stepY) {
				hBuf[hitsInBuf].diag = posX[i].pos - posY[j].pos;
				hBuf[hitsInBuf].posX = posX[i].pos;
				hBuf[hitsInBuf].seqX = posX[i].seq;
				hBuf[hitsInBuf].posY = posY[j].pos;
				hBuf[hitsInBuf].seqY = posY[j].seq;

				hitsInBuf++;
				if (hitsInBuf == MAXBUF - 1) {
					fwrite(hBuf, sizeof(hit), hitsInBuf, fOut);
					hitsInBuf = 0;
				}
			}

		nHits += ((heX.num / stepX) * (heY.num / stepY));

		if(!firstMatch)readHashEntry(&heX, fX1);
		readHashEntry(&heY, fY1);

	}

	//Closing dictionary files
	fclose(fX1);
	fclose(fY1);
	fclose(fX2);
	fclose(fY2);

	//Checking if there is something still at the buffer
	if (hitsInBuf != 0) {
		fwrite(hBuf, sizeof(hit), hitsInBuf, fOut);
	}
	fclose(fOut);

	fprintf(stdout, "HITS: matches=%" PRIu64 " Tot HITS=%" PRIu64 "\n", wordMatches, nHits);

	return 0;
}

