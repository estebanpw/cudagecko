/* DEBUG: list the fragments file
 *
 * Sintax: ./leeFrags fragsFILE [lowUpThresh.file] 
 * 
 * optionally it can include a low/Up filter scheme
 *
 *                           O.Trelles >ortrelles @ uma.es>
 *  -----------------------------------------------May 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

#define SCORE   4  // it depends on the score matrix
#define MAXLENGTH  1000

void loadThresholds(int** low, int** up, char* thresholdFile) {
	int length, lowSimilarity, upSimilarity;
	FILE *fThreshold;

	if ((*low = (int*) calloc(MAXLENGTH, sizeof(int))) == NULL)
		terror("enought memory for LowThreshold array");

	if ((*up = (int*) calloc(MAXLENGTH, sizeof(int))) == NULL)
		terror("enought memory for UpperThreshold array");

	if ((fThreshold = fopen(thresholdFile, "rt")) == NULL)
		terror("open Low-Upper threshold file");

	if(fscanf(fThreshold, "%d\t%d\t%d\n", &length, &lowSimilarity, &upSimilarity)!=3){
		terror("Error reading the threshold file. Wrong format intTABintTABint");
	}
	while (!feof(fThreshold)) {
		if (length >= MAXLENGTH)
			terror("LowUp threshold out of range (maxL)");
		*low[length] = lowSimilarity;
		*up[length] = upSimilarity;
		if(fscanf(fThreshold, "%d\t%d\t%d\n", &length, &lowSimilarity, &upSimilarity)!=3){
			terror("Error reading the threshold file. Wrong format intTABintTABint");
		}
	}
	fclose(fThreshold);
}

int main(int ac, char** av) {
	FILE* fFrags;
	uint64_t n1, n2, nFrags = 0;
	struct FragFile frag;
	int *low, *up, length;
	int similarity = 0, zone;
	int flagT = 0;

	if (ac != 3 && ac != 2)
		terror("USE: ./leeFrags fragsFILE [lowUpThresh.file (Optional)]");
	if (ac == 3)
		flagT = 1;

	// Load Thresholds---------------------------------------
	if (flagT == 1) {
		loadThresholds(&low, &up, av[2]);
	}

	// prepared for multiple files
	if ((fFrags = fopen(av[1], "rb")) == NULL)
		terror("Opening Frags binary file");

	readSequenceLength(&n1, fFrags);
	readSequenceLength(&n2, fFrags);
	fprintf(stderr, "working with fragsFile=%s SeqX=%" PRIu64 " seqY=%" PRIu64 "\n", av[1], n1,
			n2);

	readFragment(&frag, fFrags);
	while (!feof(fFrags)) {
		length = frag.length;
		if (length >= MAXLENGTH)
			length = MAXLENGTH - 1;
		fprintf(stdout,
				"d=%" PRId64 "\tx=%" PRIu64 "\ty=%" PRIu64 "\tL=%" PRIu64 "\tsco=%" PRIu64 "\tseqX=%" PRIu64 "\tseqY=%" PRIu64 "\tSIM=%lf\tident=%" PRIu64 "\tstrand=%c\tblock=%" PRIu64 ,
				frag.diag, frag.xStart, frag.yStart, frag.length, frag.score, frag.seqX, frag.seqY,
				frag.similarity, frag.ident, frag.strand,frag.block);
		if (flagT) {
			if (similarity > up[length])
				zone = 2;
			else if (similarity > low[length])
				zone = 1;
			else
				zone = 0;
			fprintf(stdout, "\tLow=%d\tUp=%d\tZone=%d", low[length], up[length],
					zone);

		}
		fprintf(stdout, "\n");
		nFrags++;
		readFragment(&frag, fFrags);
	}
	fprintf(stdout, "nFrags:%" PRIu64 "\n", nFrags);

	fclose(fFrags);

	return 0;
}
