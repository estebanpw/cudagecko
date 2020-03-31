/* FragHITS  Fragments detection betweem two bioogical seqeunces

 Search for all By-Identitity ungapped fragments with length > MinLen

 This program is an enhancement of AllFrags (that searches the full space)
 In this case, start-Fragments (HITS) has been detected in previous steps and
 this hits are used as seeds to extend the fragment. A backward step is also
 included to extend the fragment in the second opposite direction from the hit
 
 Syntax:

 AllFragsHits SeqX.file SeqY.file HitsFile Out.file Lmin SimThr

 Some new parameters have been included:
 - HitsFile (binary) long Diag/posX/posY
 - Lmin     Minimum length in the case of fixed length and percentage otherwise
 - SimThr   Similarity (identity) Threshold

 - SeqX & SeqY files are in fasta format
 - Out file (binary) will save a fragment with the format specified in structs.h.
 
 oscart@uma.es
 -------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "commonFunctions.h"
#include "structs.h"
#include "comparisonFunctions.h"

#define POINT 4

int HitsControl(char **av);
int FragFromHit(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
		uint64_t n0, struct Sequence *sY,
		uint64_t n1, uint64_t nSeqs1, uint64_t LenOrPer, uint64_t SimTh, int WL,
		int fixedL, char strand);

int main(int ac, char **av) {

	if (ac != 10)
		terror(
				"usage: FragsHits SeqX SeqY Hits.File Outfile Lmin SimThr WL fixedLen(1/0) strand(f/r)");

	if (av[9][0] != 'f' && av[9][0] != 'r')
		terror("strand argument must be 'f' or 'r'");

	HitsControl(av);

	return 0;
}

int HitsControl(char **av) {
	struct Sequence *sX, *sY;
	uint64_t n0, n1, nSeqs0, nSeqs1;
	uint64_t Lmin, SimTh;
	int WL, fixedL, i, j;
	int newFrag;
	uint64_t nFrags = 0, nHitsUsed = 0, nHits = 0;
	int64_t lastOffset, lastDiagonal;
	struct FragFile myF;
	char infoFileName[200], matFileName[200];
	char strand;
	FILE *f, *fH, *fOut, *fInf, *fMat;

	hit h;

	//MAT
	long coverage[1000][100];
	for(i=0; i<1000; i++){
		for(j=0; j<100; j++){
			coverage[i][j]=0;
		}
	}
	//---

	//Initialize values
	Lmin = atoi(av[5]);
	SimTh = atoi(av[6]);
	WL = atoi(av[7]);
	fixedL = atoi(av[8]);
	strand = av[9][0];

	//Open files
	if ((f = fopen(av[1], "rt")) == NULL)
		terror("opening seqX file");
	sX = LeeSeqDB(f, &n0, &nSeqs0, 0);
	fclose(f);

	if ((f = fopen(av[2], "rt")) == NULL)
		terror("opening seqY file");
	sY = LeeSeqDB(f, &n1, &nSeqs1, 0);
	fclose(f);

	if ((fH = fopen(av[3], "rb")) == NULL)
		terror("opening HITS file");

	if ((fOut = fopen(av[4], "wb")) == NULL)
		terror("opening output file");

	//Write sequence lengths
	writeSequenceLength(&n0, fOut);
	writeSequenceLength(&n1, fOut);

	// read Hits
	if(fread(&h, sizeof(hit), 1, fH)!=1){
		terror("Empty hits file");
	}
	lastDiagonal = h.diag;
	lastOffset = h.posX - 1;

	while (!feof(fH)) {
		nHits++;
		if (lastDiagonal != h.diag) {
			//Different diagonal update the variables
			lastDiagonal = h.diag;
			lastOffset = h.posX - 1;
		}
		//We have a casting here because of a funny error
		//where the program was saying that -1 > 0
		if (lastOffset > ((int64_t) h.posX)) {
			//The hit is covered by a previous frag in the same diagonal
			if(fread(&h, sizeof(hit), 1, fH)!=1){
				if(ferror(fH))terror("Error reading from hits file");
			}
			continue;
		}

		nHitsUsed++;
		newFrag = FragFromHit(coverage, &myF, &h, sX, n0, sY, n1, nSeqs1, Lmin, SimTh,
				WL, fixedL, strand);
		if (newFrag) {
			writeFragment(&myF, fOut);
			lastOffset = h.posX + myF.length;
			nFrags++;
		}
		if(fread(&h, sizeof(hit), 1, fH)!=1){
			if(ferror(fH))terror("Error reading from hits file");
		}
	}

	fclose(fH);
	fclose(fOut);

	sprintf(matFileName, "%s.MAT", av[4]);
	if ((fMat = fopen(matFileName, "wt")) == NULL)
		terror("opening MAT file");

	for(i=0; i<1000; i++){
		for(j=0; j<100; j++){
			fprintf(fMat, "%ld\t", coverage[i][j]);
		}
		fprintf(fMat, "\n");
	}

	fclose(fMat);

	// metadata info 
	sprintf(infoFileName, "%s.INF", av[4]);
	if ((fInf = fopen(infoFileName, "wt")) == NULL)
		terror("opening INFO file");

	fprintf(fInf, "All by-Identity Ungapped Fragments (Hits based approach)\n");
	fprintf(fInf, "[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>\n");
	fprintf(fInf, "SeqX filename        : %s\n", av[1]);
	fprintf(fInf, "SeqY filename        : %s\n", av[2]);
	fprintf(fInf, "SeqX name            : %s\n", sX->ident);
	fprintf(fInf, "SeqY name            : %s\n", sY->ident);
	fprintf(fInf, "SeqX length          : %" PRIu64 "\n", n0);
	fprintf(fInf, "SeqY length          : %" PRIu64 "\n", n1);
	fprintf(fInf, "Min.fragment.length  : %" PRIu64 "\n", Lmin);
	fprintf(fInf, "Min.Identity         : %" PRIu64 "\n", SimTh);
	fprintf(fInf, "Tot Hits (seeds)     : %" PRIu64 "\n", nHits);
	fprintf(fInf, "Tot Hits (seeds) used: %" PRIu64 "\n", nHitsUsed);
	fprintf(fInf, "Total fragments      : %" PRIu64 "\n", nFrags);
	fprintf(fInf, "========================================================\n");
	fclose(fInf);

	return nHitsUsed;
}

/**
 * Compute a fragments from one seed point
 * Similarirty thershold and length > mimL
 */
int FragFromHit(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
		uint64_t n0, struct Sequence *sY,
		uint64_t n1, uint64_t nSeqs1, uint64_t Lm, uint64_t SimTh, int WL,
		int fixedL, char strand) {
	int64_t ldiag, ldiag2;
	int64_t xfil, ycol;
	/* for version with backward search */
	int64_t xfil2, ycol2;
	int fragmentLength = 0;
	/* for version Maximum global---*/
	int64_t xfilmax, ycolmax;
	/* for version with backward search */
	int64_t xfilmax2, ycolmax2;
	int nIdentities = 0, maxIdentities = 0;
	char valueX, valueY;
	int fscore = 0, fscoreMax = 0; // full score

	uint64_t minLength =
			(fixedL) ?
					Lm :
					(uint64_t) (min(getSeqLength(sX),
							getSeqLength(sY)) * (Lm / 100.0));

	// Initialize values
	ldiag = min(n0 - H->posX, n1 - H->posY);
	//var to know how much we have until we reach the origin of coordinates
	ldiag2 = min(H->posX, H->posY);
	xfil = H->posX + WL;
	xfil2 = H->posX - 1;
	ycol = H->posY + WL;
	ycol2 = H->posY - 1;
	fragmentLength += WL;
	xfilmax = xfil;
	xfilmax2 = xfil2;
	ycolmax = ycol;
	ycolmax2 = ycol2;
	nIdentities = maxIdentities = WL;
	fscore = POINT * WL;
	fscoreMax = fscore;

	// now, looking for end_frag---
	while (fragmentLength < ldiag) {
		valueX = getValue(sX, xfil);
		valueY = getValue(sY, ycol);
		if (valueX == '*' || valueY == '*') {
			//separator between sequences
			break;
		}

		if (valueX == 'N' || valueY == 'N') {
			fscore -= 1;
		} else {
			if (valueX == valueY) {
				// match
				fscore += POINT;
				nIdentities++;
				if (fscoreMax < fscore) {
					fscoreMax = fscore;
					xfilmax = xfil;
					ycolmax = ycol;
					maxIdentities = nIdentities;
				}
			} else {
				fscore -= POINT;
			}
		}

		xfil++;
		ycol++;
		fragmentLength++;
		if (fscore < 0)
			break;
	}

	/**
	 * Backward search --- Oscar (Sept.2013)
	 **/
	fragmentLength = 0;
	fscore = fscoreMax;
	xfilmax2 = H->posX;
	ycolmax2 = H->posY;
	nIdentities = maxIdentities;
	if (xfil2 >= 0 && ycol2 >= 0)
		while (fragmentLength < ldiag2) {
			valueX = getValue(sX, xfil2);
			valueY = getValue(sY, ycol2);
			if (valueX == '*' || valueY == '*') {
				//separator between sequences
				break;
			}

			if (valueX == 'N' || valueY == 'N') {
				fscore -= 1;
			} else {
				if (valueX == valueY) {
					// matches----
					fscore += POINT;
					nIdentities++;
					if (fscoreMax < fscore) {
						fscoreMax = fscore;
						xfilmax2 = xfil2;
						ycolmax2 = ycol2;
						maxIdentities = nIdentities;
					}
				} else {
					fscore -= POINT;
				}
			}

			xfil2--;
			ycol2--;
			fragmentLength++;
			if (fscore < 0)
				break;
		}

	// Set the values of the FragFile
	myF->diag = H->diag;
	myF->xStart = (uint64_t) xfilmax2 - H->seqX;
	myF->yStart = (uint64_t) ycolmax2 - H->seqY;
	myF->xEnd = (uint64_t) xfilmax - H->seqX;
	myF->yEnd = (uint64_t) ycolmax - H->seqY;;
	myF->length = myF->xEnd - myF->xStart + 1;
	myF->score = fscoreMax;
	myF->ident = maxIdentities;
	myF->similarity = myF->score * 100.0
			/ scoreMax(&sX->datos[myF->xStart], &sY->datos[myF->yStart],
					myF->length, POINT);
	myF->seqX = H->seqX;
	myF->seqY = (strand=='f')? H->seqY : nSeqs1 - H->seqY - 1;
	myF->block = 0;
	myF->strand = strand;

	M[min(myF->length, 999)][(int)myF->similarity]++;

	if (myF->length > minLength && myF->similarity > SimTh)
		return 1;
	else
		return 0;
}

