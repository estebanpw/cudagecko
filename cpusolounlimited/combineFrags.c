/* reduce 2 fragments files produced in parallel

 Syntax: reduceFrags PrefixFRAGSfiles nChunks fileOUT

 --------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

void writeINFHeader(FILE *fIn, FILE *fOut);
void getNumberOfHits(FILE *infFile, uint64_t *nHits, uint64_t *nHitsUsed, uint64_t *nFrags);

int main(int ac, char** av) {
	uint64_t n1, n2, nHitsForward, nHitsForwardUsed, nHitsReverse, nHitsReverseUsed, nFragsForward, nFragsReverse;
	char infoFileName[200], matFileName[200];
	FILE *fIn, *fInfIn1, *fInfIn2, *fInf, *fOut, *fMat;

	struct FragFile frag;

	long coverage[1000][100], valor, i, j;
	for(i=0; i<1000; i++){
                for(j=0; j<100; j++){
                        coverage[i][j]=0;
                }
        }

	if (ac != 4)
		terror("Syntax: combineFrags forwardStrandFrags reverseStrandFrags fileOUT");

	if ((fIn = fopen(av[1], "rb")) == NULL)
		terror("Input forwardStrandFrags open error");

	sprintf(infoFileName, "%s.INF", av[1]);
	if ((fInfIn1 = fopen(infoFileName, "rt")) == NULL)
		terror("Input forwardStrandFrags.INF open error");

	//MAT
	sprintf(matFileName, "%s.MAT", av[1]);
        if ((fMat = fopen(matFileName, "rt")) == NULL)
                terror("opening forwardStrandFrags.MAT file");
	for (i=0;i<1000;i++){
		for (j=0;j<100;j++){
			if(fscanf(fMat,"%ld\t",&valor)!=1){
				terror("Error reading forwardStrandFrags.MAT file");
			}
			coverage[i][j]+=valor;
		}
	}

	fclose(fMat);

	if ((fOut = fopen(av[3], "wb")) == NULL)
		terror("OUTput file open error");

	// sequences lengths
	readSequenceLength(&n1, fIn);
	readSequenceLength(&n2, fIn);

	writeSequenceLength(&n1, fOut);
	writeSequenceLength(&n2, fOut);

	//First file...
	readFragment(&frag, fIn);

	while (!feof(fIn)) {
		writeFragment(&frag, fOut);
		readFragment(&frag, fIn);
	}

	fclose(fIn);

	if ((fIn = fopen(av[2], "rb")) == NULL)
		terror("Input reverseStrandFrags open error");

	sprintf(infoFileName, "%s.INF", av[2]);
	if ((fInfIn2 = fopen(infoFileName, "rt")) == NULL)
		terror("Input reverseStrandFrags.INF open error");

	//MAT
	sprintf(matFileName, "%s.MAT", av[2]);
	if ((fMat = fopen(matFileName, "rt")) == NULL)
                terror("opening reverseStrandFrags.MAT file");
	for (i=0;i<1000;i++){
		for (j=0;j<100;j++){
			if(fscanf(fMat,"%ld\t",&valor)!=1)
				terror("Error reading reverseStrandFrags.MAT file");
			coverage[i][j]+=valor;
		}
	}

	fclose(fMat);

	readSequenceLength(&n1, fIn);
	readSequenceLength(&n2, fIn);

	//Second file...
	readFragment(&frag, fIn);

	while (!feof(fIn)) {
		writeFragment(&frag, fOut);
		readFragment(&frag, fIn);
	}

	//MAT
	sprintf(matFileName, "%s.MAT", av[3]);
        if ((fMat = fopen(matFileName, "wt")) == NULL)
                terror("opening MAT file");
	for (i=0;i<1000;i++){
		for (j=0;j<100;j++){
			fprintf(fMat,"%ld\t",coverage[i][j]);
		}
		fprintf(fMat, "\n");
	}

	fclose(fMat);

	// metadata info 
	sprintf(infoFileName, "%s.INF", av[3]);
	if ((fInf = fopen(infoFileName, "wt")) == NULL)
		terror("opening INFO file");

	writeINFHeader(fInfIn1, fInf);
	if(fseek(fInfIn1, 0, SEEK_SET) != 0)
		terror("Error rewinding inf file");

	getNumberOfHits(fInfIn1, &nHitsForward, &nHitsForwardUsed, &nFragsForward);
	getNumberOfHits(fInfIn2, &nHitsReverse, &nHitsReverseUsed, &nFragsReverse);

	fprintf(fInf, "Tot Hits (seeds)     : %" PRIu64 "\n", nHitsForward + nHitsReverse);
	fprintf(fInf, "Tot Hits (seeds) used: %" PRIu64 "\n", nHitsForwardUsed + nHitsReverseUsed);
	fprintf(fInf, "Total fragments      : %" PRIu64 "\n", nFragsForward + nFragsReverse);
	fprintf(fInf, "========================================================\n");
	fclose(fInf);

	fclose(fIn);
	fclose(fInfIn1);
	fclose(fInfIn2);

	fclose(fOut);

	return 0;

}

void writeINFHeader(FILE *fIn, FILE *fOut){
	char tmp[1024];
	int i=0;
	while((i < 10) && (fgets(tmp, 1024, fIn)!=NULL)){
		fprintf(fOut, "%s", tmp);		
		i++;
	}
}

void getNumberOfHits(FILE *infFile, uint64_t *nHits, uint64_t *nHitsUsed, uint64_t *nFrags){
	char tmp[1024];
	int i=0;
	while((i < 10) && (fgets(tmp, 1024, infFile)!=NULL)){
		i++;
	}

	if(fgets(tmp, 1024, infFile)==NULL)
		terror("Wrong format in INF file, Tot Hits (seeds) line expected");

	if(sscanf(tmp, "Tot Hits (seeds)     : %" PRIu64 "\n", nHits)!=1)
		terror("Error reading: Tot Hits (seeds)");

	if(fgets(tmp, 1024, infFile)==NULL)
		terror("Wrong format in INF file, Tot Hits (seeds) used line expected");

	if(sscanf(tmp, "Tot Hits (seeds) used: %" PRIu64 "\n", nHitsUsed)!=1)
		terror("Error reading: Tot Hits (seeds) used");

	if(fgets(tmp, 1024, infFile)==NULL)
		terror("Wrong format in INF file, Total fragments line expected");

	if(sscanf(tmp, "Total fragments      : %" PRIu64 "\n", nFrags)!=1)
		terror("Error reading: Total fragments");
}
