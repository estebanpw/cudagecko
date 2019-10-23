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
	uint64_t n1, n2;
	FILE *fIn, *fOut;

	struct FragFile frag;

	
	if (ac != 4)
		terror("Syntax: combineFrags forwardStrandFrags reverseStrandFrags fileOUT");

	if ((fIn = fopen(av[1], "rb")) == NULL)
		terror("Input forwardStrandFrags open error");



	if ((fOut = fopen(av[3], "wb")) == NULL)
		terror("OUTput file open error");

	// sequences lengths
	readSequenceLength(&n1, fIn);
	readSequenceLength(&n2, fIn);

	writeSequenceLength(&n1, fOut);
	writeSequenceLength(&n2, fOut);

	//First file...

	while (!feof(fIn)) {
		if(0 == fread(&frag, sizeof(struct FragFile), 1, fIn)) fprintf(stdout, "Read 0 items or finished?\n");
		fwrite(&frag, sizeof(struct FragFile), 1, fOut);
	}

	fclose(fIn);

	if ((fIn = fopen(av[2], "rb")) == NULL)
		terror("Input reverseStrandFrags open error");



	readSequenceLength(&n1, fIn);
	readSequenceLength(&n2, fIn);

	while (!feof(fIn)) {
		if(0 == fread(&frag, sizeof(struct FragFile), 1, fIn)) fprintf(stdout, "Read 0 items or finished?\n");
		fwrite(&frag, sizeof(struct FragFile), 1, fOut);
	}

	
	fclose(fIn);

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
