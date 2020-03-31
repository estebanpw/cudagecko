/* 
 *
 * Sintax: ./frags2text fragsFILE.frags fastaX fastaY fragsFILE.txt
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

#define TAB_INSERT 70
#define READING_FRAG_BUFFER 10000

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

void csv_frag_to_struct_frag(char * l, struct FragFile * f){

	//0		1		2		3		4		5			6		7		8		9		10			11		12		13
	//Type	xStart	yStart	xEnd	yEnd	strand(f/r)	block	length	score	ident	similarity	%ident	SeqX	SeqY
	//Frag	3147493	3006054	3154663	2998884	f	0	7171	20868	6194	72.75	0.86	0	0

	
	float bin;

	sscanf(l, "%*s %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %c %"PRId64" %"PRIu64" %"PRIu64" %"PRIu64" %f %f %"PRIu64" %"PRIu64, &f->xStart, &f->yStart, &f->xEnd, &f->yEnd, &f->strand, &f->block, &f->length, &f->score, &f->ident, &f->similarity, &bin, &f->seqX, &f->seqY);

	//printf("Read %d items\n", items);

}

struct rIndex2 * loadReadsIndex(char * filename, uint64_t * nReads){
	struct rIndex2 * RR;
	uint64_t nR=0,i;
	FILE *f;
	uint64_t fsize;

    if ((f=fopen(filename,"rb"))==NULL) terror("Could not open index input file");

	fseeko(f, 0, SEEK_END);
	fsize = ftello(f);
	rewind(f);
	nR = fsize/sizeof(struct rIndex2);

    if ((RR =(struct rIndex2*) calloc(nR,sizeof(struct rIndex2)))==NULL) terror("Could not allocate index");

	for (i=0; i<nR; i++){
		if(0 == fread(&RR[i],sizeof(struct rIndex2),1,f)) break;
	}
	fclose(f);
	(*nReads) = nR;
	return RR;
}

void write_headers(FILE * f1, FILE * f2, FILE * fO, uint64_t pos1, uint64_t pos2, struct FragFile * f){
	fseek(f1, pos1, SEEK_SET);
	char c = fgetc(f1);
	while(c != '\n'){ fprintf(fO, "%c", c); c = fgetc(f1); }
	fprintf(fO, " ALIGNED WITH ");
	fseek(f2, pos2, SEEK_SET);
	c = fgetc(f2);
	while(c != '\n'){ fprintf(fO, "%c", c); c = fgetc(f2); }
	fprintf(fO, " LENGTH: %"PRIu64" IDENT: %"PRIu64" STRAND: %c @(%"PRIu64", %"PRIu64")\n", f->length, f->ident, f->strand, f->xStart, f->yStart);
}

void get_both_seqs(char * fastaX, char * fastaY, uint64_t iniX, uint64_t finX, uint64_t iniY, uint64_t finY, uint64_t posX, uint64_t posY, uint64_t LacX, uint64_t LacY, uint64_t lX, uint64_t lY, FILE * fO){
	char copyX[TAB_INSERT+1];
	char copyY[TAB_INSERT+1];
	memset(copyX, 0x0, TAB_INSERT+1);
	memset(copyY, 0x0, TAB_INSERT+1);
	uint64_t atcopyX = 0;
	uint64_t atcopyY = 0;
	uint64_t i;

	uint64_t pos_used_x, pos_used_y;
	
	// The X one
	//fseek(fastaX, posX, SEEK_SET);
	pos_used_x = posX;
	char cX = fastaX[pos_used_x++];
	while(cX != '\n'){ cX = fastaX[pos_used_x++]; }
	uint64_t gpX = iniX - LacX;
	uint64_t currX = 0, tab_counter;
	while(currX < gpX){
		cX = fastaX[pos_used_x++];
		if(cX != '\n') ++currX;
	}

	// the other Y
	//fseek(fastaY, posY, SEEK_SET);
	pos_used_y = posY;
	char cY = fastaY[pos_used_y++];
	while(cY != '\n'){ cY = fastaY[pos_used_y++]; }
	uint64_t gpY = iniY - LacY;
	uint64_t currY = 0;
	while(currY < gpY){
		cY = fastaY[pos_used_y++];
		if(cY != '\n') ++currY;
	}


	

	// Reached the region to show
	currX = 0;
	currY = 0;
	cX = fastaX[pos_used_x++];
	cY = fastaY[pos_used_y++];
	//fprintf(fO, "\t");
	tab_counter = 0;
	while(currX < lX && currY < lY){
		
		if(cX == 'A' || cX == 'C' || cX == 'G' || cX == 'T' || cX == 'N'){ copyX[atcopyX++] = cX; ++currX; }
		cX = fastaX[pos_used_x++];
		while(cX != 'A' && cX != 'C' && cX != 'G' && cX != 'T' && cX != 'N') cX = fastaX[pos_used_x++];


		if(cY == 'A' || cY == 'C' || cY == 'G' || cY == 'T' || cY == 'N'){ copyY[atcopyY++] = cY; ++currY; }
		cY = fastaY[pos_used_y++];
		while(cY != 'A' && cY != 'C' && cY != 'G' && cY != 'T' && cY != 'N') cY = fastaY[pos_used_y++];

		while(currX > currY){
			if(cY == 'A' || cY == 'C' || cY == 'G' || cY == 'T' || cY == 'N'){ copyY[atcopyY++] = cY; ++currY; }
			cY = fastaY[pos_used_y++];
			while(cY != 'A' && cY != 'C' && cY != 'G' && cY != 'T' && cY != 'N') cY = fastaY[pos_used_y++];
		}
		while(currX < currY){
			if(cX == 'A' || cX == 'C' || cX == 'G' || cX == 'T' || cX == 'N'){ copyX[atcopyX++] = cX; ++currX; }
			cX = fastaX[pos_used_x++];
			while(cX != 'A' && cX != 'C' && cX != 'G' && cX != 'T' && cX != 'N') cX = fastaX[pos_used_x++];
		}

		++tab_counter;
		


		if(tab_counter >= TAB_INSERT && currX == currY){ 
			
			copyX[TAB_INSERT] = '\0';
			copyY[TAB_INSERT] = '\0';
			fprintf(fO, "X:\t%.*s\n\t", TAB_INSERT, copyX);
			for(i=0; i<TAB_INSERT; i++){
				if(copyX[i] == copyY[i]) fprintf(fO, "|"); else fprintf(fO, " ");
			}
			fprintf(fO, "\nY:\t%.*s\n\n", TAB_INSERT, copyY);
			tab_counter = 0;
			atcopyX = 0; atcopyY = 0;
		}
	}
	if(atcopyX > 0){
		copyX[atcopyX] = '\0';
		copyY[atcopyY] = '\0';
		fprintf(fO, "X:\t%.*s\n\t", (int)atcopyX, copyX);
		for(i=0; i<atcopyX; i++){
			if(copyX[i] == copyY[i]) fprintf(fO, "|"); else fprintf(fO, " ");
		}
		fprintf(fO, "\nY:\t%.*s\n\n", (int)atcopyY, copyY);
		atcopyX = 0; atcopyY = 0;
	} 
	fprintf(fO, "\n");
}


void get_seq_from_to(FILE * fasta, FILE * output, uint64_t ini, uint64_t fin, uint64_t pos, uint64_t Lac, uint64_t seqNum, uint64_t l, FILE * fO){
	fseek(fasta, pos, SEEK_SET);
	char c = fgetc(fasta);
	while(c != '\n'){ c = fgetc(fasta); }
	uint64_t gp = ini - Lac;
	uint64_t curr = 0, tab_counter;
	while(curr < gp){
		c = fgetc(fasta);
		if(c != '\n') ++curr;
	}
	// Reached the region to show
	curr = 0;
	c = fgetc(fasta);
	fprintf(fO, "\t");
	tab_counter = 0;
	while(curr < l){
		if(c != '\n') fprintf(fO, "%c", c);
		c = fgetc(fasta);
		if(c != '\n' || feof(fasta)){ ++curr; ++tab_counter; }
		if(tab_counter == TAB_INSERT){ fprintf(fO, "\n\t"); tab_counter = 0;}
	}
	fprintf(fO, "\n");
}


void get_seq_from_to_rev(FILE * fasta, FILE * output, uint64_t ini, uint64_t fin, uint64_t pos, uint64_t Lac, uint64_t seqNum, uint64_t l, FILE * fO){
	fseek(fasta, pos, SEEK_SET);
	char c = fgetc(fasta);
	while(c != '\n'){ c = fgetc(fasta); }
	uint64_t gp = ini - Lac;
	uint64_t curr = 0, tab_counter;
	while(curr < gp){
		c = fgetc(fasta);
		if(c != '\n') ++curr;
	}
	// Reached the region to show
	curr = 0;
	c = fgetc(fasta);
	fprintf(fO, "\t");
	tab_counter = 0;
	while(curr < l){
		if(c != '\n') fprintf(fO, "%c", c);
		c = fgetc(fasta);
		if(c != '\n' || feof(fasta)){ ++curr; ++tab_counter; }
		if(tab_counter == TAB_INSERT){ fprintf(fO, "\n\t"); tab_counter = 0;}
	}
	fprintf(fO, "\n");
}



int main(int ac, char** av) {
	FILE* fFrags;
	struct FragFile frag;

	//fprintf(stdout, "[WARNING] - Remember that if using a CSV make sure that the reverse y coordinates are transformed\n");

	if (ac != 10)
		terror("USE: ./frags2text fragsFILE.csv fastaX fastaY fastaYrev indexX indexY indexYrev fragsFILE.txt borderSize");

	// prepared for multiple files
	if ((fFrags = fopen(av[1], "rt")) == NULL)
		terror("Opening Frags csv file");

	// Open fastas

	FILE * fX = NULL, * fY = NULL, * fYrev = NULL, * fO = NULL;
	fX = fopen(av[2], "rt");
	if(fX == NULL) terror("Could not open fasta X file");
	fY = fopen(av[3], "rt");
	if(fY == NULL) terror("Could not open fasta Y file");
	fYrev = fopen(av[4], "rt");
	if(fYrev == NULL) terror("Could not open fasta Y-rev file");


	// Get file lengths
    fseek(fX, 0, SEEK_END);
    uint64_t aprox_lenX = ftell(fX);
    rewind(fX);
	char * strfastaX = (char *) malloc(aprox_lenX*sizeof(char));
	fseek(fY, 0, SEEK_END);
    uint64_t aprox_lenY = ftell(fY);
    rewind(fY);
	char * strfastaY = (char *) malloc(aprox_lenY*sizeof(char));
	fseek(fYrev, 0, SEEK_END);
    uint64_t aprox_lenYrev = ftell(fYrev);
    rewind(fYrev);
	char * strfastaYrev = (char *) malloc(aprox_lenYrev*sizeof(char));

	if(strfastaX == NULL || strfastaY == NULL || strfastaYrev == NULL) terror("Could not allocate string sequences");

	if(aprox_lenX != fread(strfastaX, sizeof(char), aprox_lenX, fX)) terror("Read wrong number of chars at X sequence");
	if(aprox_lenY != fread(strfastaY, sizeof(char), aprox_lenY, fY)) terror("Read wrong number of chars at Y sequence");
	if(aprox_lenYrev != fread(strfastaYrev, sizeof(char), aprox_lenYrev, fYrev)) terror("Read wrong number of chars at Y reversed sequence");
	

	struct rIndex2 * RI_X, * RI_Y, * RI_Yrev;

	
	uint64_t nReads_X, nReads_Y;
	RI_X = loadReadsIndex(av[5], &nReads_X);
	RI_Y = loadReadsIndex(av[6], &nReads_Y);
	RI_Yrev = loadReadsIndex(av[7], &nReads_Y);

	fO = fopen(av[8], "wt");
	if(fO == NULL) terror("Could not open output alignments file");

	uint64_t border_size = (uint64_t) atoi(av[9]);


	//readSequenceLength(&n1, fFrags);
	//readSequenceLength(&n2, fFrags);
	// Skip first lines
	char buffer[READING_FRAG_BUFFER];
	// Lengths are in lines 6 and 7
	//SeqX length          : 249250621
	int64_t xlen, ylen;
	int skipper = 0; while(skipper < 18){ 
		fgets(buffer, READING_FRAG_BUFFER, fFrags); 

		if(skipper == 6){
			sscanf(buffer, "%*s %*s : %"PRId64, &xlen);
		}
		
		if(skipper == 7){
			sscanf(buffer, "%*s %*s : %"PRId64, &ylen);
		}
		skipper++; 
	}

	fprintf(stdout, "Sequence lengths: (%"PRId64", %"PRId64")\n", xlen, ylen);
	

	//readFragment(&frag, fFrags);
	
	
	csv_frag_to_struct_frag(buffer, &frag);

	// RI_X 	is the forward index for fasta file X
	// RI_Y 	is the forward index for fasta file Y
	// RI_Yrev 	is the reverse index for fasta file Y

	int exit = 0;

	while (!feof(fFrags) && exit == 0) {
		
		//RI[id].pos is position in file of >
		//RI[id].Lac is the sum of the reads length prior

		if(frag.strand == 'f'){
			write_headers(fX, fY, fO, RI_X[frag.seqX].pos, RI_Y[frag.seqY].pos, &frag);
		}else{
			write_headers(fX, fYrev, fO, RI_X[frag.seqX].pos, RI_Yrev[(nReads_Y - frag.seqY) - 1].pos, &frag);
		}
		

		//get_seq_from_to(fX, fO, frag.xStart, frag.xEnd, RI_X[frag.seqX].pos, RI_X[frag.seqX].Lac, frag.seqX, frag.length, fO);



		if(frag.strand == 'f'){

			int64_t leftx = MAX(0, (int64_t)frag.xStart - (int64_t)border_size);
			int64_t rightx = MIN(xlen, (int64_t)frag.xEnd + (int64_t)border_size);
			int64_t lefty = MAX(0, (int64_t)frag.yStart - (int64_t)border_size);
			int64_t righty = MIN(ylen, (int64_t)frag.yEnd + (int64_t) border_size);

			int64_t border_left_size = MIN(frag.xStart - leftx, frag.yStart - lefty);
			int64_t border_right_size = MIN(rightx - frag.xEnd, righty - frag.yEnd);

			get_both_seqs(strfastaX, strfastaY, frag.xStart - border_left_size, frag.xEnd + border_right_size, frag.yStart - border_left_size, frag.yEnd + border_right_size, RI_X[frag.seqX].pos, RI_Y[frag.seqY].pos, RI_X[frag.seqX].Lac, RI_Y[frag.seqY].Lac, frag.length + border_left_size + border_right_size, frag.length + border_left_size + border_right_size, fO);
			//get_seq_from_to(fY, fO, frag.yStart, frag.yEnd, RI_Y[frag.seqY].pos, RI_Y[frag.seqY].Lac, frag.seqY, frag.length, fO);
		}else{
			uint64_t seqYnew;
			seqYnew = (nReads_Y - frag.seqY) - 1;

			frag.yStart = ylen - frag.yStart - 1;
			frag.yEnd = ylen - frag.yEnd - 1;

			int64_t leftx = MAX(0, (int64_t)frag.xStart - (int64_t)border_size);
            int64_t rightx = MIN(xlen, (int64_t)frag.xEnd + (int64_t)border_size);
            int64_t lefty = MAX(0, (int64_t)frag.yStart - (int64_t)border_size);
            int64_t righty = MIN(ylen, (int64_t)frag.yEnd + (int64_t) border_size);

            int64_t border_left_size = MIN(frag.xStart - leftx, frag.yStart - lefty);
            int64_t border_right_size = MIN(rightx - frag.xEnd, righty - frag.yEnd);



			get_both_seqs(strfastaX, strfastaYrev,  frag.xStart - border_left_size, frag.xEnd + border_right_size, frag.yStart - border_left_size, frag.yEnd + border_right_size, RI_X[frag.seqX].pos, RI_Yrev[seqYnew].pos, RI_X[frag.seqX].Lac, RI_Yrev[seqYnew].Lac, frag.length + border_left_size + border_right_size, frag.length + border_left_size + border_right_size, fO);
			//get_seq_from_to_rev(fYrev, fO, frag.yStart, frag.yEnd, RI_Yrev[seqYnew].pos, RI_Yrev[seqYnew].Lac, seqYnew, frag.length, fO);
		}

		//readFragment(&frag, fFrags);
		if(!feof(fFrags)){
			fgets(buffer, READING_FRAG_BUFFER, fFrags);
			csv_frag_to_struct_frag(buffer, &frag);

		
		}else{
			exit = 1;
		}
		//printf("line read: %s\n", buffer);
		//printf("HI im a frag: %"PRIu64", %"PRIu64" - %"PRIu64", %"PRIu64"\n", frag.xStart, frag.xEnd, frag.yStart, frag.yEnd);
		//getchar();
	}

	fclose(fFrags);
	fclose(fX);
	fclose(fY);
	fclose(fO);

	free(RI_X);
	free(RI_Y);
	free(RI_Yrev);

	free(strfastaX);
	free(strfastaY);
	free(strfastaYrev);

	return 0;
}
