/*
 * fragmentv2.c
 *
 *  Created on: 17/07/2013
 *      Author: jarjonamedina
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>


#include <inttypes.h>
#include "structs.h"
#include "comparisonFunctions.h"

#include "fragmentv2.h"

//#include "error.h"


int writeFragments (struct FragFile* f,char* s,int nf, uint64_t xtotal, uint64_t ytotal){
/*
	FILE* fs;
	int n;

	if((fs=fopen(s,"wb"))==NULL){
		printf("***ERROR Opening output file");
		exit(-1);
	}

//
//	fwrite(&xtotal,sizeof(int),1,fs);
//	fwrite(&ytotal,sizeof(int),1,fs);
	
	fwrite(&xtotal,sizeof(uint64_t),1,fs);
	fwrite(&ytotal,sizeof(uint64_t),1,fs);

	for(n=0;n<nf;n++){
		fwrite(&f[n],sizeof(Fragment),1,fs);
//		printf("%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%d\t%c\n",f[n].xIni,f[n].yIni,f[n].xFin,f[n].yFin,f[n].length,f[n].ident,f[n].score,f[n].seqX,f[n].seqY,f[n].block,f[n].strand);
	}
	fclose(fs);
*/
return 0;

}

//
struct FragFile* readFragmentsv2(char* s,int* nf,uint64_t *xtotal,uint64_t *ytotal){
//Fragment* readFragments(char* s,int* nf,int *xtotal,int *ytotal){

	FILE* fe;

	struct FragFile* fs;
	int n;

	if((fe=fopen(s,"rb"))==NULL){
		printf("***ERROR Opening input file");
		exit(-1);
	}

/*
	if((fsalida=fopen("fragmentosv2.txt","w"))==NULL){
		printf("Fallo al leer Fragmentos");
		exit(-1);
	}

*/
	n=0;
//	fread(xtotal,sizeof( uint64_t),1,fe);
//	fread(ytotal,sizeof( uint64_t),1,fe);

	readSequenceLength(xtotal, fe);
	readSequenceLength(ytotal, fe);

	long int longFile;
	fseek(fe,0,SEEK_END);
	longFile=ftell(fe);
	n=(int)(longFile-2*sizeof(uint64_t))/sizeof(struct FragFile);
//	n=(int)(longFile-2*sizeof(int))/sizeof(Fragment);
	//printf("\n ReadFragments Complete\nnum: %d\n",n);

	rewind(fe);

	fs=(struct FragFile*)malloc(sizeof(struct FragFile)*(n+1));
	if(fs==NULL){
		printf("****ERROR: Out of memory\n");
		exit(-1);
	}
	//*nf=n-1;
	*nf=n;
	n=0;
	
//	fread(xtotal,sizeof( uint64_t),1,fe);
//	fread(ytotal,sizeof( uint64_t),1,fe);


	readSequenceLength(xtotal, fe);
	readSequenceLength(ytotal, fe);
	

//	while(!feof(fe)){
	for(n=0;n<*nf;n++){
		readFragment(&fs[n], fe);
		//fprintf(stdout,"d=%" PRId64 "\tx=%" PRIu64 "\ty=%" PRIu64 "\tL=%" PRIu64 "\tsco=%" PRIu64 "\tseqX=%" PRIu64 "\tseqY=%" PRIu64 "\tSIM=%f" "\tstrand=%c\n",fs[n].diag, fs[n].xStart, fs[n].yStart, fs[n].length, fs[n].score, fs[n].seqX, fs[n].seqY,fs[n].similarity, fs[n].strand);
//		fprintf(stdout,"%" PRId64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%c" "\t%" PRIu64 "\n",fs[n].xStart, fs[n].yStart, fs[n].xEnd, fs[n].yEnd, fs[n].length, fs[n].strand, fs[n].ident);
		
//		n++;
	}
	

	
	fclose(fe);
//	fclose(fsalida);
	return fs;
}




