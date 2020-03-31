/*
 * Fragmentv3.c
 *
 *  Created on: 23/02/14
 *      Author: jarjonamedina
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

#include "fragmentv3.h"
//#include "error.h"


int writeFragmentsv3 (Fragmentv3* f,char* s,int nf, int xtotal, int ytotal){

	FILE* fs;
	int n;

	if((fs=fopen(s,"wb"))==NULL){
		printf("***ERROR Opening output file");
		exit(-1);
	}

//
	fwrite(&xtotal,sizeof(int),1,fs);
	fwrite(&ytotal,sizeof(int),1,fs);

	for(n=0;n<nf;n++){
		fwrite(&f[n],sizeof(Fragmentv3),1,fs);
#ifdef debug
//		printf("%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%d\t%c\n",f[n].xIni,f[n].yIni,f[n].xFin,f[n].yFin,f[n].length,f[n].ident,f[n].score,f[n].seqX,f[n].seqY,f[n].block,f[n].strand);
#endif 
	}
	fclose(fs);

return 0;
}

//
Fragmentv3* readFragmentsv3(char* s,int* nf,int *xtotal,int *ytotal){

	FILE* fe;

	Fragmentv3* fs;
	int n;

	if((fe=fopen(s,"rb"))==NULL){
		printf("***ERROR Opening input file: %s\n",s);
		exit(-1);
	}


	n=0;
	fread(xtotal,sizeof(int),1,fe);
	fread(ytotal,sizeof(int),1,fe);

#ifdef debug
//	printf("readFragments\nxIni \tyIni \txFin \tyFin \tlength \tident \tscore \tseqX \tseqY \tblock \tstrand\n");
#endif 

	long int longFile;
	fseek(fe,0,SEEK_END);
	longFile=ftell(fe);
	n=(int)(longFile-2*sizeof(int))/sizeof(Fragmentv3);
//	printf("\n ReadFragments Complete\nnum: %d\n",n);
	rewind(fe);

	//Incorporamos fragmento punto y final
	fs=(Fragmentv3*)malloc(sizeof(Fragmentv3)*(n));
	if(fs==NULL){
		printf("****ERROR: Out of memory\n");
		exit(-1);
	}
	
	*nf=n;
	n=0;
	
	fread(xtotal,sizeof(int),1,fe);
	fread(ytotal,sizeof(int),1,fe);

	while(!feof(fe)){
		fread(&fs[n],sizeof(Fragmentv3),1,fe);
		n++;
	}

	fclose(fe);

	return fs;
}




