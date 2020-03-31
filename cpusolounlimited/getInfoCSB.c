/*

getInfoCSB.c

This program takes a fragment file with CSB marked and extract CSB composition

Use: getInfoCSB file.frags fragment_composition(0 no, 1 yes) \n

example:
./getInfo S3chr2R-S2chr2R-L200-S40-K8.csb.frags 1 

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

#include "fragmentv3.h"
#include "fragmentv2.h"

#include "lista.h"
#define  MAX_LEVELS  900000


int main(int ac,char** av){

	
	if(ac<2){
		printf("Use: getInfo file.frags\n");
		exit(-1);
	}
	
	// Read fragments
	struct FragFile* f;
	int nf; // Number of fragments
	uint64_t xtotal,ytotal;
	nf=0;
	f=readFragmentsv2(av[1],&nf,&xtotal,&ytotal);
	
	/******************************************/
	int i;
	
	// Print header. Frags file info
	printf("Total CSB: 0\n");
	printf("========================================================\n");
	
	printf("Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");
	
	double likeness;
				
	for(i=0;i<nf;i++){
		//similarity=100.0 * (((double)f[i].score)/((double)f[i].length*4.0));
		likeness=100.0 * (((double)f[i].ident)/((double)f[i].length));
		
		if(f[i].strand=='r'){
			f[i].yStart = ytotal - f[i].yStart - 1;
			f[i].yEnd = ytotal - f[i].yEnd - 1;
		}

		printf("Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%.2f,%.2f,%"PRIu64",%"PRIu64"\n",(uint64_t)f[i].xStart,(uint64_t)f[i].yStart,(uint64_t)f[i].xEnd,(uint64_t)f[i].yEnd,f[i].strand,(uint64_t)f[i].block,(uint64_t)f[i].length,(uint64_t)f[i].score,(uint64_t)f[i].ident,f[i].similarity,likeness,(uint64_t)f[i].seqX,(uint64_t)f[i].seqY);
				
	}
	
	return 0;
}
