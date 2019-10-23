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
	
	double similarity,likeness;
				
	for(i=0;i<nf;i++){
		similarity=(((double)f[i].score)/((double)f[i].length*4.0));
		likeness=(((double)f[i].ident)/((double)f[i].length));
		
		if(f[i].strand=='r'){
			f[i].yStart = ytotal - f[i].yStart - 1;
			f[i].yEnd = ytotal - f[i].yEnd - 1;
		}

		printf("Frag,%d,%d,%d,%d,%c,%d,%d,%d,%d,%.2f,%.2f,%d,%d\n",(int)f[i].xStart,(int)f[i].yStart,(int)f[i].xEnd,(int)f[i].yEnd,f[i].strand,(int)f[i].block,(int)f[i].length,(int)f[i].score,(int)f[i].ident,similarity,likeness,(int)f[i].seqX,(int)f[i].seqY);
				
	}
	
	return 0;
}
