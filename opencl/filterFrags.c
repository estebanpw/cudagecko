/*

filterFrags.c

Exports frags to CSV filtering by thresholds

Use: filterFrags file.frags length similarity

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"

int main(int ac,char** av){

	char seqXname[2000] = "SEQ-X";
	char seqYname[2000] = "SEQ-Y";
	
	if(ac<4){
		printf("Use: filterFrags <file.frags> <length> <similarity> [seqXname] [seqYname]\n");
		exit(-1);
	}
	if(ac==6){
		seqXname[0] = '\0';
		strcat(seqXname, av[4]);
		seqYname[0] = '\0';
		strcat(seqYname, av[5]);
	}
	
	// Read fragments
	struct FragFile frag;
	uint64_t xtotal,ytotal;
	//f=readFragmentsv2(av[1],&nf,&xtotal,&ytotal);

	FILE * fFrags = fopen(av[1], "rb");
	if(fFrags == NULL){ fprintf(stderr, "Could not open input file\n"); exit(-1); }

	read_sequence_length(&xtotal, fFrags);
	read_sequence_length(&ytotal, fFrags);
	

	

	uint64_t min_l = (uint64_t) atoi(av[2]);
	double min_sim = (double) atof(av[3]);
	
	/******************************************/
	fprintf(stdout, "All by-Identity Ungapped Fragments (Hits based approach)\n");
	fprintf(stdout, "[Abr.98/Apr.2010/Dec.2011/Jun.2018 -- BITLAB\n");
	fprintf(stdout, "SeqX filename        : %s\n", seqXname);
	fprintf(stdout, "SeqY filename        : %s\n", seqYname);
	fprintf(stdout, "SeqX name            : %s\n", seqXname);
	fprintf(stdout, "SeqY name            : %s\n", seqYname);
	fprintf(stdout, "SeqX length          : %" PRIu64 "\n", xtotal);
	fprintf(stdout, "SeqY length          : %" PRIu64 "\n", ytotal);
	fprintf(stdout, "Min.fragment.length  : %" PRIu64 "\n", min_l);
	fprintf(stdout, "Min.Identity         : %" PRIu64 "\n", (uint64_t) (100*min_sim));
	fprintf(stdout, "Tot Hits (seeds)     : Not in use\n");
	fprintf(stdout, "Tot Hits (seeds) used: Not in use\n");
	fprintf(stdout, "Total fragments      : Not in use\n");
	fprintf(stdout, "========================================================\n");

	
	
	printf("Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");
	
	//double similarity,likeness;
				
	while(!feof(fFrags)){
		if(0 == fread(&frag, sizeof(struct FragFile), 1, fFrags)) { fprintf(stderr, "Read zero frags\n"); exit(-1); }

		
		if(frag.strand=='r'){
			frag.yStart = ytotal - frag.yStart - 1;
			frag.yEnd = ytotal - frag.yEnd - 1;
		}
		
		
		
		if(frag.similarity >= min_sim && (uint64_t)frag.length >= min_l){
			//fprintf(stdout, "Wtf frg ident is : %"PRIu64" %"PRIu64"\n", frag.ident, frag.length);		
			//getchar();
			fprintf(stdout, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%.2f,%.2f,%"PRIu64",%"PRIu64"\n",frag.xStart,frag.yStart,frag.xEnd,frag.yEnd,frag.strand,frag.block,frag.length,frag.score,frag.ident,frag.similarity,(float)(100*frag.ident)/frag.length,frag.seqX,frag.seqY);
		}
				
	}

	fclose(fFrags);
	
	return 0;
}
