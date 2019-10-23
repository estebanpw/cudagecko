/* makefasta.c
 Cuts out dna region

 -----------------------------------------------------
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>


int main(int ac, char** av){
	if(ac!=4){
		terror("USE: makefasta guide seqIn ratio num");
	}
	
	
	FILE * fasta_in = NULL;
	FILE * out = NULL;

	uint64_t guide = (uint64_t) atoi(av[1]);

	fasta_in = fopen(av[2], "rt");
	if(fasta_in == NULL) terror("Could not open fasta file");

	uint64_t ratio = (uint64_t) atoi(av[3]);

	char outname[128] = "temp/";
	strcat(outname, av[4]);
	strcat(outname, ".fasta");
	out = fopen(outname, "wt");
	if(out == NULL) terror("Could not open output file");

	fseek(fasta_in, guide*ratio, SEEK_SET);


	
	fclose(fasta_in);
	fclose(out);
	
	
	return 0;
}

