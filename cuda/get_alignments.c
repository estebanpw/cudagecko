#include "common.h"

int main(int argc, char ** av){
   
    // Check parameters 
    if(argc < 4) {
        fprintf(stderr, "USE: get_alignments csv fasta_x fasta_y [min align length]\n"); 
        exit(-1);
    }
    
    
    // Open files for reading
    FILE * csv = NULL;
    csv = fopen(av[1], "rt"); if(csv == NULL) { fprintf(stderr, "Could not open input csv file\n"); exit(-1);}  

	FILE * fastax = NULL;
	fastax = fopen(av[2], "rt"); if(fastax == NULL) { fprintf(stderr, "Could not open input fasta X file\n"); exit(-1);}

	FILE * fastay = NULL;
	fastay = fopen(av[3], "rt"); if(fastay == NULL) { fprintf(stderr, "Could not open input fasta Y file\n"); exit(-1);}

	uint32_t min_len = 0;
	if(argc == 5) min_len = (uint32_t) atoi(av[4]);

    // Load sequence
    uint32_t l_fastax, l_fastay;
	std::vector<uint64_t> index_x, index_y, index_r;
    char * s_x = load_seq(fastax, &l_fastax, &index_x);
    char * s_y = load_seq(fastay, &l_fastay, &index_y);

	// Generate reverse index	
	for (std::vector<uint64_t>::reverse_iterator it=index_y.rbegin(); it!=index_y.rend(); ++it)
		index_r.push_back(l_fastay - *it);
	

    // Close input sequences since they are already in memory
    fclose(fastax);
	fclose(fastay);

	// Get reverse complement sequence
	char * r_y = reverse_complement_sequence(s_y, l_fastay);


	/*
	// Debug
	FILE * out_x = fopen("x.fasta", "wt");
	FILE * out_y = fopen("y.fasta", "wt");
	FILE * out_yr = fopen("yrev.fasta", "wt");

	fwrite(s_x, 1, l_fastax, out_x); 
	fwrite(s_y, 1, l_fastay, out_y); 
	fwrite(r_y, 1, l_fastay, out_yr); 

	fclose(out_x);
	fclose(out_y);
	fclose(out_yr);
	*/

	// Get alignments
	get_alignments(s_x, s_y, r_y, l_fastax, l_fastay, &index_x, &index_y, &index_r, csv, min_len);
	

	free(s_x);
	free(s_y);
	free(r_y);
}
