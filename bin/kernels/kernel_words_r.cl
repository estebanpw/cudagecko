#define BYTES_IN_WORD 8

// Parameters for GPU words
typedef struct parameters{
    ulong kmer_size;
	ulong seq_length;
	ulong t_work_items;
	ulong kmers_per_work_item;
    ulong offset;
} Parameters;

//Struct for GPU words
typedef struct {
    unsigned char b[8];     // The actual word
    ulong pos;           // The position in the sequence file
} wordGPU;

void shift_word(unsigned char * b) {
    int i;
    for (i = 0; i < BYTES_IN_WORD - 1; i++) {
        b[i] <<= 2;
        b[i] |= (b[i + 1] >> 6);
    }
    b[BYTES_IN_WORD - 1] <<= 2;
}

void shift_word_right(unsigned char * b) {
    int i;
    for (i = BYTES_IN_WORD - 1; i > 0; i--) {
        b[i] >>= 2;
        b[i] |= (b[i - 1] << 6);
    }
    b[i] >>= 2;
}


__kernel void kernel_words(__global wordGPU * words, __global Parameters * params, __global const char * sequence) {
 

    // Get the index of the current element to be processed
	ulong global_id = get_global_id(0);
	ulong kmers_in_work_item = params->kmers_per_work_item;
	ulong t_work_items = params->t_work_items;
	ulong offset = params->offset;
    ulong kmer_size = params->kmer_size;
	ulong j, k;
	
	// Until reaching end of sequence
	for(j=0; j<kmers_in_work_item; j++){
		
		// Coalescent
		ulong pos = global_id + (j * t_work_items);

		// Non coalescent (Naive approach)
		//ulong pos = (global_id * kmers_in_work_item) + j; 
		
		// Completely not coalescent
		//uint seed = global_id;
		//uint t = seed ^ (seed << 11);  
		//ulong pos = (ulong) ((local_id ^ (local_id >> 19) ^ (t ^ (t >> 8))) % params->seq_length);
		
        unsigned char bad = 0;

		for(k=0; k<kmer_size; k++){

            //shift_word(&words[pos].b);
            shift_word_right(&words[pos].b);
            switch(sequence[pos+k]){
                case 'A':
                    words[pos].b[0] |= 192;
                    break;
                case 'C':
                    words[pos].b[0] |= 128;
                    break;
                case 'G':
                    words[pos].b[0] |= 64;
                    break;
                case 'T':
                    break;
                default:
                    bad = 1;
                    break;
                    
            }
        }

		

		if(bad == 0){
			// seq will be added afterwards
            words[pos].pos = pos;
		}else{
            for(i=0; i<8; i++) words[pos].b[i] = 0xFF; //Put at max so it easy to dispose off
            words[pos].pos = 0xFFFFFFFFFFFFFFFF;
        }
	}

}