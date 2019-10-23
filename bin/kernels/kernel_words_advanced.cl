
__constant ulong pow4[33]={1L, 4L, 16L, 64L, 256L, 1024L, 4096L, 16384L, 65536L, 
    262144L, 1048576L, 4194304L, 16777216L, 67108864L, 268435456L, 1073741824L, 4294967296L, 
    17179869184L, 68719476736L, 274877906944L, 1099511627776L, 4398046511104L, 17592186044416L, 
    70368744177664L, 281474976710656L, 1125899906842624L, 4503599627370496L, 18014398509481984L, 
    72057594037927936L, 288230376151711744L, 1152921504606846976L, 4611686018427387904L};

//Struct for GPU words
typedef struct {
    ulong h;        // The actual word
    ulong pos;      // The position in the sequence file
} word_hash_GPU;

typedef struct {
    ulong offset;
    ulong end;
    ulong kmer_size;
    ulong kmers_in_work_item;
    ulong t_work_items;
} param_words_advanced;

__kernel void kernel_words_advanced(__global word_hash_GPU * words, __global param_words_advanced * params, __global char * sequence) {
 
    
    // Get the index of the current element to be processed
	ulong global_id = get_global_id(0);
	//ulong local_id = get_local_id(0);
	ulong k_size = params->kmer_size;
    ulong offset = params->offset;
    ulong end = params->end;
    ulong t_work_items = params->t_work_items;
	ulong j, k;
	
	// Until reaching end of sequence
	//for(j=0; j<params->kmers_in_work_item; j++){
		
    // Coalescent
    ulong pos = offset + global_id;
    if(pos < end-k_size){


		ulong hash_full = 0;
		
		unsigned char bad = 0;

        if(pos < end-k_size){
        
            for(k=0; k<k_size; k++){
                
                switch(sequence[pos+k]){
                    case 'A': {}
                    break;
                    case 'C': hash_full += pow4[k]; 
                    break;
                    case 'G': hash_full += pow4[k] * 2; 
                    break;
                    case 'T': hash_full += pow4[k] * 3;
                    break;
                    case '\n': {}
                    break;
                    default: { bad = 1; }
                    break;
                }
            }
            
            
            if(bad == 0){
                words[pos-offset].h = hash_full;
                words[pos-offset].pos = pos;
            }else{
                words[pos-offset].h = 0xFFFFFFFFFFFFFFFF;
                words[pos-offset].pos = 0xFFFFFFFFFFFFFFFF;
            }
        }
	
    }else{
        words[pos-offset].h = 0xFFFFFFFFFFFFFFFF;
        words[pos-offset].pos = 0xFFFFFFFFFFFFFFFF;
    }

}