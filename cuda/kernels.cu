#include "kernels.cuh"

__constant__ uint64_t pow4[33]={1L, 4L, 16L, 64L, 256L, 1024L, 4096L, 16384L, 65536L,
    262144L, 1048576L, 4194304L, 16777216L, 67108864L, 268435456L, 1073741824L, 4294967296L,
    17179869184L, 68719476736L, 274877906944L, 1099511627776L, 4398046511104L, 17592186044416L,
    70368744177664L, 281474976710656L, 1125899906842624L, 4503599627370496L, 18014398509481984L,
	72057594037927936L, 288230376151711744L, 1152921504606846976L, 4611686018427387904L};
	

__constant__ uint64_t pow4_G[33]={2*1L, 2*4L, 2*16L, 2*64L, 2*256L, 2*1024L, 2*4096L, 2*16384L, 2*65536L,
    (uint64_t)2*262144L, (uint64_t)2*1048576L,(uint64_t)2*4194304L, (uint64_t)2*16777216L, (uint64_t)2*67108864L, (uint64_t)2*268435456L, (uint64_t)2*1073741824L, (uint64_t)2*4294967296L,
    (uint64_t)2*17179869184L, (uint64_t)2*68719476736L, (uint64_t)2*274877906944L, (uint64_t)2*1099511627776L, (uint64_t)2*4398046511104L, (uint64_t)2*17592186044416L,
    (uint64_t)2*70368744177664L, (uint64_t)2*281474976710656L, (uint64_t)2*1125899906842624L, (uint64_t)2*4503599627370496L, (uint64_t)2*18014398509481984L,
    (uint64_t)2*72057594037927936L, (uint64_t) 2*288230376151711744L, (uint64_t) 2*1152921504606846976L, (uint64_t) 2*4611686018427387904L};

__constant__ uint64_t pow4_T[33]={3*1L, 3*4L, 3*16L, 3*64L, 3*256L, 3*1024L, 3*4096L, 3*16384L, 3*65536L,
    (uint64_t)3*262144L, (uint64_t) 3*1048576L, (uint64_t)3*4194304L, (uint64_t)3*16777216L, (uint64_t)3*67108864L, (uint64_t)3*268435456L, (uint64_t)3*1073741824L, (uint64_t)3*4294967296L,
    (uint64_t)3*17179869184L, (uint64_t)3*68719476736L, (uint64_t)3*274877906944L, (uint64_t)3*1099511627776L, (uint64_t)3*4398046511104L, (uint64_t)3*17592186044416L,
    (uint64_t)3*70368744177664L, (uint64_t)3*281474976710656L, (uint64_t)3*1125899906842624L, (uint64_t)3*4503599627370496L, (uint64_t)3*18014398509481984L,
    (uint64_t)3*72057594037927936L, (uint64_t) 3*288230376151711744L, (uint64_t) 3*1152921504606846976L, (uint64_t) 3*4611686018427387904L};

#define MIN_P_IDENT 75

__global__ void kernel_frags_forward_register(uint64_t * h_p1, uint64_t * h_p2, uint64_t * left_offset, uint64_t * right_offset, const char * seq_x, const char * seq_y, uint64_t query_len, uint64_t ref_len){


	// LEFT ALIGNMENT
	int64_t warp_pos_x_left = h_p1[blockIdx.x];
	int64_t warp_pos_y_left = h_p2[blockIdx.x];
	int64_t thre_pos_x, thre_pos_y;
	int score = 0;
	int p_ident = 100;
	int cell_score;

	//if(threadIdx.x == 0) printf("checking hit %"PRIu64", %"PRIu64"\n", h_p1[blockIdx.x], h_p2[blockIdx.x]);

	while(p_ident > MIN_P_IDENT && (warp_pos_x_left - 32) >= 0 && (warp_pos_y_left - 32) >= 0)
	{
		warp_pos_x_left -= 32;
		warp_pos_y_left -= 32;
		thre_pos_x = warp_pos_x_left + threadIdx.x;
		thre_pos_y = warp_pos_y_left + threadIdx.x;
		char v_x = seq_x[thre_pos_x];
		char v_y = seq_y[thre_pos_y];

		cell_score = (v_x == v_y);

		//printf("I am %d :: %c and %c -> %d\n", threadIdx.x, v_x, v_y, cell_score);

		for (int offset = 16; offset > 0; offset = offset >> 1)
			cell_score += __shfl_down_sync(0xFFFFFFFF, cell_score, offset);

		//if(threadIdx.x == 0) printf("Sum is %d\n", cell_score);
		int idents = __shfl_sync(0xFFFFFFFF, cell_score, 0);
		score = score + idents;
		p_ident = (100 * score) / (int) (h_p1[blockIdx.x] - thre_pos_x);
		
		//if(threadIdx.x == 0 && blockIdx.x < 20) printf("Pident %d\n", p_ident);



		//if(threadIdx.x == 3) printf("Final score %d\n", score);

	}

	// RIGHT ALIGNMENT
	int64_t warp_pos_x_right = h_p1[blockIdx.x] + 32;
	int64_t warp_pos_y_right = h_p2[blockIdx.x] + 32;
	score = 0;
	p_ident = 100;

	while(p_ident > MIN_P_IDENT && (warp_pos_x_right + 32) < query_len && (warp_pos_y_right + 32) < ref_len)
	{
		thre_pos_x = warp_pos_x_right + threadIdx.x;
		thre_pos_y = warp_pos_y_right + threadIdx.x;
		char v_x = seq_x[thre_pos_x];
		char v_y = seq_y[thre_pos_y];

		cell_score = (v_x == v_y);

		for (int offset = 16; offset > 0; offset = offset >> 1)
			cell_score += __shfl_down_sync(0xFFFFFFFF, cell_score, offset);
		
		int idents = __shfl_sync(0xFFFFFFFF, cell_score, 0);
		score = score + idents;
		p_ident = (100 * score) / (int) (thre_pos_x - (h_p1[blockIdx.x] + 32));

		warp_pos_x_right += 32;
		warp_pos_y_right += 32;

	}

	// Save at the end
	if(threadIdx.x == 0){
		left_offset[blockIdx.x] = h_p1[blockIdx.x] - (uint64_t) warp_pos_x_left;
		right_offset[blockIdx.x] = (uint64_t) (warp_pos_x_right + 32) - h_p1[blockIdx.x];
	}

}


__global__ void kernel_register_fast_hash_rotational(uint64_t * hashes, uint64_t * positions, const char * sequence, uint64_t offset) {
	


	int i;
	ULLI hash = 0;

	// Each will do 4 kmers
	// So 32 threads do 32*4 = 128 kmers
	// So last kmer processed starts at position 127+32 = we need 160 bytes
	// 160 / 8 = 20
	

	ULLI value = ((ULLI *)sequence)[threadIdx.x + blockIdx.x * 16 ]; // 16*8 = 128 are the bytes 
	ULLI temp_value = value;
	ULLI segment[2] = {0, 0};
	char byte;
	ULLI bad;
	

	i = 0;

	int index = threadIdx.x >> 1; // This is its ID / 2
	int mod2 = (threadIdx.x & 1); // This is threadIdx.x % 2
	int mod2_times_32 = (mod2 << 5); // This is (threadIdx.x % 2) * 32, i.e. a 0 if threadIdx is even and a 32 if its odd

	hash = 0;
	bad = 0xFFFFFFFFFFFFFFFF;

	// Compute the last kmer so that we can make a collaborative algorithm
	// The last INIT kmer starts on 124, located on thread 15 on the fourth byte
	temp_value = __shfl_sync(0xFFFFFFFF, value, (124 + threadIdx.x) >> 3);
	byte = (char) (temp_value >> ((threadIdx.x & 7) << 3));

	if(byte == 'C') hash += pow4[threadIdx.x];
	if(byte == 'G') hash += pow4_G[threadIdx.x];
	if(byte == 'T') hash += pow4_T[threadIdx.x];

	// Sum segments and store them in second array 
	//sum (ID+1)*4 mod 32 … + 0,1,2,3

	int offset_segment = (((threadIdx.x + 1) << 2) & 31);

	segment[1] += __shfl_sync(0xFFFFFFFF, hash, offset_segment);
	segment[1] = segment[1] << 2;
	segment[1] += __shfl_sync(0xFFFFFFFF, hash, offset_segment + 1);
	segment[1] = segment[1] << 2;
	segment[1] += __shfl_sync(0xFFFFFFFF, hash, offset_segment + 2);
	segment[1] = segment[1] << 2;
	segment[1] += __shfl_sync(0xFFFFFFFF, hash, offset_segment + 3);
	segment[1] = segment[1] << 2;

	hash = 0;
	
	
	

	// The 3 is because it is 0011 in binary so it only selects the destination and source thread with the mask
	temp_value = __shfl_sync(0xFFFFFFFF, value, index);

	
	//byte = ((char*) temp_value)[0 + mod2_times_four]; hash = hash << 2;
	byte = (char) (temp_value >> (24 + mod2_times_32)); hash = hash << 2;
	if((char) byte == 'C') hash = hash + 1; if((char) byte == 'G') hash = hash + 2;
	if((char) byte == 'T') hash = hash + 3; if((char) byte == 'N') bad = 0;

	byte = (char) (temp_value >> (16 + mod2_times_32)); hash = hash << 2;
	if((char) byte == 'C') hash = hash + 1; if((char) byte == 'G') hash = hash + 2;
	if((char) byte == 'T') hash = hash + 3; if((char) byte == 'N') bad = 0;

	byte = (char) (temp_value >> (8 + mod2_times_32)); hash = hash << 2;
	if((char) byte == 'C') hash = hash + 1; if((char) byte == 'G') hash = hash + 2;
	if((char) byte == 'T') hash = hash + 3; if((char) byte == 'N') bad = 0;

	byte = (char) (temp_value >> (0 + mod2_times_32)); hash = hash << 2;
	if((char) byte == 'C') hash = hash + 1; if((char) byte == 'G') hash = hash + 2;
	if((char) byte == 'T') hash = hash + 3; if((char) byte == 'N') bad = 0;

	segment[0] = hash;

	// To this point we have calculated a 4-byte segment in each thread (along with the complete hash of initial-last kmer (no. 124))
	// Now we can begin the "rotations"
	// Each thread will receive the hash value of the next thread

	
	//offset_segment = 0 + ((threadIdx.x + 1) >> 5);
	
	// THis results in a local (global) load because offset is not known at compile time
	// I think the rotational algorithm is better suited for shared memory?
	// With the IFs it doesnt, because it knows the memory address (thread number does not matter)
	
	if(threadIdx.x < 31){

		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[0], (threadIdx.x + 1) & 31)) << (8)); // 64 
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[0], (threadIdx.x + 2) & 31)) << (16));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[0], (threadIdx.x + 3) & 31)) << (24));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[0], (threadIdx.x + 4) & 31)) << (32));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[0], (threadIdx.x + 5) & 31)) << (40));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[0], (threadIdx.x + 6) & 31)) << (48));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[0], (threadIdx.x + 7) & 31)) << (56));

	}else{
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[1], (threadIdx.x + 1) & 31)) << (8)); // 64 
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[1], (threadIdx.x + 2) & 31)) << (16));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[1], (threadIdx.x + 3) & 31)) << (24));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[1], (threadIdx.x + 4) & 31)) << (32));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[1], (threadIdx.x + 5) & 31)) << (40));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[1], (threadIdx.x + 6) & 31)) << (48));
		hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[1], (threadIdx.x + 7) & 31)) << (56));
	}
	
		
	/*
	hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[offset_segment], (threadIdx.x + 1) & 31)) << (8)); // 64 
	hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[offset_segment], (threadIdx.x + 2) & 31)) << (16));
	hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[offset_segment], (threadIdx.x + 3) & 31)) << (24));
	hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[offset_segment], (threadIdx.x + 4) & 31)) << (32));
	hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[offset_segment], (threadIdx.x + 5) & 31)) << (40));
	hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[offset_segment], (threadIdx.x + 6) & 31)) << (48));
	hash = hash + ((__shfl_sync(0xFFFFFFFF, segment[offset_segment], (threadIdx.x + 7) & 31)) << (56));
	*/
	
	
	//printf("PREV at %"PRIu64" -> %.32s\n", 0, &sequence[0]);
	
	//table[threadIdx.x + 32*i + 192 * blockIdx.x] = hash & bad;
	hashes[threadIdx.x + 128 * blockIdx.x] = hash & bad;
	positions[threadIdx.x + 128 * blockIdx.x] = (threadIdx.x + blockIdx.x * 128 + offset) | (~bad);
	//table[threadIdx.x + 96 * blockIdx.x] = hash & bad;
	

	unsigned kmer_start = (32 + (threadIdx.x << 2)); // 4 because of 4 kmers per thread
	unsigned int_pos = kmer_start >> 3; // 8 because bytes per register

	for(i=1; i<4; i++){
		
		bad = 0xFFFFFFFFFFFFFFFF;
		temp_value = __shfl_sync(0xFFFFFFFF, value, int_pos);
		byte = (char) (temp_value >> ((kmer_start & 7) << 3));
		++kmer_start;
		int_pos = kmer_start >> 3;

		// Remember - the kmer is not well built right now but first focus on performance

		
		hash = hash << 2;

		
		if((char) byte == 'C') hash = hash + 1;
		if((char) byte == 'G') hash = hash + 2;
		if((char) byte == 'T') hash = hash + 3;
		if((char) byte == 'N') bad = 0;
		

		/*
		hash += (byte & 6) >> 1;
		if(byte == 'N') bad = 0;
		*/

		//table[threadIdx.x + blockIdx.x * blockDim.x*8 + blockDim.x * k] = hash & bad;
		hashes[threadIdx.x + (i << 5) + 128 * blockIdx.x] = hash & bad;
		positions[threadIdx.x + (i << 5) + 128 * blockIdx.x] = (threadIdx.x + (i << 5) + blockIdx.x * 128 + offset) | (~bad);

		
	}

	if(threadIdx.x == 0 && blockIdx.x == 0) printf("POST at %"PRIu64" -> %.32s @ %"PRIu64"\n", threadIdx.x + (32) + 128 * blockIdx.x, &sequence[threadIdx.x + (32) + 128 * blockIdx.x], hashes[threadIdx.x + (32) + 128 * blockIdx.x]);
	
}

__global__ void kernel_index_global32(uint64_t * hashes, uint64_t * positions, const char * sequence, uint64_t offset) {
		

	uint64_t k, hash = 0;
		
	uint64_t bad = 0xFFFFFFFFFFFFFFFF;

	for(k=0; k<32; k++){

		char c = sequence[threadIdx.x + k + blockIdx.x * blockDim.x];
		//if(threadIdx.x == 0) printf("%c", c);
		//char c = sequence[0];

		if(c == 'A') hash += 0;
		if(c == 'C') hash += pow4[k];
		if(c == 'G') hash += pow4_G[k];
		if(c == 'T') hash += pow4_T[k];
		if(c == 'N') bad = 0;
		
	}

	// [0 - 32] * [0-N]* [32]
	hashes[threadIdx.x + blockIdx.x * blockDim.x] = hash & bad;
	//hashes[0] = hash & bad;
	positions[threadIdx.x + blockIdx.x * blockDim.x] = (threadIdx.x + blockIdx.x * blockDim.x + offset) | (~bad);
}
