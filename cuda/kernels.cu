#include "kernels.cuh"

#define MIN_P_IDENT 80

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

//__device__ int printid = 0;

__device__ void left_extend(int32_t warp_pos_x_left, int32_t warp_pos_y_left, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t * best_offset_left, const char * seq_x, const char * seq_y)
{

	int32_t hp1 = warp_pos_x_left;
	*best_offset_left = (uint32_t) warp_pos_x_left;
	int32_t score = 32, best_score = 32, total_idents = 16; // Half of the hit for each side
	int p_ident = 100;
	int cell_score;
	int32_t thre_pos_x = warp_pos_x_left + threadIdx.x - (int32_t) x_seq_off;
	int32_t thre_pos_y = warp_pos_y_left + threadIdx.x - (int32_t) y_seq_off;

	while(score > 0 && (warp_pos_x_left - 32) >= (int32_t) x_seq_off && (warp_pos_y_left - 32) >= (int32_t) y_seq_off)
	{

		warp_pos_x_left -= 32;
		warp_pos_y_left -= 32;

		thre_pos_x -= 32;
		thre_pos_y -= 32;


		char v_x = seq_x[thre_pos_x];
		char v_y = seq_y[thre_pos_y];

		cell_score = (v_x == v_y);
		if(v_x == 'N' || v_y == 'N') cell_score = 0;

		for (int offset = 16; offset > 0; offset = offset >> 1)
			cell_score += __shfl_down_sync(0xFFFFFFFF, cell_score, offset);


		int idents = __shfl_sync(0xFFFFFFFF, cell_score, 0);
		score = score + (int32_t) idents;
		total_idents += (int32_t) idents;
		score = score - (int32_t) (32 - idents);
		p_ident = ((100 * total_idents) / (int) (hp1 - warp_pos_x_left + 16));
		if(score > best_score && p_ident >= MIN_P_IDENT){ best_score = score; *best_offset_left = warp_pos_x_left; }

	}

}

__device__ void right_extend(int32_t warp_pos_x_right, int32_t warp_pos_y_right, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t * best_offset_right, const char * seq_x, const char * seq_y, int32_t hp1, uint32_t x_lim, uint32_t y_lim)
{
	int32_t thre_pos_x = warp_pos_x_right + threadIdx.x - (int32_t) x_seq_off;
	int32_t thre_pos_y = warp_pos_y_right + threadIdx.x - (int32_t) y_seq_off;
	int32_t score = 32;
	int32_t total_idents = 16;
	int32_t best_score = 32;
	int p_ident = 100;

	while(score > 0 && (warp_pos_x_right + 32) < (int32_t) x_lim && (warp_pos_y_right + 32) < (int32_t) y_lim)
	{
		char v_x = seq_x[thre_pos_x];
		char v_y = seq_y[thre_pos_y];

		int cell_score = (v_x == v_y);
		//(v_x == 'N' || v_y == 'N') ? (cell_score = 0) : (0);
		if(v_x == 'N' || v_y == 'N') cell_score = 0;

		for (int offset = 16; offset > 0; offset = offset >> 1)
			cell_score += __shfl_down_sync(0xFFFFFFFF, cell_score, offset);

		int idents = __shfl_sync(0xFFFFFFFF, cell_score, 0);
		score = score + (int32_t) idents;
		total_idents += (int32_t) idents;
		score = score - (int32_t) (32 - idents);
		p_ident = (int) ((100 * total_idents) / (int32_t) (warp_pos_x_right - (hp1) + 16));

		warp_pos_x_right += 32;
		warp_pos_y_right += 32;
		thre_pos_x += 32;
		thre_pos_y += 32;

		if(score > best_score && p_ident >= MIN_P_IDENT){ best_score = score; *best_offset_right = warp_pos_x_right; }

	}
}

__global__ void kernel_frags_forward_register(uint32_t * h_p1, uint32_t * h_p2, uint32_t * left_offset, uint32_t * right_offset, const char * seq_x, const char * seq_y, uint32_t query_len, uint32_t ref_len, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t x_lim, uint32_t y_lim, uint32_t n_hits_kept, uint32_t n_frags_per_block){

	__shared__ int32_t hp1_shared[32];
	__shared__ int32_t hp2_shared[32];
	__shared__ int32_t diags_shared[32];
	__shared__ int32_t l_o[32];
	__shared__ int32_t r_o[32];


    int32_t curr_frag = 0;
    int32_t real_id = blockIdx.x * 32;
    int32_t prev_d = 0xFFFFFFFF, x_begin = 0xFFFFFFFF, x_finish = 0xFFFFFFFF;


	hp1_shared[threadIdx.x]   = (int32_t) h_p1[real_id + threadIdx.x];
	hp2_shared[threadIdx.x]   = (int32_t) h_p2[real_id + threadIdx.x];
	diags_shared[threadIdx.x] = hp1_shared[threadIdx.x] - hp2_shared[threadIdx.x];


    while(curr_frag < 32)
    {

		if(real_id >= n_hits_kept) break;

		int32_t hp1 = hp1_shared[curr_frag];
		uint32_t best_offset_left = (uint32_t) hp1;

        if(prev_d == diags_shared[curr_frag] && x_begin <= hp1 && hp1 <= x_finish) {
			l_o[curr_frag] = 0;
	        r_o[curr_frag] = 0;
			++real_id; ++curr_frag;  continue; 
		}

		prev_d = diags_shared[curr_frag];
		uint32_t best_offset_right = (uint32_t) hp1 + 32;


		left_extend(hp1, hp2_shared[curr_frag], x_seq_off, y_seq_off, &best_offset_left, seq_x, seq_y);

		right_extend(best_offset_right, hp2_shared[curr_frag] + 32, x_seq_off, y_seq_off, &best_offset_right, seq_x, seq_y, hp1, x_lim, y_lim);

		x_begin  = best_offset_left;
		x_finish = best_offset_right;

		l_o[curr_frag] = hp1 - (uint32_t) best_offset_left;
		r_o[curr_frag] = (uint32_t) (best_offset_right) - hp1;

		++curr_frag;
		++real_id;

	}

	real_id = blockIdx.x * 32;
	left_offset[real_id + threadIdx.x]  = l_o[threadIdx.x];
	right_offset[real_id + threadIdx.x] = r_o[threadIdx.x];


}

__global__ void kernel_frags_forward_per_thread(uint32_t * h_p1, uint32_t * h_p2, uint32_t * left_offset, uint32_t * right_offset, const char * seq_x, const char * seq_y, uint32_t query_len, uint32_t ref_len, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t x_lim, uint32_t y_lim, uint32_t n_hits){


		// LEFT ALIGNMENT

		int32_t id = threadIdx.x + blockIdx.x * blockDim.x;
		if(id >= n_hits) return;


		int32_t hp1 = (int32_t) h_p1[id];
		int32_t warp_pos_x_left = hp1;
		int32_t warp_pos_y_left = (int32_t) h_p2[id];
		uint32_t best_offset_left = (uint32_t) hp1;
		int32_t score = 32, best_score = 32, total_idents = 16; // Half of the hit for each side
		int p_ident = 100;
		int cell_score;

		
		while(score > 0 && (warp_pos_x_left-1) >= (int32_t) x_seq_off && (warp_pos_y_left-1) >= (int32_t) y_seq_off)
		{
			
			--warp_pos_x_left;
			--warp_pos_y_left;
		

			char v_x = seq_x[warp_pos_x_left - (int32_t) x_seq_off];
			char v_y = seq_y[warp_pos_y_left - (int32_t) y_seq_off];

			cell_score = (v_x == v_y);

			if(cell_score == 0) --score; else ++score;
		total_idents += (int32_t) cell_score;
		p_ident = ((100 * total_idents) / (int) (hp1 - warp_pos_x_left + 16));
		if(score > best_score && p_ident >= MIN_P_IDENT){ best_score = score; best_offset_left = warp_pos_x_left; }
		
		
	}
	
	// RIGHT ALIGNMENT
	int32_t warp_pos_x_right = hp1 + 32;
	int32_t warp_pos_y_right = (int32_t) h_p2[id] + 32;
	uint32_t best_offset_right = (uint32_t) warp_pos_x_right;
	score = 32;
	total_idents = 16;
	best_score = 32;
	p_ident = 100;

	
	while(score > 0 && (warp_pos_x_right + 1) < (int32_t) x_lim && (warp_pos_y_right + 1) < (int32_t) y_lim)
	{
		++warp_pos_x_right;
		++warp_pos_y_right;
		char v_x = seq_x[warp_pos_x_right - (int32_t) x_seq_off];
		char v_y = seq_y[warp_pos_y_right - (int32_t) y_seq_off];

		//cell_score = (v_x == v_y && v_x != '\0' && v_y != '\0');
		cell_score = (v_x == v_y);

		if(cell_score == 0) --score; else ++score;

		total_idents += (int32_t) cell_score;
		p_ident = (int) ((100 * total_idents) / (int32_t) (warp_pos_x_right - (hp1) + 16));

		if(score > best_score && p_ident >= MIN_P_IDENT){ best_score = score; best_offset_right = warp_pos_x_right; }

	}
	
	
    
	left_offset[id] = hp1 - (uint32_t) best_offset_left;
	//printf("BLOCK %d %d left offset:  %u\n", blockIdx.x, counter++, left_offset[blockIdx.x]);
	right_offset[id] = (uint32_t) (best_offset_right) - hp1;
	//printf("BLOCK %d %d right offset: %u\n", blockIdx.x, counter++, right_offset[blockIdx.x]);

}


__global__ void kernel_frags_reverse_register(uint32_t * h_p1, uint32_t * h_p2, uint32_t * left_offset, uint32_t * right_offset, const char * seq_x, const char * seq_y, uint32_t query_len, uint32_t ref_len, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t x_lim, uint32_t y_lim, uint32_t n_hits_kept, uint32_t n_frags_per_block){


	int32_t per_block = (int32_t) n_frags_per_block;
	int32_t curr_frag = 0;
    int32_t real_id = blockIdx.x * per_block;

    int32_t prev_d = 0xFFFFFFFF, x_begin = 0xFFFFFFFF, x_finish = 0xFFFFFFFF;

    while(curr_frag < per_block)
    {

        if(real_id >= n_hits_kept) return;


		// LEFT ALIGNMENT
		int32_t hp1 = (int32_t) h_p1[real_id];
		int32_t warp_pos_x_left = hp1;
		int32_t warp_pos_y_left = (int32_t) h_p2[real_id];
		//int32_t thre_pos_x, thre_pos_y;
		int32_t thre_pos_x = warp_pos_x_left + threadIdx.x - (int32_t) x_seq_off;
		int32_t thre_pos_y = warp_pos_y_left + threadIdx.x - (int32_t) y_seq_off;

		uint32_t best_offset_left = (uint32_t) hp1;
		int32_t score = 32, best_score = 32, total_idents = 16;
		int p_ident = 100;
		int cell_score;

		// Warp_pos_y_left is hp2[blockIdx.x]
        if(prev_d == hp1-warp_pos_y_left && x_begin <= hp1 && hp1 <= x_finish) { ++real_id; ++curr_frag;  continue; }

        // Save info on this frag
        prev_d = hp1 - warp_pos_y_left;

		
		while(score > 0 && (warp_pos_x_left - 32) >= (int32_t) x_seq_off && (warp_pos_y_left - 32) >= (int32_t) y_seq_off)
		{
			
			warp_pos_x_left -= 32;
			warp_pos_y_left -= 32;
			//thre_pos_x = warp_pos_x_left + threadIdx.x;
			//thre_pos_y = warp_pos_y_left + threadIdx.x;

			thre_pos_x -= 32;
			thre_pos_y -= 32;

			char v_x = seq_x[thre_pos_x];
			char v_y = seq_y[thre_pos_y];	
			//char v_x = seq_x[thre_pos_x - (int32_t) x_seq_off];
			//char v_y = seq_y[thre_pos_y - (int32_t) y_seq_off];

			
			//if(blockIdx.x == 0) printf("T%d -> %c - %c\n", threadIdx.x, v_x, v_y);
			

			//cell_score = (v_x == v_y && v_x != '\0' && v_y != '\0');
			cell_score = (v_x == v_y);
			(v_x == 'N' || v_y == 'N') ? (cell_score = 0) : (0);

			for (int offset = 16; offset > 0; offset = offset >> 1)
				cell_score += __shfl_down_sync(0xFFFFFFFF, cell_score, offset);

			
			int idents = __shfl_sync(0xFFFFFFFF, cell_score, 0);

			score = score + (int32_t) idents;
			total_idents += (int32_t) idents;
			score = score - (int32_t) (32 - idents);
			p_ident = ((100 * total_idents) / (int) (hp1 - warp_pos_x_left + 16));
			if(score > best_score && p_ident >= MIN_P_IDENT){ best_score = score; best_offset_left = warp_pos_x_left; }
		
			
			

		}
		
		
		
		//if(threadIdx.x==0) printf("BLOCK %d %d end left al\n", blockIdx.x, counter++);

		// RIGHT ALIGNMENT
		int32_t warp_pos_x_right = hp1 + 32;
		int32_t warp_pos_y_right = (int32_t) h_p2[real_id] + 32;
		uint32_t best_offset_right = (uint32_t) warp_pos_x_right;

		thre_pos_x = warp_pos_x_right + threadIdx.x - (int32_t) x_seq_off;
		thre_pos_y = warp_pos_y_right + threadIdx.x - (int32_t) y_seq_off;

		score = 32;
		best_score = 32;
		p_ident = 100;
		total_idents = 16;

		
		//while(p_ident > MIN_P_IDENT && (warp_pos_x_right + 32) < (int64_t) x_lim && (warp_pos_y_right + 32) < (int64_t) y_lim)
		while(score > 0 && (warp_pos_x_right + 32) < (int32_t) x_lim && (warp_pos_y_right + 32) < (int32_t) y_lim)
		{
			//thre_pos_x = warp_pos_x_right + threadIdx.x;
			//thre_pos_y = warp_pos_y_right + threadIdx.x;
			//char v_x = seq_x[thre_pos_x - (int32_t) x_seq_off];
			//char v_y = seq_y[thre_pos_y - (int32_t) y_seq_off];

			char v_x = seq_x[thre_pos_x];
			char v_y = seq_y[thre_pos_y];

			//cell_score = (v_x == v_y && v_x != '\0' && v_y != '\0');
			cell_score = (v_x == v_y);
			(v_x == 'N' || v_y == 'N') ? (cell_score = 0) : (0);

			for (int offset = 16; offset > 0; offset = offset >> 1)
				cell_score += __shfl_down_sync(0xFFFFFFFF, cell_score, offset);
			
			int idents = __shfl_sync(0xFFFFFFFF, cell_score, 0);
			score = score + (int32_t) idents;
			total_idents += (int32_t) idents;
			score = score - (int32_t) (32 - idents);
			p_ident = (int) ((100 * total_idents) / (int32_t) (warp_pos_x_right - (hp1) + 16));

			warp_pos_x_right += 32;
			warp_pos_y_right += 32; 
			thre_pos_x += 32;
			thre_pos_y += 32;

			if(score > best_score && p_ident >= MIN_P_IDENT){ best_score = score; best_offset_right = warp_pos_x_right; }

		}
		
		//if(threadIdx.x==0) printf("BLOCK %d %d end right al\n", blockIdx.x, counter++);

        x_begin  = best_offset_left;
        x_finish = best_offset_right;
        ++curr_frag;
		

		// Save at the end
		if(threadIdx.x == 0){
			left_offset[real_id] = hp1 - (uint32_t) best_offset_left;
			//right_offset[blockIdx.x] = (uint64_t) (best_offset_right + 32) - h_p1[blockIdx.x];
			right_offset[real_id] = (uint32_t) (best_offset_right) - hp1;
		}

		++real_id;

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
	// The last INIT kmer starts on 124 (that is 128-4), located on thread 15 on the fourth byte (124/8 = 15.5, 0.5 * 8 = 4)
	// So we distribute the next bytes to the other threads by shuffle
	// Thread 0 will get bytes from thread (124+0)/8 = 15.5 -> thread 15, 
	// and so will the following 4 threads until thread 16 is reached and so on
	temp_value = __shfl_sync(0xFFFFFFFF, value, (124 + threadIdx.x) >> 3);
	// We fetched the full 8 bytes from the value variable, now we need to select only bytes 124
	// Do so by fetching byte i.e. thread 0 has fetched bytes from ULLI from thread 15
	// That is 15*8 = 120 to 127 bytes, but it needs to stick to 124
	//int remainder = ((124 + threadIdx.x) & 2) << 3; // This is % 3 and multiply by 8
	// Which results in 4 for thread 0

	byte = (char) (temp_value >> ((threadIdx.x & 7) << 3));
	//byte = (char) (temp_value >> (8*(remainder-1)));
	

	//if(blockIdx.x == 0) printf("%d -> %c (my val: %d)\n", threadIdx.x, byte, (124 + threadIdx.x) >> 3);

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
	
	
	//printf("PREV at %" PRIu64" -> %.32s\n", 0, &sequence[0]);
	
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

	//if(threadIdx.x == 0 && blockIdx.x == 0) printf("POST at %" PRIu64" -> %.32s @ %" PRIu64"\n", threadIdx.x + (32) + 128 * blockIdx.x, &sequence[threadIdx.x + (32) + 128 * blockIdx.x], hashes[threadIdx.x + (32) + 128 * blockIdx.x]);
	
}


__global__ void kernel_index_global32(uint64_t * hashes, uint32_t * positions, const char * sequence, uint32_t offset) {

    uint64_t hash = 0, k = 0;
    uint64_t bad = 0xFFFFFFFFFFFFFFFF;
    for(k=0; k<32; k++){
        char c = sequence[threadIdx.x + k + blockIdx.x * blockDim.x];
        // IF-binary
        
		/*	
        hash = hash << 2;
        if(c == 'A') hash += 0;
        if(c == 'C') hash += 1;
        if(c == 'G') hash += 2;
        if(c == 'T') hash += 3;
        if(c == 'N') bad = 0;
		*/
       
	
		
        // IF-cached
		
        if(c == 'C') hash += pow4[k];
        if(c == 'G') hash += pow4_G[k];
        if(c == 'T') hash += pow4_T[k];
        if(c == 'N') bad = 0;
		

       	 
        // Pure arithmetic
		/*
        hash = hash << 2;
        hash += (c & 6) >> 1;
        if(c == 'N') bad = 0;
		*/
        

    }
    // [0 - 32] * [0-N]* [32]
    hashes[threadIdx.x + blockIdx.x * blockDim.x] = hash & bad;
    positions[threadIdx.x + blockIdx.x * blockDim.x] = (threadIdx.x + blockIdx.x * blockDim.x + offset) | (~bad);
}

__global__ void kernel_index_global32_advanced(uint64_t * hashes, uint32_t * positions, const uchar4 * sequence, uint32_t offset) {
		

	uint64_t hash = 0;
    int k;

    int access = (threadIdx.x >> 2) + (blockIdx.x >> 2) * (blockDim.x >> 2);
		
	uint64_t bad = 0xFFFFFFFFFFFFFFFF;

	for(k=0; k<8; k++){

		//char c = sequence[threadIdx.x + k + blockIdx.x * blockDim.x];
        int superk = k << 2;
        uchar4 c = sequence[access + k];

        //if(blockIdx.x == 0 && threadIdx.x == 1) printf("%c %c %c %c\n", c.x, c.y, c.z, c.w);

        if(c.x == 'A') hash += 0;
        if(c.x == 'C') hash += pow4[superk];
        if(c.x == 'G') hash += pow4_G[superk];
        if(c.x == 'T') hash += pow4_T[superk];
        if(c.x == 'N') bad = 0;
        if(c.y == 'A') hash += 0;
        if(c.y == 'C') hash += pow4[superk + 1];
        if(c.y == 'G') hash += pow4_G[superk + 1];
        if(c.y == 'T') hash += pow4_T[superk + 1];
        if(c.y == 'N') bad = 0;
        if(c.z == 'A') hash += 0;
        if(c.z == 'C') hash += pow4[superk + 2];
        if(c.z == 'G') hash += pow4_G[superk + 2];
        if(c.z == 'T') hash += pow4_T[superk + 2];
        if(c.z == 'N') bad = 0;
        if(c.w == 'A') hash += 0;
        if(c.w == 'C') hash += pow4[superk + 3];
        if(c.w == 'G') hash += pow4_G[superk + 3];
        if(c.w == 'T') hash += pow4_T[superk + 3];
        if(c.w == 'N') bad = 0;

	}

	// [0 - 32] * [0-N]* [32]
	hashes[threadIdx.x + blockIdx.x * blockDim.x] = hash & bad;
	positions[threadIdx.x + blockIdx.x * blockDim.x] = (threadIdx.x + blockIdx.x * blockDim.x + offset) | (~bad);
}

__global__ void kernel_reverse_complement(const char * sequence, char * reverse_sequence, uint32_t seq_len) {
	
	int32_t id = threadIdx.x + blockIdx.x * blockDim.x;
	int32_t sslen = (int32_t) seq_len;

	if(id < seq_len){
    // This if is not bad since only the last block will not work it in lockstep

		//uint64_t lookup = seq_len - id;
		int32_t lookup = (sslen - 1) - id;
		char original = sequence[id];

        
		char complement = 'N';
		if(original == 'A') complement = 'T';
		if(original == 'C') complement = 'G';
		if(original == 'G') complement = 'C';
		if(original == 'T') complement = 'A';
        
        
		reverse_sequence[lookup] = complement;
	}
}

