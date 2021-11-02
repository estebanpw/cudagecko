#include "kernels.cuh"

#define MIN_P_IDENT 80

// LOOK UP TABLE

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

__global__ void kernel_filter_hits_parallel(uint64_t * diagonals_merged, uint32_t ref_len, uint32_t n_hits_found)
{
    uint64_t mask_d  = 0xFFFFFFFF00000000; 
    uint64_t mask_px = 0x00000000FFFFFFFF;

    uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index >= n_hits_found) return;

    uint32_t d_prev, px_prev;

    uint64_t v = diagonals_merged[index];
    uint32_t d_curr = (uint32_t) ((v & mask_d) >> 32);
    uint32_t px = (uint32_t) (v & mask_px);
    uint32_t py = ref_len + px - d_curr;
    uint64_t res = ((uint64_t) px << (uint64_t) 32) + (uint64_t) py;

    int32_t start = (threadIdx.x < 32) ? 0 : 32;
    //for(int32_t i=0; i<blockDim.x-1; i++){
    for(int32_t i=start; i<start+32; i++){

        d_prev  = __shfl_sync(0xFFFFFFFF, d_curr, i);
        px_prev = __shfl_sync(0xFFFFFFFF, px, i); 

        if(threadIdx.x > i && d_prev == d_curr && px < (px_prev+63)){
            res = 0xFFFFFFFFFFFFFFFF; // Last diagonal and last hash value means filtered
            px = 0xFFFFFFFF; // reset the position and diagonal so that it doesnt break other hits
            d_curr = 0xFFFFFFFF;
        }
    }
    diagonals_merged[index] = res;
}

__global__ void kernel_filter_hits(uint64_t * diagonals_merged, uint32_t ref_len, uint32_t n_hits_found)
{
    uint64_t mask_d  = 0xFFFFFFFF00000000; 
    uint64_t mask_px = 0x00000000FFFFFFFF;

    uint32_t index = blockIdx.x * 1024 + threadIdx.x * blockDim.x;
    if(index >= n_hits_found) return;

    uint64_t v = diagonals_merged[index];
    uint32_t d_curr = (uint32_t) ((v & mask_d) >> 32);
    uint32_t px = (uint32_t) (v & mask_px);
    uint32_t py = ref_len + px - d_curr; //values_y[index]; 
    uint64_t res = ((uint64_t) px << (uint64_t) 32) + (uint64_t) py;

    diagonals_merged[index] = res;

    for(int i=1; i<blockDim.x; i++){

        if(index + i >= n_hits_found) return;

        v = diagonals_merged[index + i];
        uint32_t d_next = (uint32_t) ((v & mask_d) >> 32);
        uint32_t px_next = (uint32_t) (v & mask_px);
        uint32_t py_next = ref_len + px_next - d_next;//values_y[index + i]; 

        if(d_curr != d_next || px_next > (px+63)){
            res = ((uint64_t) px_next << (uint64_t) 32) + (uint64_t) py_next;
            diagonals_merged[index + i] = res;
            px = px_next;
            py = py_next;
            d_curr = d_next;
        }else{
            diagonals_merged[index + i] = 0xFFFFFFFFFFFFFFFF; // Last diagonal and last hash value means filtered
        }
    }
}

__device__ uint32_t binary_search_hits(uint64_t key_x, uint64_t * hashes_y, uint32_t limit_y)
{
    uint32_t l = 0, r = limit_y;
    while(l < r){
        uint32_t m = (l + r) >> 1;
        if(hashes_y[m] < key_x)
            l = m + 1;
        else
            r = m;
    }
    return l;
}


__global__ void kernel_find_leftmost_items(uint64_t * hashes_x, uint32_t * pos_x, uint64_t * hashes_y, uint32_t * pos_y, uint32_t limit_x, uint32_t limit_y)
{
    *pos_x = binary_search_hits(0xFFFFFFFFFFFFFFFF, hashes_x, limit_x);
    *pos_y = binary_search_hits(0xFFFFFFFFFFFFFFFF, hashes_y, limit_y);
}

/*

 * kernel_compact_hits -> Removes gaps between diagonal hits

 * Constraints:
    * 

 * @param hashes_x The hash keys for the x or query sequence in ascending order
 * @return Nothing
 */


__global__ void kernel_compact_hits(uint64_t * ptr_device_diagonals, uint32_t * ptr_hits_log, uint32_t * ptr_accum_log, uint32_t mem_block, uint64_t * ptr_copy_place_diagonals, uint32_t offset_remover)
{
    /*
    if(blockIdx.x == 8191 && threadIdx.x == 0){
        printf("Hi, im from block %d. I will copy %u hits from position %u to position %u\n", blockIdx.x, ptr_hits_log[blockIdx.x], blockIdx.x * mem_block, ptr_accum_log[blockIdx.x] - offset_remover);
    }
    */

    for(uint32_t i=0; i<ptr_hits_log[blockIdx.x]; i+=blockDim.x){
        
        if(i+threadIdx.x < ptr_hits_log[blockIdx.x]){
            ptr_copy_place_diagonals[ptr_accum_log[blockIdx.x] - offset_remover + i + threadIdx.x] = ptr_device_diagonals[blockIdx.x * mem_block + i + threadIdx.x];
        }
    }

}


/*

 * kernel_hits -> Computes the hits between two lists of sorted words

 * Constraints:
    * All "bad" words must have value 0xFF...FF and position value 0xFF..FF and therefore their sorted positions must be at the very end of the list

 * @param hashes_x The hash keys for the x or query sequence in ascending order
 * @param hashes_y The hash keys for the y or reference sequence in ascending order
 * @param positions_x The positional value paired to each query key
 * @param positions_y The positional value paired to each reference key
 * @param hit_write_section The memory block to store all matched hits in the form of a diagonal d=(rlen + x - y)<<32 + x
 * @param mem_block The size of memory section from @hit_write_section reserved for each block (in number of hits, not bytes)
 * @param limit_* The number of properly formed words for the query and reference as retrieved by the binary search (i.e. index where 0xFF..FF begins)
 * @param error Variable to store function-specific error and return it to the host control
 * @param ref_len The length of the reference sequence (required to calculate the diagonal value)
 * @param hits_log The resulting number of hits created by each block
 * @param atomic_distributer An integer variable in global memory that should be accessed atomically and used to distribute the extra large memory blocks of hits
 * @param auxiliary_hit_memory A special section of memory used in conjunction with an atomic operation to allow some blocks to write a larger amount of hits
 * @param extra_large_memory_block The size (in number of hits) that a block can have in the additional hit space
 * @param max_extra_sections The maximum number of memory sections that can be given additionally with the auxiliary_hit_memory
 * @param hits_log_extra Same as hits_log but for the extra space
 * @return Nothing, its a cuda kernel

*/

__global__ void kernel_hits(uint64_t * hashes_x, uint64_t * hashes_y, uint32_t * positions_x, uint32_t * positions_y, 
    uint64_t * hit_write_section, int32_t mem_block, uint32_t limit_x, uint32_t limit_y, int32_t * error, uint32_t ref_len,
    uint32_t * hits_log, int32_t * atomic_distributer, uint64_t * auxiliary_hit_memory, uint32_t extra_large_memory_block, 
    uint32_t max_extra_sections, uint32_t * hits_log_extra)//, uint64_t * messages)
{
    // !!!!!!!!!!!!!
    // TODO: make the x words also fully aligned to 256 
    // !!!!!!!!!!!!!

    uint32_t current_mem_block = (uint32_t) mem_block;
    uint32_t accumulated = 0;

    uint32_t batch = (uint32_t) limit_x / gridDim.x;
    uint32_t word_x_id = blockIdx.x * batch;
    
    __shared__ uint64_t cached_keys_x[32]; // Needs to change if blockDim > 32
    __shared__ uint64_t cached_keys_y[32]; 
    __shared__ uint32_t cached_pos_x[32];
    __shared__ uint32_t cached_pos_y[32];

    // Keeps track of where we are writing in the additional memory
    int32_t atomic_position = -1;

    uint64_t * ptr_to_write_section = &hit_write_section[blockIdx.x * current_mem_block]; // To be changed if needed in case of more hits

    uint32_t word_y_id = binary_search_hits(hashes_x[word_x_id], hashes_y, limit_y);
    //if(threadIdx.x == 0) printf("[%d] I found my hash [%" PRIu64"] [I was looking for %" PRIu64" at p=%u] at p=%u, and i work from [%u to %u]\n", blockIdx.x, hashes_y[word_y_id], hashes_x[word_x_id], word_x_id, word_y_id, word_x_id, word_x_id+batch);
    
    uint32_t misaligned = word_y_id % 256;
    if(misaligned < word_y_id) word_y_id = word_y_id - misaligned; // Make it aligned

    uint32_t i = 0;

    uint64_t current_save;

    while(i < batch)
    {
        cached_keys_x[threadIdx.x] = hashes_x[word_x_id+threadIdx.x];
        cached_pos_x[threadIdx.x]  = positions_x[word_x_id+threadIdx.x];

        for(uint32_t current_word_x=0; current_word_x<blockDim.x; current_word_x++)
        {
            uint64_t partial_diag = (((uint64_t)ref_len + (uint64_t)cached_pos_x[current_word_x]) << 32) + (uint64_t)cached_pos_x[current_word_x];

            uint32_t current_word_y_id = word_y_id;
            while (current_word_y_id < limit_y)
            {
                // Cache current block of y-words
                cached_keys_y[threadIdx.x] = hashes_y[current_word_y_id+threadIdx.x];
                cached_pos_y[threadIdx.x]  = positions_y[current_word_y_id+threadIdx.x];

                uint32_t hit = 0;
                // Match keys
                
                //if(cached_pos_x[current_word_x] != 0xFFFFFFFF && cached_pos_y[threadIdx.x] != 0xFFFFFFFF) // This IF is no longer necessary because of leftmost searches
                if(cached_keys_x[current_word_x] == cached_keys_y[threadIdx.x])
                {
                    hit = 1;
                    //ptr_to_write_section[accumulated + threadIdx.x] = partial_diag - (((uint64_t) cached_pos_y[threadIdx.x]) << 32 );
                    current_save = partial_diag - (((uint64_t) cached_pos_y[threadIdx.x]) << 32 );
                }
                uint32_t mingap = threadIdx.x;
                if(hit == 0) mingap = 0xFFFFFFFF;
                // Distribute where first match happens
                for (int offset = 16; offset > 0; offset = offset >> 1)
                    mingap = min(mingap, __shfl_down_sync(0xFFFFFFFF, mingap, offset));
                mingap = __shfl_sync(0xFFFFFFFF, mingap, 0);

                // And write hits
                if(hit == 1)
                    ptr_to_write_section[accumulated + threadIdx.x - mingap] = current_save;

                // Distribute number of matches
                for (int offset = 16; offset > 0; offset = offset >> 1)
                    hit += __shfl_down_sync(0xFFFFFFFF, hit, offset);
                hit = __shfl_sync(0xFFFFFFFF, hit, 0);

                accumulated += hit;

                if(accumulated >= (uint32_t) current_mem_block - blockDim.x) {
                    // Assign one of the large blocks
                    if(threadIdx.x == 0 && atomic_position == -1) hits_log[blockIdx.x] = (uint32_t) accumulated;
                    if(threadIdx.x == 0 && atomic_position > -1) hits_log_extra[atomic_position] = (uint32_t) accumulated;
                    if(threadIdx.x == 0) atomic_position = atomicAdd(atomic_distributer, 1);
                    atomic_position = __shfl_sync(0xFFFFFFFF, atomic_position, 0);

                    accumulated = 0;
                    ptr_to_write_section = &auxiliary_hit_memory[(uint32_t) atomic_position * extra_large_memory_block];

                    current_mem_block = extra_large_memory_block;

                    if(atomic_position >= max_extra_sections){
                        *error = - ((int32_t) blockIdx.x + 1); 
                        return;
                    }
                }

                
                // Control flow
                if(cached_keys_x[current_word_x] < cached_keys_y[0])
                {
                    if(current_word_x < blockDim.x-1 && cached_keys_x[current_word_x+1] > cached_keys_x[current_word_x]
                        && cached_keys_x[current_word_x+1] >= cached_keys_y[0]){
                            word_y_id = current_word_y_id;
                        }
                    break;
                }
                if(cached_keys_x[current_word_x] > cached_keys_y[0]){ // this condition has increased performance
                    word_y_id = current_word_y_id;
                }
                current_word_y_id += blockDim.x;
                

            }

        }
        //if(threadIdx.x == 0 && messageID<messageMAX) messages[messageID++] = (((uint64_t)3 << 32) | word_x_id); //printf("%d Fetched next x block %d!\n", messageID++, word_x_id);
        word_x_id += blockDim.x;
        i += blockDim.x;
    }

    if(threadIdx.x == 0 && atomic_position > -1) hits_log_extra[atomic_position] = (uint32_t) accumulated;
    if(threadIdx.x == 0 && atomic_position == -1) hits_log[blockIdx.x] = (uint32_t) accumulated;



    /*
                // Very fast on some, on others it goes very slow [T3]
                // I know it is somehow related to length, but not to number of hits (on bosta it goes well, despite having more hits than galga)
                // so it must be linked to the actual distribution of hits
                if(cached_keys_x[current_word_x] < cached_keys_y[0])
                {
                    if(current_word_x < blockDim.x-1 && cached_keys_x[current_word_x+1] > cached_keys_x[current_word_x]
                        && cached_keys_x[current_word_x+1] >= cached_keys_y[0]){
                            word_y_id = current_word_y_id;
                            if(threadIdx.x == 0 && messageID<messageMAX) messages[messageID++] = (((uint64_t)0 << 32) | word_y_id); //printf("%d Strong y increment %d!\n", messageID++, word_y_id);
                        }
                    if(threadIdx.x == 0 && messageID<messageMAX) messages[messageID++] = (((uint64_t)1 << 32) | word_x_id); //printf("%d Increment x %d!\n", messageID++, word_x_id);
                    break;
                }
                if(cached_keys_x[current_word_x] > cached_keys_y[0]){ // this condition has increased performance
                    word_y_id = current_word_y_id;
                    if(threadIdx.x == 0 && messageID<messageMAX) messages[messageID++] = (((uint64_t)4 << 32) | word_y_id); //printf("%d Major Strong y increment %d!\n", messageID++, word_x_id);
                }

                if(threadIdx.x == 0 && messageID<messageMAX) messages[messageID++] = (((uint64_t)2 << 32) | current_word_y_id); //printf("%d Basic y increment %d!\n", messageID++, current_word_y_id);
                current_word_y_id += blockDim.x;
                */

                /*
                // A bit faster because of lesser ifs (equivalent to the next one) [T2]
                if(cached_keys_x[current_word_x] <= cached_keys_y[0])
                    break;
                current_word_y_id += blockDim.x;
                */

                /*
                // Slow but perfect [T1]
                if(hit == blockDim.x)
                    current_word_y_id += blockDim.x;
                else if(hit < blockDim.x && cached_keys_x[current_word_x] > cached_keys_y[0])
                    current_word_y_id += blockDim.x;
                else if(cached_keys_x[current_word_x] > cached_keys_y[0])
                {
                    word_y_id += blockDim.x;
                    current_word_y_id = word_y_id; 
                }
                else
                    break;
                */ 

    


    /*
    int32_t word_x_id = blockIdx.x;
    int32_t word_y_id = 0;
    int32_t current_word_y_id = word_y_id;

    int32_t accumulated = blockIdx.x * mem_block;
    int32_t highest_winner = 0, total_hits = 0;

    // Loop to generate hits
    while(word_x_id < limit_x && word_y_id + threadIdx.x < limit_y) //limit_x)
    {

        if(positions_x[word_x_id] == 0xFFFFFFFF || positions_y[word_y_id] == 0xFFFFFFFF) break;

        // One read to check whether we have to update x
        if(hashes_x[word_x_id] < hashes_y[word_y_id]){
            word_x_id += gridDim.x;
        
        }else if(hashes_x[word_x_id] > hashes_y[word_y_id]){
            //++word_y_id;
            word_y_id += max(1, highest_winner);
            //word_y_id += blockDim.x;
        }else{

            current_word_y_id = word_y_id;

            while(hashes_x[word_x_id] == hashes_y[current_word_y_id]){
            
                if(hashes_x[word_x_id] == hashes_y[current_word_y_id+threadIdx.x]){
                    // Perform the write
                    hit_write_section[accumulated + threadIdx.x] = (((uint64_t)ref_len + (uint64_t)positions_x[word_x_id] - (uint64_t)positions_y[current_word_y_id+threadIdx.x]) << 32 ) + (uint64_t)positions_x[word_x_id];
                    highest_winner = threadIdx.x + 1;
                }else{
                    highest_winner = 0;
                }

                // Distribute highest threadIdx.x among threads
                
                for (int offset = 16; offset > 0; offset = offset >> 1)
                    highest_winner += __shfl_down_sync(0xFFFFFFFF, highest_winner, offset);

                highest_winner = __shfl_sync(0xFFFFFFFF, highest_winner, 0);
                    
                accumulated += highest_winner;
                total_hits += highest_winner;

                if(accumulated >= (blockIdx.x+1) * mem_block) { *error = -1; return; }

                current_word_y_id += blockDim.x;
            }
            word_x_id += gridDim.x;
        }

        //if(threadIdx.x == 0) printf("Thr Id %d from block %d ::: Outcome: %d -> indices [%d] [%d] (hits so far %d)\n", threadIdx.x, blockIdx.x, highest_winner, word_x_id, word_y_id, total_hits);

    }

    if(threadIdx.x == 0) hits_log[blockIdx.x] = (uint32_t) total_hits;

    */
    

}

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
        //if(threadIdx.x==0) printf("Accesing %d\n", thre_pos_x);
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

        /*
        // This section skips filtered hits
        if(hp1 == 0xFFFFFFFF && hp2_shared[curr_frag] == 0xFFFFFFFF){
            l_o[curr_frag] = 0;
            r_o[curr_frag] = 0;
            ++real_id; ++curr_frag;  continue; 
        }
        */

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


__global__ void kernel_index_global32(uint64_t * hashes, uint32_t * positions, const char * sequence, uint32_t offset, uint32_t seq_lim) {

    uint64_t hash = 0, k = 0;
    uint64_t bad = 0;
    
    if(threadIdx.x + 32 + blockIdx.x * blockDim.x > seq_lim) return;

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
        if(c == 'N') bad = 0xFFFFFFFFFFFFFFFF;
        

            
        // Pure arithmetic
        /*
        hash = hash << 2;
        hash += (c & 6) >> 1;
        if(c == 'N') bad = 0;
        */
        

    }
    // [0 - 32] * [0-N]* [32]
    hashes[threadIdx.x + blockIdx.x * blockDim.x] = hash | bad;//& bad; // reconsider changing to | ~bad
    positions[threadIdx.x + blockIdx.x * blockDim.x] = (threadIdx.x + blockIdx.x * blockDim.x + offset) | bad;
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

