#include <cuda.h>
#include <stdio.h>
#include <inttypes.h>

#define WARP_SIZE 32
#define KMER_SIZE 32
#define KMERS_PER_THREAD 3
#define BYTES_PER_REGISTER 4
#define ULLI unsigned long long int
#define LLI  long long int

__global__ void kernel_filter_hits_parallel_shared(uint64_t * diagonals_merged, uint32_t ref_len, uint32_t n_hits_found);

__global__ void kernel_filter_hits_parallel(uint64_t * diagonals_merged, uint32_t ref_len, uint32_t n_hits_found);

__global__ void kernel_filter_hits(uint64_t * diagonals_merged, uint32_t ref_len, uint32_t n_hits_found);

__global__ void kernel_find_leftmost_items(uint64_t * hashes_x, uint32_t * pos_x, uint64_t * hashes_y, uint32_t * pos_y, uint32_t limit_x, uint32_t limit_y);

__global__ void kernel_compact_hits(uint64_t * ptr_device_diagonals, uint32_t * ptr_hits_log, uint32_t * ptr_accum_log, uint32_t mem_block, uint64_t * ptr_copy_place_diagonals, uint32_t offset_remover);

__global__ void kernel_hits_load_balancing(uint64_t * hashes_x, uint64_t * hashes_y, uint32_t * positions_x, uint32_t * positions_y, uint64_t * hit_write_section, int32_t mem_block, uint32_t limit_x, uint32_t limit_y, int32_t * error, uint32_t ref_len, uint32_t * hits_log, int32_t * atomic_distributer, uint64_t * auxiliary_hit_memory, uint32_t extra_large_memory_block, uint32_t max_extra_sections, uint32_t * hits_log_extra, int32_t * atomic_queue);

__global__ void kernel_hits(uint64_t * hashes_x, uint64_t * hashes_y, uint32_t * positions_x, uint32_t * positions_y, uint64_t * hit_write_section, int32_t mem_block, uint32_t limit_x, uint32_t limit_y, int32_t * error, uint32_t ref_len, uint32_t * hits_log, int32_t * atomic_distributer, uint64_t * auxiliary_hit_memory, uint32_t extra_large_memory_block, uint32_t max_extra_sections, uint32_t * hits_log_extra);//, uint64_t * messages);

__global__ void kernel_frags_forward_register(uint32_t * h_p1, uint32_t * h_p2, uint32_t * left_offset, uint32_t * right_offset, const char * seq_x, const char * seq_y, uint32_t query_len, uint32_t ref_len, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t x_lim, uint32_t y_lim, uint32_t n_hits_kept, uint32_t n_frags_per_block);

__global__ void kernel_frags_forward_per_thread(uint32_t * h_p1, uint32_t * h_p2, uint32_t * left_offset, uint32_t * right_offset, const char * seq_x, const char * seq_y, uint32_t query_len, uint32_t ref_len, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t x_lim, uint32_t y_lim, uint32_t n_hits);

__global__ void kernel_frags_reverse_register(uint32_t * h_p1, uint32_t * h_p2, uint32_t * left_offset, uint32_t * right_offset, const char * seq_x, const char * seq_y, uint32_t query_len, uint32_t ref_len, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t x_lim, uint32_t y_lim, uint32_t n_hits_kept, uint32_t n_frags_per_block);

__global__ void kernel_register_fast_hash_rotational(uint64_t * hashes, uint64_t * positions, const char * sequence, uint64_t offset);

__global__ void kernel_index_global32(uint64_t * hashes, uint32_t * positions, const char * sequence, uint32_t offset, uint32_t seq_lim);

__global__ void kernel_index_global32_advanced(uint64_t * hashes, uint32_t * positions, const uchar4 * sequence, uint32_t offset);

__global__ void kernel_reverse_complement(const char * sequence, char * reverse_sequence, uint32_t seq_len);

