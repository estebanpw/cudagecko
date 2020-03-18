#include <cuda.h>
#include <stdio.h>
#include <inttypes.h>

#define WARP_SIZE 32
#define KMER_SIZE 32
#define KMERS_PER_THREAD 3
#define BYTES_PER_REGISTER 4
#define ULLI unsigned long long int
#define LLI  long long int

// Placeholder hit - I had some trouble with compiling order

typedef struct ph_hit{
    uint32_t p1;
    uint32_t p2;
} ph_Hit;

__global__ void kernel_frags_forward_register(uint32_t * h_p1, uint32_t * h_p2, uint32_t * left_offset, uint32_t * right_offset, const char * seq_x, const char * seq_y, uint32_t query_len, uint32_t ref_len, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t x_lim, uint32_t y_lim);

__global__ void kernel_frags_reverse_register(uint32_t * h_p1, uint32_t * h_p2, uint32_t * left_offset, uint32_t * right_offset, const char * seq_x, const char * seq_y, uint32_t query_len, uint32_t ref_len, uint32_t x_seq_off, uint32_t y_seq_off, uint32_t x_lim, uint32_t y_lim);

__global__ void kernel_register_fast_hash_rotational(uint64_t * hashes, uint64_t * positions, const char * sequence, uint64_t offset);

__global__ void kernel_index_global32(uint64_t * hashes, uint32_t * positions, const char * sequence, uint32_t offset);

__global__ void kernel_reverse_complement(const char * sequence, char * reverse_sequence, uint32_t seq_len);

__global__ void kernel_generate_hits(uint64_t * hash_x, uint64_t * hash_y, uint32_t * pos_x, uint32_t * pos_y, uint64_t * diagonals, ph_Hit * h);