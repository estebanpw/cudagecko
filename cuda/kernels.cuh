#include <cuda.h>
#include <stdio.h>
#include <inttypes.h>

#define WARP_SIZE 32
#define KMER_SIZE 32
#define KMERS_PER_THREAD 3
#define BYTES_PER_REGISTER 4
#define ULLI unsigned long long int
#define LLI  long long int

__global__ void kernel_frags_forward_register(uint64_t * h_p1, uint64_t * h_p2, uint64_t * left_offset, uint64_t * right_offset, const char * seq_x, const char * seq_y, uint64_t query_len, uint64_t ref_len, uint64_t x_seq_off, uint64_t y_seq_off, uint64_t x_lim, uint64_t y_lim);

__global__ void kernel_register_fast_hash_rotational(uint64_t * hashes, uint64_t * positions, const char * sequence, uint64_t offset);

__global__ void kernel_index_global32(uint64_t * hashes, uint64_t * positions, const char * sequence, uint64_t offset);