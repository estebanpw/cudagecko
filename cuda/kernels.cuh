#include <cuda.h>
#include <inttypes.h>

#define WARP_SIZE 32
#define KMER_SIZE 32
#define KMERS_PER_THREAD 3
#define BYTES_PER_REGISTER 4
#define ULLI unsigned long long int


__global__ void kernel_register_fast_hash_rotational(uint64_t * hashes, uint64_t * positions, const char * sequence, ULLI offset);