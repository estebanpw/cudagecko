#include "structs.h"



__global__ void kernel_register_fast_hash_rotational(uint64_t * hashes, uint64_t * positions, const char * sequence, ULLI offset);