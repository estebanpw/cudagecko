#include <cuda.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#define BUFFER_SIZE 2048
#define CORES_PER_COMPUTE_UNIT 32
#define KMER_SIZE 32

// Warning: Does not prevent double evaluation
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct word{
    uint64_t hash;
    uint64_t pos;
} Word;

typedef struct hit{
    uint64_t p1;
    uint64_t p2;
} Hit;