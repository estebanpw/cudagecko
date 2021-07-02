#include <cuda.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "vcl/vectorclass.h"

#define BUFFER_SIZE 2048
#define CORES_PER_COMPUTE_UNIT 32
#define KMER_SIZE 32
#define MIN_P_IDENT 0.8 // THIS MUST BE CHANGED IN KERNELS.CU as well
#define OUTPUT_P_IDENT 80.00

// Warning: Does not prevent double evaluation
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct word{
    uint64_t hash;
    uint32_t pos;
} Word;

typedef struct hit{
    uint32_t p1;
    uint32_t p2;
} Hit;
