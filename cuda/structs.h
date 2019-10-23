#include <cuda.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

// Warning: Does not prevent double evaluation
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define ULLI unsigned long long int

typedef struct word{
    uint64_t hash;
    uint64_t pos;
} Word;