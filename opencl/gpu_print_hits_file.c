
// Standard utilities and common systems includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "structs.h"


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char ** av)
{


    hitGPU hgpu;
    FILE * hits = fopen(av[1], "rb");
    if(hits == NULL) { fprintf(stderr, "Could not read input\n"); exit(-1); }

    uint64_t i = 0, v;
    while(!feof(hits)){
        v = fread(&hgpu, sizeof(hitGPU), 1, hits);
        if(v == 0) v = 0; // So that gcc doesnt complain
        
        if(i % atoi(av[2]) == 0){
            fprintf(stdout, "[%"PRIu64"] @%"PRId64" (%"PRIu64", %"PRIu64") from (seq X: %"PRIu64", seq Y: %"PRIu64")\n", i, (int64_t) hgpu.pos_x - (int64_t) hgpu.pos_y, hgpu.pos_x, hgpu.pos_y, hgpu.seq_x, hgpu.seq_y);
            getchar();
        }
        ++i;
    }

    fclose(hits);

    return 0;
}