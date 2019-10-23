
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
void showWordGPU(wordGPU * w, char * ws) {
    char Alf[] = {'A', 'C', 'G', 'T'};
    int i;
    int wsize = 8;
    unsigned char c;
    for (i = 0; i < wsize; i++) {
        c = w->b[i];
        c = c >> 6;
        ws[4 * i] = Alf[(int) c];
        c = w->b[i];
        c = c << 2;
        c = c >> 6;
        ws[4 * i + 1] = Alf[(int) c];
        c = w->b[i];
        c = c << 4;
        c = c >> 6;
        ws[4 * i + 2] = Alf[(int) c];
        c = w->b[i];
        c = c << 6;
        c = c >> 6;
        ws[4 * i + 3] = Alf[(int) c];
    }
    ws[33] = '\0';
}


int main(int argc, char ** av)
{


    wordGPU wgpu;
    FILE * words = fopen(av[1], "rb");
    if(words == NULL) { fprintf(stderr, "Could not read input\n"); exit(-1); }

    uint64_t i = 0, v;
    char buff[33];
    while(!feof(words)){
        v = fread(&wgpu, sizeof(wordGPU), 1, words);
        if(v == 0) v = 0; // So that gcc doesnt complain
        showWordGPU(&wgpu, buff);
        if(i % atoi(av[2]) == 0){
            fprintf(stdout, "[%"PRIu64"] %s @%"PRIu64"\n", i, buff, wgpu.pos);
            getchar();
        }
        ++i;
    }

    fclose(words);

    return 0;
}