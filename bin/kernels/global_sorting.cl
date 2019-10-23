

typedef struct param_sort{
    ulong N;
    ulong comparators_per_wi;
    ulong stage;
    ulong step;
    ulong kmer_size;
} parameter_sort;

typedef struct {
    unsigned char b[8];  // The actual word
    ulong pos;           // The position in the sequence file
} wordGPU;

int wordcmp(__global wordGPU * w1, __global wordGPU * w2) {
    int i = 0, limit;
    int n = 32;

    if (n % 4 != 0) {
        w1->b[n / 4] = w1->b[n / 4] >> (2 * (3 - ((n - 1) % 4)));
        w1->b[n / 4] = w1->b[n / 4] << (2 * (3 - ((n - 1) % 4)));
        w2->b[n / 4] = w2->b[n / 4] >> (2 * (3 - ((n - 1) % 4)));
        w2->b[n / 4] = w2->b[n / 4] << (2 * (3 - ((n - 1) % 4)));
        limit = (n / 4) + 1;
    } else {
        limit = n / 4;
    }

    for (i = 0; i < limit; i++) {
        if (w1->b[i] < w2->b[i]) return -1;
        if (w1->b[i] > w2->b[i]) return +1;
    }
    return 0;
}

void copy_wordGPU(__global wordGPU * origin, __global wordGPU * destination){
    int i;
    for(i=0;i<8;i++) destination->b[i] = origin->b[i];
    destination->pos = origin->pos;
}

void swap_wordGPUs(__global wordGPU * w1, __global wordGPU * w2){
    
    int i;
    unsigned char byte;
    for(i=0; i<8; i++){
        byte = w1->b[i];
        w1->b[i] = w2->b[i];
        w2->b[i] = byte;
    }


    w1->pos = w1->pos ^ w2->pos;
    w2->pos = w2->pos ^ w1->pos;
    w1->pos = w1->pos ^ w2->pos;
}


__kernel void kernel_sort(__global wordGPU * in, __global parameter_sort * params)
{
    ulong stage = params->stage;
    ulong step = params->step;
    ulong N = params->N;
    ulong kmer = params->kmer_size;


    ulong k;
    for(k=0; k<params->comparators_per_wi; k++){

        ulong i = get_global_id(0) * params->comparators_per_wi + k;

        ulong j = i ^ (1 << (stage-1)); // Generates sibling

        int dir;
        if((i/(2 << (step-1)) % 2) == 0) dir = 0; else dir = 1; // 0 is ascending, 1 is descending

        if(j < N && i < j){
            // Load values at I and J
            //unsigned char * ikey = in[i].b;
            //unsigned char * jkey = in[j].b;

            // Compare j with i
            //bool is_j_smaller_than_i = (jkey <= ikey);

            int is_j_smaller_than_i = wordcmp(&in[j], &in[i]);
            
            if(is_j_smaller_than_i == -1 && dir == 0){
                // Store
                swap_wordGPUs(&in[i], &in[j]);

                //in[i] = jkey;
                //in[j] = ikey;
            }else if(is_j_smaller_than_i != -1 && dir == 1){
                swap_wordGPUs(&in[i], &in[j]);
                //in[i] = jkey;
                //in[j] = ikey;
            }
        }
    }
}