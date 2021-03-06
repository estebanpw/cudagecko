

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


typedef struct {
    ulong pos_x;
    ulong pos_y;
} hit_advanced;



__kernel void kernel_sort(__global hit_advanced * in, __global parameter_sort * params)
{
    ulong stage = params->stage;
    ulong step = params->step;
    ulong N = params->N;


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
            

            int is_j_smaller_than_i = (int) (((long)in[j].pos_x-(long)in[j].pos_y) < ((long)in[i].pos_x-(long)in[i].pos_y));
            
            if(is_j_smaller_than_i == 1 && dir == 0){
                // Store
                //swap_wordGPUs(&in[i], &in[j]);

                /*
                X := X XOR Y
                Y := Y XOR X
                X := X XOR Y
                */
                /*
                ulong aux;
                aux = in[i].pos_x;
                in[i].pos_x = in[j].pos_x;
                in[j].pos_x = aux;

                aux = in[i].pos_y;
                in[i].pos_y = in[j].pos_y;
                in[j].pos_y = aux;
                */

                
                if(in[i].pos_x != in[j].pos_x){ // XOR with its own will produce 0
                    in[i].pos_x = in[i].pos_x ^ in[j].pos_x;
                    in[j].pos_x = in[j].pos_x ^ in[i].pos_x;
                    in[i].pos_x = in[i].pos_x ^ in[j].pos_x;
                }

                if(in[i].pos_y != in[j].pos_y){
                    in[i].pos_y = in[i].pos_y ^ in[j].pos_y;
                    in[j].pos_y = in[j].pos_y ^ in[i].pos_y;
                    in[i].pos_y = in[i].pos_y ^ in[j].pos_y;
                }

                
                

                //in[i] = jkey;
                //in[j] = ikey;
            }else if(is_j_smaller_than_i != 1 && dir == 1){
                //swap_wordGPUs(&in[i], &in[j]);
                //in[i] = jkey;
                //in[j] = ikey;

                /*
                ulong aux;
                aux = in[i].pos_x;
                in[i].pos_x = in[j].pos_x;
                in[j].pos_x = aux;

                aux = in[i].pos_y;
                in[i].pos_y = in[j].pos_y;
                in[j].pos_y = aux;
                */

                
                if(in[i].pos_x != in[j].pos_x){
                    in[i].pos_x = in[i].pos_x ^ in[j].pos_x;
                    in[j].pos_x = in[j].pos_x ^ in[i].pos_x;
                    in[i].pos_x = in[i].pos_x ^ in[j].pos_x;
                }

                if(in[i].pos_y != in[j].pos_y){
                    in[i].pos_y = in[i].pos_y ^ in[j].pos_y;
                    in[j].pos_y = in[j].pos_y ^ in[i].pos_y;
                    in[i].pos_y = in[i].pos_y ^ in[j].pos_y;
                }

                
                
            }
        }
    }
}