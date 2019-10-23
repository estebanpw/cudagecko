

typedef struct param_sort{
    ulong N;
    ulong comparators_per_wi;
    ulong stage;
    ulong step;
    ulong kmer_size;
} parameter_sort;

typedef struct {
    ulong h;        // The actual word
    ulong pos;      // The position in the sequence file
} word_hash_GPU;




__kernel void kernel_words_sort(__global word_hash_GPU * in, __global parameter_sort * params)
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
            

            int is_j_smaller_than_i = (int) (((ulong)in[j].h) < ((ulong)in[i].h));
            
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
                

                
                if(in[i].h != in[j].h){ // XOR with its own will produce 0
                    in[i].h = in[i].h ^ in[j].h;
                    in[j].h = in[j].h ^ in[i].h;
                    in[i].h = in[i].h ^ in[j].h;
                }

                if(in[i].pos != in[j].pos){
                    in[i].pos = in[i].pos ^ in[j].pos;
                    in[j].pos = in[j].pos ^ in[i].pos;
                    in[i].pos = in[i].pos ^ in[j].pos;
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

                
                if(in[i].h != in[j].h){
                    in[i].h = in[i].h ^ in[j].h;
                    in[j].h = in[j].h ^ in[i].h;
                    in[i].h = in[i].h ^ in[j].h;
                }

                if(in[i].pos != in[j].pos){
                    in[i].pos = in[i].pos ^ in[j].pos;
                    in[j].pos = in[j].pos ^ in[i].pos;
                    in[i].pos = in[i].pos ^ in[j].pos;
                }
                
            }
        }
    }
}