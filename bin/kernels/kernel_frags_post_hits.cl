#define POINT 4


typedef struct {
    ulong kmer_size;
    ulong size_x;
    ulong size_y;
    ulong t_hits;
} parameter_frags;

typedef struct {
    unsigned char b[8];  // The actual word
    ulong pos;           // The position in the sequence file
} wordGPU;

typedef struct {
    ulong pos_x;
    ulong pos_y;
    ulong seq_x;
    ulong seq_y;
} hitGPU;

typedef struct {
    ulong pos_x;
    ulong pos_y;
    ulong length;
    ulong identities;
    // Notice that seq_x and seq_y can be retrieved from the hits when CPU is writing them
} reduced_frag;

typedef struct {
    long diag;
    //Ocurrence position in sequence X
    ulong posX;
    //Ocurrence position in sequence Y
    ulong posY;
    //For multiple sequence files this var
    //reflects in what sequence of X file
    //occurs the word
    ulong seqX;
    //For multiple sequence files this var
    //reflects in what sequence of Y file
    //occurs the word
    ulong seqY;
} hit;

__kernel void kernel_frags(__global char * seq_x, __global char * seq_y, __global hitGPU * in, __global reduced_frag * rf, __global parameter_frags * p)
{

    // Each work item is a hit
    ulong id = get_global_id(0);

    

    if(id < p->t_hits && in[id].pos_x != 0xFFFFFFFFFFFFFFFF && in[id].pos_y != 0xFFFFFFFFFFFFFFFF){

        long kmer_size = (long) p->kmer_size;
        long score = kmer_size * POINT;
        long score_right = score;
        long score_left = score;
        long identities = (long) kmer_size;
        long last_identities = identities;
        
        long max_len_x = p->size_x;
        long max_len_y = p->size_y;
        
        long start_x = in[id].pos_x;
        long end_x = start_x + kmer_size;
        long start_y = in[id].pos_y;
        long end_y = start_y + kmer_size;

        long curr_x = end_x;
        long curr_y = end_y;
        long curr_score = score_right;

        // To the right
        
        while(curr_score > 0){

            ++curr_x;
            ++curr_y;
            
            if(curr_x == max_len_x || curr_y == max_len_y) break;
            if(seq_x[curr_x] == '*') break; // Remember to have sequences converted
            if(seq_y[curr_y] == '*') break;

            if(seq_x[curr_x] == seq_y[curr_y] && seq_x[curr_x] != 'N'){
                // It is a match
                curr_score += POINT;
                ++identities;
                if(curr_score >= score_right){
                    // Update best position
                    score_right = curr_score;
                    end_x = curr_x;
                    end_y = curr_y;
                    last_identities = identities;
                }
            }else{
                // It is a missmatch or N
                curr_score -= POINT;
            }
        }

        // To the left

        curr_x = start_x;
        curr_y = start_y;
        curr_score = score_left;
        identities = last_identities;


        while(curr_score > 0){

            --curr_x;
            --curr_y;
            
            if(curr_x < 0 || curr_y < 0) break;
            if(seq_x[curr_x] == '*') break; // Remember to have sequences converted
            if(seq_y[curr_y] == '*') break;

            if(seq_x[curr_x] == seq_y[curr_y] && seq_x[curr_x] != 'N'){
                // It is a match
                curr_score += POINT;
                ++identities;
                if(curr_score >= score_left){
                    // Update best position
                    score_left = curr_score;
                    start_x = curr_x;
                    start_y = curr_y;
                    last_identities = identities;
                }
            }else{
                // It is a missmatch or N
                curr_score -= POINT;
            }
        }
        
        // Save frag
        rf[id].pos_x = (ulong) start_x;
        rf[id].pos_y = (ulong) start_y;
        rf[id].length = (ulong) end_x - start_x;
        rf[id].identities = (ulong) last_identities;

    }else{
        rf[id].pos_x = 0xFFFFFFFFFFFFFFFF;
    }
    

}