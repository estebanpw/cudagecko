#define POINT 4

typedef struct {
    ulong pos_x;
    ulong pos_y;
    ulong length;
    ulong identities;
    // Notice that seq_x and seq_y can be retrieved from the hits when CPU is writing them
} reduced_frag;


__kernel void kernel_craze(__global char * seq_x, __global char * seq_y, __global reduced_frag * rf)
{

    // Each work item is a hit
    ulong id = get_global_id(0);

    // Make 16 comparisons
    ulong i;
    int abort = 0;
    for(i=id; i<id+32; i++){
        if(seq_x[i] != seq_y[i]){ abort = 1; break; }
    }

}