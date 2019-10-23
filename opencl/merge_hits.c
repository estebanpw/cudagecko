
// Standard utilities and common systems includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "structs.h"


#define BUFFER_SIZE 2048
#define HITS_BUFFER 400000  // times 40 and for 300 frames aprox it is 4.8 GB
#define MAX_KERNEL_SIZE BUFFER_SIZE*100
#define CORES_PER_COMPUTE_UNIT 32
#define MIN(a,b) (((a)<(b))?(a):(b))

char * get_dirname(char * path);
char * get_basename(char * path);
void init_args(int argc, char ** av, FILE ** hits_file, FILE ** hits_file_headers, cl_uint * selected_device, FILE ** hits_filtered, ulong * fixed_power_of_two, cl_ulong * queue_size);

typedef struct{
    hit * h;
    uint64_t i;
    uint64_t curr_pos;
    uint64_t limit;

} buffer_hit;

typedef struct{
    hit * h;
} hit_pointer;

class fixed_priority_queue{
    private:
    uint64_t size, t_inserted, t_liberated;
    FILE * writer;
    hit_pointer * queue;
    void reset() { memset(this->queue, 0x0, this->size); this->t_inserted = 0; }
    

    public:

    fixed_priority_queue(uint64_t size, FILE * write){
        this->queue = (hit_pointer *) calloc(size, sizeof(hit_pointer));
        if(this->queue == NULL){ fprintf(stderr, "Could not allocate priority queue\n"); exit(-1); }
        this->size = size;
        this->writer = write;
        this->t_inserted = 0;
        this->t_liberated = 0;
    }
    void insert(hit * h){
        int64_t L = 0;
        int64_t R = (int64_t) this->t_inserted - 1;
        int64_t m = 0;
        while(L <= R){
            m = (L + R) / 2;
            if(this->queue[m].h == NULL){ 
                R = m - 1; 
                //printf("m: %"PRId64" L: %"PRId64", R: %"PRId64" \n", m, L, R);
                continue; 
            } else { 
                //printf("m: %"PRId64" L: %"PRId64", R: %"PRId64" (comparing: %"PRId64" vs %"PRId64") \n", m, L, R, this->queue[m].h->diag, h->diag);
            }

            if (this->queue[m].h->diag < h->diag){
                L = m + 1;
                m++;
            }else if (this->queue[m].h->diag > h->diag){
                R = m - 1;
            }else if(this->queue[m].h->diag == h->diag){
                break;
            }
        }
        //printf("adding %"PRId64" total inserted: %"PRId64" at m: %"PRId64"\n", h->diag, t_inserted+1, m);
        memmove(&this->queue[m+1], &this->queue[m], (t_inserted-m)*sizeof(hit_pointer));
        this->queue[m].h = h;
        this->t_inserted++;
        if(this->t_inserted == this->size){ fprintf(stderr, "Reached max merge list\n"); exit(-1); }
        /*
        uint64_t r;
        for(r=0; r<10; r++){
            if(this->queue[r].h == NULL) printf("at %"PRIu64" -> diag: NULL\n", r);
            else printf("at %"PRIu64" -> diag: %"PRId64"\n", r, this->queue[r].h->diag);
        }
        */
    }

    void liberate_hits(int64_t boundary){ 
        uint64_t i = 0;
        //fprintf(stdout, "The head is \t%"PRId64"\n", this->queue[0].h->diag);
        while(this->queue[i].h->diag <= boundary && i<this->size){
            fwrite(this->queue[i].h, sizeof(hit), 1, this->writer);
            ++i;
        }
        t_liberated += i;
        if(t_liberated % 10000 < 1000) printf("Liberating %"PRIu64" hits. Curr load=%"PRIu64". T-liber = %"PRIu64"\n", i, t_inserted-i, t_liberated);
        //getchar();
        memmove(&this->queue[0], &this->queue[i], (t_inserted-i)*sizeof(hit_pointer));
        t_inserted -= i;
    }

};

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
    cl_uint selected_device = 0;
    
    ////////////////////////////////////////////////////////////////////////////////
    // Get info of devices
    ////////////////////////////////////////////////////////////////////////////////
    cl_platform_id platform_id = NULL;
    cl_device_id * devices = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_uint compute_units;
    char device_info[BUFFER_SIZE]; device_info[0] = '\0';
    cl_bool device_available;
    cl_ulong device_RAM, global_device_RAM[MAX_DEVICES];
    cl_ulong queue_size = 10000, fixed_power_of_two = 0;
    size_t work_group_size[3], work_group_size_global;
    cl_int ret;
    char * path_kernels = get_dirname(argv[0]);
    fprintf(stdout, "[INFO] Working on directory %s\n", path_kernels);

    FILE * hits_file = NULL, * hits_filtered = NULL, * hits_file_headers = NULL;

    init_args(argc, argv, &hits_file, &hits_file_headers, &selected_device, &hits_filtered, &fixed_power_of_two, &queue_size);

    if(CL_SUCCESS != (ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms))){ fprintf(stderr, "Failed to get platforms\n"); exit(-1); }
    fprintf(stdout, "Detected %d platform(s)\n", ret_num_platforms);
    

    // Query how many devices there are
    if(CL_SUCCESS != (ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, 0, NULL, &ret_num_devices))){ fprintf(stderr, "Failed to query number of devices\n"); exit(-1); }
    if ((devices = (cl_device_id *) malloc(sizeof(cl_device_id) * ret_num_devices)) == NULL){ fprintf(stderr, "Could not allocate devices information\n"); exit(-1); }

    // Query information about each device
    if(CL_SUCCESS != (ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, ret_num_devices, devices, &ret_num_devices))){ fprintf(stderr, "Failed to get devices information\n"); exit(-1); }
    
    fprintf(stdout, "Found %d device(s)\n", ret_num_devices);
    fprintf(stdout, "These are:\n");
    cl_uint i;
    for(i=0; i<ret_num_devices; i++){

        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, BUFFER_SIZE, device_info, NULL))) fprintf(stderr, "Failed to get device name %d\n", i);
        fprintf(stdout, "\tDevice [%d]: %s\n", i, device_info);
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_VERSION, BUFFER_SIZE, device_info, NULL))) fprintf(stderr, "Failed to get device profile %d\n", i);
        fprintf(stdout, "\t\tProfile      : %s\n", device_info);
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_AVAILABLE, sizeof(cl_bool), &device_available, NULL))){ fprintf(stderr, "Failed to get device availability %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tis available?: %d\n", (int)device_available);
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &global_device_RAM[i], NULL))){ fprintf(stderr, "Failed to get device global memory %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tGlobal mem   : %"PRIu64" (%"PRIu64" MB)\n", (uint64_t) global_device_RAM[i], (uint64_t) global_device_RAM[i] / (1024*1024));
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &device_RAM, NULL))){ fprintf(stderr, "Failed to get device local memory %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tLocal mem    : %"PRIu64" (%"PRIu64" KB)\n", (uint64_t) device_RAM, (uint64_t) device_RAM / (1024));
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &compute_units, NULL))){ fprintf(stderr, "Failed to get device local memory %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tCompute units: %"PRIu64"\n", (uint64_t) compute_units);
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size_global, NULL))){ fprintf(stderr, "Failed to get device global work items size %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tMax work group size: %zu\n", work_group_size_global);
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, 3*sizeof(size_t), &work_group_size, NULL))){ fprintf(stderr, "Failed to get device work items size %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tWork size items: (%zu, %zu, %zu)\n", work_group_size[0], work_group_size[1], work_group_size[2]);
    }
    //selected_device = 3;
    fprintf(stdout, "[INFO] Using device %d\n", selected_device);
    
    

    ////////////////////////////////////////////////////////////////////////////////
    // Calculate how many words we can load at once
    ////////////////////////////////////////////////////////////////////////////////

    ulong N_at_once = (global_device_RAM[selected_device]-100*1024*1024) / sizeof(hitGPU); // Minus 100 MB
    
    // Find closest power of two to use
    fprintf(stdout, "[INFO] Maximum hits to load at once: %"PRIu64"\n", (uint64_t) N_at_once);
    if(fixed_power_of_two == 0){
        N_at_once = (ulong) floor(log2((float) N_at_once));
        N_at_once = (ulong) pow(2, (float) N_at_once);
    }else{
        N_at_once = (ulong) pow(2, (float) fixed_power_of_two);
        if(N_at_once * sizeof(hitGPU) > (global_device_RAM[selected_device]-100*1024*1024)){ fprintf(stderr, "Exceeding GPU RAM limits\n"); exit(-1); }
    }
    fprintf(stdout, "[INFO] Reading frames of %"PRIu64" hits, i.e. %"PRIu64" MB.\n", (uint64_t) N_at_once, (uint64_t) ((N_at_once * sizeof(hitGPU)) / (1024*1024)));

    // Determine number of frames
    uint64_t n_frames = 0;
    char c = 'z';
    while(!feof(hits_file_headers)){
        c = fgetc(hits_file_headers);
        if(c == '\n') n_frames++;
    }
    rewind(hits_file_headers);

    fprintf(stdout, "[INFO] Detected %"PRIu64" frames\n", n_frames);

    // Determine frame distance and load frames
    uint64_t * frame_positions = (uint64_t *) malloc(n_frames * sizeof(uint64_t));
    uint64_t min_distance_between_frames = 0xFFFFFFFFFFFFFFFF;
    if(1 != fscanf(hits_file_headers, "%"PRIu64, &frame_positions[0])){ fprintf(stderr, "Could not read header\n"); exit(-1); }
    for(i=1; i<n_frames; i++){
        if(1 != fscanf(hits_file_headers, "%"PRIu64, &frame_positions[i])){ fprintf(stderr, "Could not read header\n"); exit(-1); }
        fprintf(stdout, "Frames: %"PRIu64"\n", frame_positions[i]);
        if((frame_positions[i] - frame_positions[i-1]) < min_distance_between_frames) min_distance_between_frames = (frame_positions[i] - frame_positions[i-1]);
    }

    // Allocate buffers
    uint64_t buffer_size = MIN(BUFFER_SIZE, min_distance_between_frames);
    buffer_hit * buffers = (buffer_hit *) malloc(n_frames * sizeof(buffer_hit));
    if(buffers == NULL) { fprintf(stderr, "Could not allocate first level loop of buffer frames\n"); exit(-1); }
    for(i=0; i<n_frames; i++){
        buffers[i].h = (hit *) malloc(buffer_size * sizeof(hit));
        if(buffers[i].h == NULL) { fprintf(stderr, "Could not allocate second level loop of buffer frames\n"); exit(-1); }
    }

    fseek(hits_file, 0, SEEK_END);
    uint64_t t_hits = ftell(hits_file)/sizeof(hit);
    rewind(hits_file);

    // Read items
    uint64_t read;
    for(i=0; i<n_frames-1; i++){
        buffers[i].curr_pos = frame_positions[i];
        buffers[i].i = 0;
        buffers[i].limit = frame_positions[i+1];
        fseek(hits_file, buffers[i].curr_pos * sizeof(hit), SEEK_SET);
        read = fread(buffers[i].h, sizeof(hit), buffer_size, hits_file);
        if(read == 0) fprintf(stdout, "Read zero items\n");
    }
    buffers[n_frames-1].curr_pos = frame_positions[n_frames-1];
    fseek(hits_file, buffers[i].curr_pos * sizeof(hit), SEEK_SET);
    read = fread(buffers[i].h, sizeof(hit), buffer_size, hits_file);
    buffers[n_frames-1].limit = t_hits;
    buffers[n_frames-1].i = 0;
    
    // Create queue
    fixed_priority_queue * fpq = new fixed_priority_queue(queue_size, hits_filtered);

    // Loop
    uint64_t buffers_finished = 0;
    while( buffers_finished < n_frames ){

        int64_t minimum_value = 9223372036854775807;
        int64_t remember_queue = -1;
        // Check if a buffer is done 
        uint64_t t;
        for(t=0; t<n_frames; t++){
            // If we reached end of buffer and there is still more to read
            if((buffers[t].i + buffers[t].curr_pos) == buffers[t].limit) buffers_finished++;
            if(buffers[t].i == buffer_size && (buffers[t].i + buffers[t].curr_pos) < buffers[t].limit){
                // Read another piece
                buffers[t].curr_pos += buffers[t].i;
                fseek(hits_file, buffers[t].curr_pos * sizeof(hit), SEEK_SET);
                read = fread(buffers[i].h, sizeof(hit), buffer_size, hits_file);
                buffers[t].i = 0;
            }

            

            if(minimum_value > buffers[t].h[buffers[t].i].diag){ minimum_value = buffers[t].h[buffers[t].i].diag; remember_queue = t; }
            fpq->insert(&buffers[t].h[buffers[t].i]);
            buffers[t].i++;
            
        }
        // Re apply to queue with the minimum
        if(remember_queue != -1){
                fpq->insert(&buffers[remember_queue].h[buffers[remember_queue].i]); buffers[remember_queue].i++;
        }
        
        // Write to disk
        fpq->liberate_hits(minimum_value);

    }


    
    free(frame_positions);

    for(i=0; i<n_frames; i++){
        free(&buffers[i]);
    }
    free(buffers);
    fclose(hits_filtered);
    
    fprintf(stdout, "[INFO] Completed processing query splits\n");
    
    return 0;
}


void init_args(int argc, char ** av, FILE ** hits_file, FILE ** hits_file_headers, cl_uint * selected_device, FILE ** hits_filtered, ulong * fixed_power_of_two, cl_ulong * queue_size){
    
    int pNum = 0;
    
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           mergehits -f [file] -o [file]\n");
            fprintf(stdout, "OPTIONAL:\n");
            
            fprintf(stdout, "           -dev        [Integer: d>=0] Selects the device to be used\n");
            fprintf(stdout, "           -power      [Integer: k>=1] Use fixed size power of 2 for sorting\n");
            fprintf(stdout, "           -q          [Integer: q>=1] Fixed size for merge queue\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }

        if(strcmp(av[pNum], "-f") == 0){
            *hits_file = fopen(av[pNum+1], "rb");
            if(*hits_file==NULL){ fprintf(stderr, "Could not open input hits file\n"); exit(-1); }
            char temp[2048]; temp[0] = '\0';
            strcpy(temp, av[pNum+1]);
            strcat(temp, ".headers");
            *hits_file_headers = fopen(temp, "rt");
            if(*hits_file_headers==NULL){ fprintf(stderr, "Could not open input hits headers file\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-o") == 0){
            *hits_filtered = fopen(av[pNum+1], "wb");
            if(*hits_filtered==NULL){ fprintf(stderr, "Could not open output filtered hits file\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-dev") == 0){
            *selected_device = (cl_uint) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 0) { fprintf(stderr, "Device must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-power") == 0){
            *fixed_power_of_two = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) <= 0) { fprintf(stderr, "Power of 2 must be >0\n"); exit(-1); }
        }
        if(strcmp(av[pNum], "-q") == 0){
            *queue_size = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) <= 0) { fprintf(stderr, "Queue size must be >0\n"); exit(-1); }
        }

        pNum++;

    }   
    
    if(*hits_file==NULL || *hits_filtered==NULL){ fprintf(stderr, "You have to include an input and output hits file!\n"); exit(-1); }
}

char * get_dirname(char * path){
    int pos_last = 0, i = 0;
    while(path[i] != '\0'){
        if(path[i] == '/'){
            pos_last = i;
        }
        ++i;
    }
    char * dirname = (char *) malloc(BUFFER_SIZE * sizeof(char));
    if(dirname == NULL){ fprintf(stderr, "Could not allocate dirname char\n"); exit(-1); }

    memcpy(&dirname[0], &path[0], pos_last);
    dirname[pos_last] = '\0';

    return dirname;
}

char * get_basename(char * path){
    char * s = strrchr(path, '/');
    if (!s) return strdup(path); else return strdup(s + 1);
}