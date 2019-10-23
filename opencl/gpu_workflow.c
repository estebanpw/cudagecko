
// Standard utilities and common systems includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "structs.h"
#include "commonFunctions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define BUFFER_SIZE 2048
#define MAX_KERNEL_SIZE BUFFER_SIZE*100
#define CORES_PER_COMPUTE_UNIT 32
#define DIMENSION 1000

uint64_t filter_hits(hit_advanced * hits_in, ulong kmer_size, ulong n_hits_found);
void bitonic_sort_words(ulong N_at_once, ulong kmer_size, ulong words_per_work_item, cl_mem * cl_words, cl_kernel * kernel_words_sort, size_t global_item_size, size_t local_item_size, cl_context * context, cl_command_queue * command_queue);
void bitonic_sort_hits(ulong n_hits_found, ulong words_per_work_item, cl_mem * cl_hits_list, cl_kernel * kernel_sort_hits, size_t local_item_size, cl_context * context, cl_command_queue * command_queue);
char * get_dirname(char * path);
char * get_basename(char * path);
void init_args(int argc, char ** av, FILE ** query, FILE ** ref, FILE ** frags, cl_uint * selected_device, ulong * kmer_size, ulong * min_len, float * min_sim, long double * expected_value, unsigned char * strand, ulong * N_at_once);
char * read_seq(FILE * f, uint64_t * l, uint64_t ** seq_index_numbers, uint64_t * n_seqs);
void quicksort_frags(struct FragFile * array, int64_t x, int64_t y);
void filter_frags(struct FragFile * array, uint64_t n_frags, FILE * outfrags, uint64_t len_query, uint64_t len_ref);

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
    cl_ulong words_per_work_item = 1, kmer_size = 32, min_len = 20, N_at_once = 0, n_seqs_x = 0, n_seqs_y = 0;
    float min_sim = 30.0;
    size_t work_group_size[3], work_group_size_global;
    cl_int ret;
    char * path_kernels = get_dirname(argv[0]);
    unsigned char strand = 'f';
    long double min_expected_value = 1.0;
    uint64_t * seq_index_x, * seq_index_y;
    fprintf(stdout, "[INFO] Working on directory %s\n", path_kernels);

    FILE * query = NULL, * ref = NULL, * frags = NULL;

    init_args(argc, argv, &query, &ref, &frags, &selected_device, &kmer_size, &min_len, &min_sim, &min_expected_value, &strand, &N_at_once);

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
    
    // Create an OpenCL context
    cl_context context = clCreateContext(NULL, 1, &devices[selected_device], NULL, NULL, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not create context. Error: %d\n", ret); exit(-1); }

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, devices[selected_device], 0, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not create command queue. Error: %d\n", ret); exit(-1); }

    

    ////////////////////////////////////////////////////////////////////////////////
    // Split into N_at_once work items (power of 2 to ease up sorting)
    ////////////////////////////////////////////////////////////////////////////////

    if(N_at_once == 0){
        N_at_once = (ulong) pow(2, 24); // It is 16 MB
    }else{
        N_at_once = (ulong) pow(2, N_at_once); // User defined
    }
    
    // Find closest power of two to use
    fprintf(stdout, "[INFO] Maximum words to load at once: %"PRIu64"\n", (uint64_t) N_at_once);
    /*
    if(fixed_power_of_two == 0){
        N_at_once = (ulong) floor(log2((float) N_at_once));
        N_at_once = (ulong) pow(2, (float) N_at_once);
    }else{
        N_at_once = (ulong) pow(2, (float) fixed_power_of_two);
        if(N_at_once * sizeof(hitGPU) > (global_device_RAM[selected_device]-100*1024*1024)){ fprintf(stderr, "Exceeding GPU RAM limits\n"); exit(-1); }
    }
    */
    fprintf(stdout, "[INFO] Reading frames of %"PRIu64" words, i.e. %"PRIu64" MB.\n", (uint64_t) N_at_once, (uint64_t) ((N_at_once * sizeof(word_hash_GPU)) / (1024*1024)));

    ////////////////////////////////////////////////////////////////////////////////
    // Load sequences into GPU ram
    ////////////////////////////////////////////////////////////////////////////////

    fprintf(stdout, "[INFO] Loading input sequences.\n");
    // First load into host
    uint64_t len_query = 0, len_ref = 0;
    char * seq_x = read_seq(query, &len_query, &seq_index_x, &n_seqs_x);
    char * seq_y = read_seq(ref, &len_ref, &seq_index_y, &n_seqs_y);
    fwrite(&len_query, sizeof(uint64_t), 1, frags);
    fwrite(&len_ref, sizeof(uint64_t), 1, frags);
    if(seq_x == NULL || seq_y == NULL){ fprintf(stderr, "Could not allocate memory for sequences\n"); exit(-1); }

    // Load into GPU
    cl_mem cl_seq_x = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, len_query * sizeof(char), seq_x, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for seq_x in device. Error: %d\n", ret); exit(-1); }

    cl_mem cl_seq_y = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, len_ref * sizeof(char), seq_y, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for seq_y in device. Error: %d\n", ret); exit(-1); }

    // No longer needed, its in GPU
    free(seq_x);
    free(seq_y);

    ////////////////////////////////////////////////////////////////////////////////
    // Create words kernel
    ////////////////////////////////////////////////////////////////////////////////

    // Load kernel
    FILE * read_kernel; 
    char kernel_temp_path[BUFFER_SIZE];
    kernel_temp_path[0] = '\0';
    strcat(kernel_temp_path, path_kernels);
    
    strcat(kernel_temp_path, "/kernels/kernel_words_advanced.cl"); 
    read_kernel = fopen(kernel_temp_path, "r");
    if(!read_kernel){ fprintf(stderr, "Failed to load kernel words.\n"); exit(-1); }
    char * source_str = (char *) malloc(MAX_KERNEL_SIZE);
    if(source_str == NULL) { fprintf(stderr, "Could not allocate kernel words\n"); exit(-1); }
    size_t source_size = fread(source_str, 1, MAX_KERNEL_SIZE, read_kernel);
    fclose(read_kernel);

    // Create a program from the kernel source
    cl_program program_words = clCreateProgramWithSource(context, 1, (const char **) &source_str, (const size_t *) &source_size, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating program words: %d\n", ret); exit(-1); }


    // Build the program
    ret = clBuildProgram(program_words, 1, &devices[selected_device], NULL, NULL, NULL);
    if(ret != CL_SUCCESS){ 
        fprintf(stderr, "Error building program words: %d\n", ret); 
        if (ret == CL_BUILD_PROGRAM_FAILURE) {
            // Determine the size of the log
            size_t log_size;
            clGetProgramBuildInfo(program_words, devices[selected_device], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
            // Allocate memory for the log
            char *log = (char *) malloc(log_size);
            // Get the log
            clGetProgramBuildInfo(program_words, devices[selected_device], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
            // Print the log
            fprintf(stdout, "%s\n", log);
            free(log);
        }
        exit(-1); 
    }

    
    // Create the OpenCL kernel
    cl_kernel kernel_words_advanced = clCreateKernel(program_words, "kernel_words_advanced", &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating kernel words: %d\n", ret); exit(-1); }
    

    ////////////////////////////////////////////////////////////////////////////////
    // Create sort words kernel
    ////////////////////////////////////////////////////////////////////////////////

    
    // Load kernel
    kernel_temp_path[0] = '\0';
    strcat(kernel_temp_path, path_kernels);
    strcat(kernel_temp_path, "/kernels/kernel_words_sort.cl"); 
    read_kernel = fopen(kernel_temp_path, "r");
    if(!read_kernel){ fprintf(stderr, "Failed to load kernel sort words.\n"); exit(-1); }
    source_size = fread(source_str, 1, MAX_KERNEL_SIZE, read_kernel);
    fclose(read_kernel);

    // Create a program from the kernel source
    cl_program program_words_sort = clCreateProgramWithSource(context, 1, (const char **) &source_str, (const size_t *) &source_size, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating program sort words: %d\n", ret); exit(-1); }


    // Build the program
    ret = clBuildProgram(program_words_sort, 1, &devices[selected_device], NULL, NULL, NULL);
    if(ret != CL_SUCCESS){ 
        fprintf(stderr, "Error building program sort words: %d\n", ret); 
        if (ret == CL_BUILD_PROGRAM_FAILURE) {
            // Determine the size of the log
            size_t log_size;
            clGetProgramBuildInfo(program_words_sort, devices[selected_device], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
            // Allocate memory for the log
            char *log = (char *) malloc(log_size);
            // Get the log
            clGetProgramBuildInfo(program_words_sort, devices[selected_device], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
            // Print the log
            fprintf(stdout, "%s\n", log);
            free(log);
        }
        exit(-1); 
    }

    
    // Create the OpenCL kernel
    cl_kernel kernel_words_sort = clCreateKernel(program_words_sort, "kernel_words_sort", &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating kernel sort words: %d\n", ret); exit(-1); }

    ////////////////////////////////////////////////////////////////////////////////
    // Create sort hits
    ////////////////////////////////////////////////////////////////////////////////

    
    // Load kernel
    kernel_temp_path[0] = '\0';
    strcat(kernel_temp_path, path_kernels);
    strcat(kernel_temp_path, "/kernels/global_sorting_hits_advanced.cl"); 
    read_kernel = fopen(kernel_temp_path, "r");
    if(!read_kernel){ fprintf(stderr, "Failed to load kernel frags.\n"); exit(-1); }
    source_size = fread(source_str, 1, MAX_KERNEL_SIZE, read_kernel);
    fclose(read_kernel);

    // Create a program from the kernel source
    cl_program program_sort_hits = clCreateProgramWithSource(context, 1, (const char **) &source_str, (const size_t *) &source_size, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating program sort hits: %d\n", ret); exit(-1); }


    // Build the program
    ret = clBuildProgram(program_sort_hits, 1, &devices[selected_device], NULL, NULL, NULL);
    if(ret != CL_SUCCESS){ 
        fprintf(stderr, "Error building program sort hits: %d\n", ret); 
        if (ret == CL_BUILD_PROGRAM_FAILURE) {
            // Determine the size of the log
            size_t log_size;
            clGetProgramBuildInfo(program_sort_hits, devices[selected_device], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
            // Allocate memory for the log
            char *log = (char *) malloc(log_size);
            // Get the log
            clGetProgramBuildInfo(program_sort_hits, devices[selected_device], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
            // Print the log
            fprintf(stdout, "%s\n", log);
            free(log);
        }
        exit(-1); 
    }

    
    // Create the OpenCL kernel
    cl_kernel kernel_sort_hits = clCreateKernel(program_sort_hits, "kernel_sort", &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating kernel sort hits: %d\n", ret); exit(-1); }

    ////////////////////////////////////////////////////////////////////////////////
    // Create frags kernel
    ////////////////////////////////////////////////////////////////////////////////

    
    // Load kernel
    kernel_temp_path[0] = '\0';
    strcat(kernel_temp_path, path_kernels);
    strcat(kernel_temp_path, "/kernels/kernel_frags_advanced.cl"); 
    read_kernel = fopen(kernel_temp_path, "r");
    if(!read_kernel){ fprintf(stderr, "Failed to load kernel frags.\n"); exit(-1); }
    source_size = fread(source_str, 1, MAX_KERNEL_SIZE, read_kernel);
    fclose(read_kernel);

    // Create a program from the kernel source
    cl_program program_frags = clCreateProgramWithSource(context, 1, (const char **) &source_str, (const size_t *) &source_size, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating program frags: %d\n", ret); exit(-1); }


    // Build the program
    ret = clBuildProgram(program_frags, 1, &devices[selected_device], NULL, NULL, NULL);
    if(ret != CL_SUCCESS){ 
        fprintf(stderr, "Error building program frags: %d\n", ret); 
        if (ret == CL_BUILD_PROGRAM_FAILURE) {
            // Determine the size of the log
            size_t log_size;
            clGetProgramBuildInfo(program_frags, devices[selected_device], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
            // Allocate memory for the log
            char *log = (char *) malloc(log_size);
            // Get the log
            clGetProgramBuildInfo(program_frags, devices[selected_device], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
            // Print the log
            fprintf(stdout, "%s\n", log);
            free(log);
        }
        exit(-1); 
    }

    
    // Create the OpenCL kernel
    cl_kernel kernel_frags = clCreateKernel(program_frags, "kernel_frags", &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating kernel frags: %d\n", ret); exit(-1); }




    ////////////////////////////////////////////////////////////////////////////////
    // Array to cache sorted words
    ////////////////////////////////////////////////////////////////////////////////

    
    uint64_t n_arrays_y = len_ref / N_at_once + 1;
    uint64_t offset_computed = 0;
    word_hash_GPU * words_cache_y = (word_hash_GPU *) malloc(n_arrays_y * N_at_once * sizeof(word_hash_GPU));
    if(words_cache_y == NULL){ fprintf(stderr, "Error allocating host cache Y words\n"); exit(-1); }





    ////////////////////////////////////////////////////////////////////////////////
    // Create parameters and working sizes
    ////////////////////////////////////////////////////////////////////////////////

    // Set working size
    size_t local_item_size = 256; // Number of work items in a work group
    size_t global_item_size = N_at_once/words_per_work_item + (local_item_size - (N_at_once/words_per_work_item) % local_item_size);

    fprintf(stdout, "[INFO] Work group size: %"PRIu64", work-items: %"PRIu64", splits go on %"PRIu64"\n", local_item_size, global_item_size, N_at_once);


    param_words_advanced pwa_x;
    pwa_x.offset = 0;
    pwa_x.end = len_query;
    pwa_x.kmer_size = kmer_size;
    pwa_x.kmers_in_work_item = words_per_work_item;
    pwa_x.t_work_items = global_item_size;

    param_words_advanced pwa_y;
    pwa_y.offset = 0;
    pwa_y.end = len_ref;
    pwa_y.kmer_size = kmer_size;
    pwa_y.kmers_in_work_item = words_per_work_item;
    pwa_y.t_work_items = global_item_size;
    

    // Allocate host memory to retrieve words
    word_hash_GPU * words_x = (word_hash_GPU *) malloc(N_at_once * sizeof(word_hash_GPU));
    if(words_x == NULL){ fprintf(stderr, "Error allocating host words\n"); exit(-1); }

    // Initialize at 0xFF... because if sequence is smaller than N_at_once it yields bad words
    // The int value is converted to unsigned char
    for(i=0; i<N_at_once; i++){
        words_x[i].h = 0xFFFFFFFFFFFFFFFF;
        words_x[i].pos = 0xFFFFFFFFFFFFFFFF;
        words_cache_y[i].h = 0xFFFFFFFFFFFFFFFF;
        words_cache_y[i].pos = 0xFFFFFFFFFFFFFFFF;
    }
    for(i=N_at_once; i<n_arrays_y * N_at_once; i++){
        words_cache_y[i].h = 0xFFFFFFFFFFFFFFFF;
        words_cache_y[i].pos = 0xFFFFFFFFFFFFFFFF;
    }

    // Allocate words list
    cl_mem cl_words_list_x = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N_at_once * sizeof(word_hash_GPU), words_x, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for word list x in device. Error: %d\n", ret); exit(-1); }

    cl_mem cl_words_list_y = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N_at_once * sizeof(word_hash_GPU), words_x, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for word list y in device. Error: %d\n", ret); exit(-1); }

    // Allocate hits
    hit_advanced * hits = (hit_advanced *) calloc(N_at_once, sizeof(hit_advanced));
    if(hits == NULL){ fprintf(stderr, "Error allocating host hits\n"); exit(-1); }


    // Allocate space on host for writing the output
    reduced_frag * frag_write = (reduced_frag *) calloc(N_at_once, sizeof(reduced_frag));
    if(frag_write == NULL){ fprintf(stderr, "Could not allocate memory for host frags"); exit(-1); }

    // Allocate frags on GPU
    cl_mem cl_frags = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N_at_once * sizeof(reduced_frag), frag_write, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for frags in device. Error: %d\n", ret); exit(-1); }

    uint64_t all_frags_forever = 0;

    struct FragFile * final_frags_array = (struct FragFile *) malloc(N_at_once * sizeof(struct FragFile));
    if(final_frags_array == NULL) { fprintf(stderr, "Could not allocate storage for frags\n"); exit(-1); }
    uint64_t reallocs_frags = 1;

    ////////////////////////////////////////////////////////////////////////////////
    // Start computing
    ////////////////////////////////////////////////////////////////////////////////



    // Read the input query in chunks
    while(pwa_x.offset < len_query) {

        // Reload params
        cl_mem cl_params_x = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(param_words_advanced), &pwa_x, &ret);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for word list x in device. Error: %d\n", ret); exit(-1); }

        // Set kernel arguments
        ret = clSetKernelArg(kernel_words_advanced, 0, sizeof(cl_mem), (void *)&cl_words_list_x);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (wx1): %d\n", ret); exit(-1); }
        ret = clSetKernelArg(kernel_words_advanced, 1, sizeof(cl_mem), (void *)&cl_params_x);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (wx2): %d\n", ret); exit(-1); }
        ret = clSetKernelArg(kernel_words_advanced, 2, sizeof(cl_mem), (void *)&cl_seq_x);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (wx3): %d\n", ret); exit(-1); }

        // Execute the OpenCL kernel on the lists
        ret = clEnqueueNDRangeKernel(command_queue, kernel_words_advanced, 1, NULL, 
                &global_item_size, &local_item_size, 0, NULL, NULL);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Error enqueuing kernel (1): %d\n", ret); exit(-1); }

        // Wait for kernel to finish
        ret = clFlush(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
        ret = clFinish(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event (X-words): %d\n", ret); exit(-1); }

        ////////////////////////////////////////////////////////////////////////////////
        // Bitonic algorithm for words sorting the X sequence
        ////////////////////////////////////////////////////////////////////////////////
        bitonic_sort_words(N_at_once, kmer_size, words_per_work_item, &cl_words_list_x, &kernel_words_sort, global_item_size, local_item_size, &context, &command_queue);
        // Retrieve sorted X words
        ret = clEnqueueReadBuffer(command_queue, cl_words_list_x, CL_TRUE, 0, N_at_once*sizeof(word_hash_GPU), words_x, 0, NULL, NULL);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not read from buffer: %d\n", ret); exit(-1); }


        while(pwa_y.offset < len_ref) {

            if(offset_computed <= pwa_y.offset){
                // Reload params
                cl_mem cl_params_y = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(param_words_advanced), &pwa_y, &ret);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for word list x in device. Error: %d\n", ret); exit(-1); }

                // Set kernel arguments
                ret = clSetKernelArg(kernel_words_advanced, 0, sizeof(cl_mem), (void *)&cl_words_list_y);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (wy1): %d\n", ret); exit(-1); }
                ret = clSetKernelArg(kernel_words_advanced, 1, sizeof(cl_mem), (void *)&cl_params_y);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (wy2): %d\n", ret); exit(-1); }
                ret = clSetKernelArg(kernel_words_advanced, 2, sizeof(cl_mem), (void *)&cl_seq_y);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (wy3): %d\n", ret); exit(-1); }

                // Execute the OpenCL kernel on the lists
                ret = clEnqueueNDRangeKernel(command_queue, kernel_words_advanced, 1, NULL, 
                        &global_item_size, &local_item_size, 0, NULL, NULL);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Error enqueuing kernel (1): %d\n", ret); exit(-1); }

                // Wait for kernel to finish
                ret = clFlush(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
                ret = clFinish(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event(Y-words): %d\n", ret); exit(-1); }                       

                ret = clReleaseMemObject(cl_params_y); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (6)\n"); exit(-1); }
            

                ////////////////////////////////////////////////////////////////////////////////
                // Bitonic algorithm for words sorting the Y sequence
                ////////////////////////////////////////////////////////////////////////////////
                bitonic_sort_words(N_at_once, kmer_size, words_per_work_item, &cl_words_list_y, &kernel_words_sort, global_item_size, local_item_size, &context, &command_queue);
                ret = clEnqueueReadBuffer(command_queue, cl_words_list_y, CL_TRUE, 0, N_at_once*sizeof(word_hash_GPU), &words_cache_y[offset_computed], 0, NULL, NULL);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Could not read from buffer: %d\n", ret); exit(-1); }
                offset_computed += N_at_once;

            }


            ////////////////////////////////////////////////////////////////////////////////
            // Generate hits for split
            ////////////////////////////////////////////////////////////////////////////////



            uint64_t id_x = 0, id_y = 0, n_hits_found = 0;
            //printf("before gen hits row %"PRIu64" pwa %"PRIu64"\n", pwa_x.offset, pwa_y.offset);
            while(id_x < N_at_once || id_y < N_at_once) {
                // Compare
                //if(pwa_x.offset == 16777216 && pwa_y.offset == 67108864) fprintf(stdout, "indices: %"PRIu64", %"PRIu64", h:%"PRIu64", %"PRIu64" at %"PRIu64", %"PRIu64"\n", id_x, id_y, words_x[id_x].h, words_cache_y[id_y + pwa_y.offset].h, words_x[id_x].pos, words_cache_y[id_y + pwa_y.offset].pos);
                if(n_hits_found > N_at_once) { fprintf(stderr, "Too many hits at once: %"PRIu64", current max: %"PRIu64"\n", n_hits_found, N_at_once); exit(-1); }
                if(words_x[id_x].pos == 0xFFFFFFFFFFFFFFFF) break;
                if(words_cache_y[id_y + pwa_y.offset].pos == 0xFFFFFFFFFFFFFFFF) break;
                if(id_x == N_at_once && words_cache_y[id_y + pwa_y.offset].h == 0xFFFFFFFFFFFFFFFF) break;
                if(id_y == N_at_once && words_x[id_x].h == 0xFFFFFFFFFFFFFFFF) break;
                if(id_x == N_at_once && id_y == N_at_once) break;
                
                
                if(words_x[id_x].h == words_cache_y[id_y + pwa_y.offset].h) {
                    // This is a hit
                    //printf("Made hit: ");
                    hits[n_hits_found].pos_x = words_x[id_x].pos;
                    hits[n_hits_found].pos_y = words_cache_y[id_y + pwa_y.offset].pos;

                    ++n_hits_found;
                    if(n_hits_found == N_at_once){ fprintf(stderr, "Reached maximum limit of hits\n"); exit(-1); }
                    if(id_x < N_at_once){ ++id_x; }
                    if(id_y < N_at_once){ ++id_y; }
                    //printf("\n");
                }
                else if(words_x[id_x].h < words_cache_y[id_y + pwa_y.offset].h){  if(id_x < N_at_once) ++id_x;  }
                else { if(id_y < N_at_once) ++id_y; }
                
            }
            //printf("after gen hits\n");

            
            
            

            

            ////////////////////////////////////////////////////////////////////////////////
            // Add some hits sorting & filtering
            ////////////////////////////////////////////////////////////////////////////////

            cl_mem cl_hits_list = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, n_hits_found * sizeof(hit_advanced), hits, &ret);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for word list x in device. Error: %d\n", ret); exit(-1); }

            bitonic_sort_hits(n_hits_found, words_per_work_item, &cl_hits_list, &kernel_sort_hits, local_item_size, &context, &command_queue);

            ret = clEnqueueReadBuffer(command_queue, cl_hits_list, CL_TRUE, 0, n_hits_found * sizeof(hit_advanced), &hits[0], 0, NULL, NULL);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Could not read from buffer: %d\n", ret); exit(-1); }

            ret = clReleaseMemObject(cl_hits_list); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (hits)\n"); exit(-1); }
            if(ret != CL_SUCCESS){ fprintf(stderr, "Could not release hits: %d\n", ret); exit(-1); }

            // Filter hits
            uint64_t hits_kept = filter_hits(hits, kmer_size, n_hits_found);

            //fprintf(stdout, "[INFO] Finished generating hits = %"PRIu64" for row %"PRIu64", cell %"PRIu64". After filtering: %"PRIu64"\n", n_hits_found, pwa_x.offset, pwa_y.offset, hits_kept);


            ////////////////////////////////////////////////////////////////////////////////
            // Go for the frags, you mad?
            ////////////////////////////////////////////////////////////////////////////////

            if(n_hits_found > 0){

                //fprintf(stdout, "[INFO] Generating frags\n");

                // Copy hits to device
                cl_mem cl_hits = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, hits_kept * sizeof(hit_advanced), hits, &ret);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for hits in device. Error: %d\n", ret); exit(-1); }

                // Params for the frags computation
                parameter_frags pfrags;
                pfrags.kmer_size = kmer_size;
                pfrags.size_x = len_query;
                pfrags.size_y = len_ref;
                pfrags.t_hits = hits_kept;

                // Adjust work items
                size_t global_item_size_hits = hits_kept + (local_item_size - (hits_kept) % local_item_size);

                // Allocate params
                cl_mem params_frags = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(parameter_frags), &pfrags, &ret);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for ranndom numbers in device. Error: %d\n", ret); exit(-1); }



                // Set the arguments of the kernel
                ret = clSetKernelArg(kernel_frags, 0, sizeof(cl_mem), (void *)&cl_seq_x);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (1): %d\n", ret); exit(-1); }

                ret = clSetKernelArg(kernel_frags, 1, sizeof(cl_mem), (void *)&cl_seq_y);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (2): %d\n", ret); exit(-1); }

                ret = clSetKernelArg(kernel_frags, 2, sizeof(cl_mem), (void *)&cl_hits);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (3): %d\n", ret); exit(-1); }

                ret = clSetKernelArg(kernel_frags, 3, sizeof(cl_mem), (void *)&cl_frags);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (4): %d\n", ret); exit(-1); }

                ret = clSetKernelArg(kernel_frags, 4, sizeof(cl_mem), (void *)&params_frags);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (5): %d\n", ret); exit(-1); }


                ret = clEnqueueNDRangeKernel(command_queue, kernel_frags, 1, NULL, 
                    &global_item_size_hits, &local_item_size, 0, NULL, NULL);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Error enqueuing kernel frags: %d\n", ret); exit(-1); }

                // Wait for kernel to finish
                ret = clFlush(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
                ret = clFinish(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event(frags): %d\n", ret); exit(-1); }


                // Release params and hits
                ret = clReleaseMemObject(cl_hits); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (hits)\n"); exit(-1); }
                if(ret != CL_SUCCESS){ fprintf(stderr, "Could not release hits: %d\n", ret); exit(-1); }
                ret = clReleaseMemObject(params_frags); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (params frags)\n"); exit(-1); }
                if(ret != CL_SUCCESS){ fprintf(stderr, "Could not release hits: %d\n", ret); exit(-1); }

                // Read frags and treat them to a nice dinner
                ret = clEnqueueReadBuffer(command_queue, cl_frags, CL_TRUE, 0, hits_kept*sizeof(reduced_frag), frag_write, 0, NULL, NULL);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Could not read from buffer: %d\n", ret); exit(-1); }

                struct FragFile frag;
                uint64_t total_frags = 0;
                for(i=0; i<hits_kept; i++){

                    frag.length = frag_write[i].length;
                    frag.ident = frag_write[i].identities;
                    frag.score = frag_write[i].identities * POINT + (frag_write[i].length - frag_write[i].identities) * (-POINT);
                    frag.similarity = 100*(float) frag.score / (float) (frag_write[i].length * POINT);
                    frag.evalue = 0.333 * len_query * len_ref * expl(-0.275 * frag.score); //0.333 is karlin and 0.275 is lambda

                    if(frag.similarity > min_sim && frag.length > min_len && frag.evalue < min_expected_value){

                        frag.diag = (int64_t) frag_write[i].pos_x - (int64_t) frag_write[i].pos_y;
                        frag.xStart = frag_write[i].pos_x - seq_index_x[frag_write[i].pos_x];
                        frag.xEnd = (frag_write[i].pos_x + frag_write[i].length) - seq_index_x[frag_write[i].pos_x];

                        

                        if(strand == 'r'){
                            frag.yStart = frag_write[i].pos_y - seq_index_y[frag_write[i].pos_y];
                            frag.yEnd = (frag_write[i].pos_y + frag_write[i].length) - seq_index_y[frag_write[i].pos_y];

                        }else{
                            frag.yStart = frag_write[i].pos_y - seq_index_y[frag_write[i].pos_y];
                            frag.yEnd = (frag_write[i].pos_y + frag_write[i].length) - seq_index_y[frag_write[i].pos_y];
                        }
                        
                        
                        frag.seqX = seq_index_x[frag_write[i].pos_x];


                        if(strand == 'f') frag.seqY = seq_index_y[frag_write[i].pos_y]; else frag.seqY = n_seqs_y - seq_index_y[frag_write[i].pos_y];


                        frag.block = 0;
                        frag.strand = strand;

                        //fprintf(stdout, "(1) Hit and (2) Frag\n");
                        //print_hit(&hgpu[i]);
                        //print_fragment(&frag);
                        fwrite(&frag, sizeof(struct FragFile), 1, frags);
                        final_frags_array[total_frags] = frag;
                        ++total_frags;
                        if(total_frags == reallocs_frags*N_at_once){
                            ++reallocs_frags;
                            final_frags_array = (struct FragFile *) realloc(final_frags_array, reallocs_frags * N_at_once * sizeof(struct FragFile));
                            if(final_frags_array == NULL) { fprintf(stderr, "Could not reallocate fragments\n"); exit(-1); }
                        }
                        //getchar();
                    }

                }
                all_frags_forever += total_frags;
                //fprintf(stdout, "[INFO] Generated %"PRIu64" frags\n", total_frags);
            }

            









            // Advance pointers
            pwa_y.offset += N_at_once;
            
            
        }
        
        fprintf(stdout, "[INFO] Finished row %"PRIu64". Current frags = %"PRIu64"\n", pwa_x.offset, all_frags_forever);
        // Advance pointers
        pwa_x.offset += N_at_once;
        pwa_y.offset = 0;

        ret = clReleaseMemObject(cl_params_x); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (6)\n"); exit(-1); }
    }

    ret = clReleaseMemObject(cl_words_list_x); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (7)\n"); exit(-1); }
    ret = clReleaseMemObject(cl_words_list_y); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (8)\n"); exit(-1); }
    ret = clReleaseMemObject(cl_seq_x); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (9)\n"); exit(-1); }
    ret = clReleaseMemObject(cl_seq_y); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (10)\n"); exit(-1); }

    fprintf(stdout, "[INFO] Completed processing splits. T frags = %"PRIu64"\n", all_frags_forever);

    //fprintf(stdout, "[INFO] Sorting frags\n");
    //quicksort_frags(final_frags_array, 0, all_frags_forever-1);
    
    /*
    for(i=0; i<all_frags_forever; i++){
        fprintf(stdout, "d=%"PRId64" sx=%"PRIu64", sy=%"PRIu64", x=%"PRIu64", y=%"PRIu64"\n", final_frags_array[i].diag, final_frags_array[i].seqX, final_frags_array[i].seqY, final_frags_array[i].xStart, final_frags_array[i].yStart);
    }
    printf("AFTER: ................\n");
    */
    //filter_frags(final_frags_array, all_frags_forever, frags, len_query, len_ref);

    free(words_cache_y);
    free(source_str);
    free(words_x);
    free(hits);
    free(frag_write);
    free(final_frags_array);
    
    return 0;
}

uint64_t filter_hits(hit_advanced * hits_in, ulong kmer_size, ulong n_hits_found){
    int64_t diagonal;
    uint64_t lastPosition = 0, t_filtered = 0, t = 1, t_kept = 0;
    diagonal = (int64_t)hits_in[0].pos_x - (int64_t)hits_in[0].pos_y;
    if(n_hits_found == 0) return 0;
    while (t < (n_hits_found-1)) {

        if(diagonal == ((int64_t)hits_in[t+1].pos_x - (int64_t)hits_in[t+1].pos_y) && hits_in[t].pos_x < (lastPosition+2*kmer_size)){
            // Delete this hit
            diagonal = 9223372036854775807;
            lastPosition = 0;
            ++t_filtered;
        }else{
            lastPosition = hits_in[t+1].pos_x + (2 * kmer_size - 1);
            diagonal = (int64_t)hits_in[t+1].pos_x - (int64_t)hits_in[t+1].pos_y;
            memmove(&hits_in[t_kept], &hits_in[t+1], sizeof(hit_advanced));
            ++t_kept;
        }
        ++t;
    }
    return t_kept;
}

void bitonic_sort_words(ulong N_at_once, ulong kmer_size, ulong words_per_work_item, cl_mem * cl_words, cl_kernel * kernel_words_sort, size_t global_item_size, size_t local_item_size, cl_context * context, cl_command_queue * command_queue){

    long logn = (long)log2((float)(N_at_once));


    ////////////////////////////////////////////////////////////////////////////////
    // Bitonic algorithm control loop
    ////////////////////////////////////////////////////////////////////////////////

    int ret;
    parameter_sort psort;
    psort.N = N_at_once;
    psort.kmer_size = kmer_size;
    psort.comparators_per_wi = words_per_work_item;

    long step, stage;
    
    for(step=1; step<=logn; step++){
        for(stage=step; stage>0; stage--){

            psort.step = (ulong) step;
            psort.stage = (ulong) stage;

            // Here goes kernel execution with new params

            cl_mem params = clCreateBuffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(parameter_sort), &psort, &ret);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for param sort in device. Error: %d\n", ret); exit(-1); }

            // Set kernel arguments
            ret = clSetKernelArg(*kernel_words_sort, 0, sizeof(cl_mem), (void *)cl_words);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (sw1): %d\n", ret); exit(-1); }
            ret = clSetKernelArg(*kernel_words_sort, 1, sizeof(cl_mem), (void *)&params);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (sw2): %d\n", ret); exit(-1); }

            // Execute the OpenCL kernel on the lists
            ret = clEnqueueNDRangeKernel(*command_queue, *kernel_words_sort, 1, NULL, 
                    &global_item_size, &local_item_size, 0, NULL, NULL);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Error enqueuing kernel (1): %d\n", ret); exit(-1); }

            // Wait for kernel to finish
            ret = clFlush(*command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
            ret = clFinish(*command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event(sort-words): %d\n", ret); exit(-1); }

            ret = clReleaseMemObject(params); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (6)\n"); exit(-1); }


        }

    }
}

void bitonic_sort_hits(ulong n_hits_found, ulong words_per_work_item, cl_mem * cl_hits_list, cl_kernel * kernel_sort_hits, size_t local_item_size, cl_context * context, cl_command_queue * command_queue){
    long logn = (long)log2((float)(n_hits_found));


    ////////////////////////////////////////////////////////////////////////////////
    // Bitonic algorithm control loop
    ////////////////////////////////////////////////////////////////////////////////

    int ret;
    parameter_sort psort;
    psort.N = n_hits_found;
    psort.comparators_per_wi = words_per_work_item;

    size_t global_item_size = n_hits_found/words_per_work_item + (local_item_size - (n_hits_found/words_per_work_item) % local_item_size);

    long step, stage;
    
    for(step=1; step<=logn; step++){
        for(stage=step; stage>0; stage--){

            psort.step = (ulong) step;
            psort.stage = (ulong) stage;

            // Here goes kernel execution with new params

            cl_mem params = clCreateBuffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(parameter_sort), &psort, &ret);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for param sort in device. Error: %d\n", ret); exit(-1); }

            // Set kernel arguments
            ret = clSetKernelArg(*kernel_sort_hits, 0, sizeof(cl_mem), (void *)cl_hits_list);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (sw1): %d\n", ret); exit(-1); }
            ret = clSetKernelArg(*kernel_sort_hits, 1, sizeof(cl_mem), (void *)&params);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (sw2): %d\n", ret); exit(-1); }

            // Execute the OpenCL kernel on the lists
            ret = clEnqueueNDRangeKernel(*command_queue, *kernel_sort_hits, 1, NULL, 
                    &global_item_size, &local_item_size, 0, NULL, NULL);
            if(ret != CL_SUCCESS){ fprintf(stderr, "Error enqueuing kernel (1): %d\n", ret); exit(-1); }

            // Wait for kernel to finish
            ret = clFlush(*command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
            ret = clFinish(*command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event(sort-hits): %d\n", ret); exit(-1); }

            ret = clReleaseMemObject(params); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (6)\n"); exit(-1); }


        }

    }
}

void init_args(int argc, char ** av, FILE ** query, FILE ** ref, FILE ** frags, cl_uint * selected_device, ulong * kmer_size, ulong * min_len, float * min_sim, long double * expected_value, unsigned char * strand, ulong * N_at_once){
    
    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           gpu_workflow -q [file] -r [file] -frags [file]\n");
            fprintf(stdout, "OPTIONAL:\n");
            
            fprintf(stdout, "           -dev        [Integer: d>=0] Selects the device to be used\n");
            fprintf(stdout, "           -kmer       [Integer: k>=1] Size of K-mer to be used\n");
            fprintf(stdout, "           -power      [Integer: p>0 ] Power of two to process splits\n");
            fprintf(stdout, "           -l          [Integer: k>=1] Minimum length\n");
            fprintf(stdout, "           -s          [Float  : 0<k<100] Minimum similarity\n");
            fprintf(stdout, "           -e          [Float  : 0<=e] Minimum attained e-value\n");
            fprintf(stdout, "           --r         It is a reverse comparison\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }

        if(strcmp(av[pNum], "-frags") == 0){
            *frags = fopen(av[pNum+1], "wb");
        }

        if(strcmp(av[pNum], "-q") == 0){
            *query = fopen(av[pNum+1], "rt");
            if(*query==NULL){ fprintf(stderr, "Could not open input query sequence file\n"); exit(-1); }
        }
        if(strcmp(av[pNum], "-r") == 0){
            *ref = fopen(av[pNum+1], "rt");
            if(*ref==NULL){ fprintf(stderr, "Could not open input ref sequence file\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-dev") == 0){
            *selected_device = (cl_uint) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 0) { fprintf(stderr, "Device must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-kmer") == 0){
            *kmer_size = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) <= 0) { fprintf(stderr, "Kmer size must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-power") == 0){
            *N_at_once = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) <= 0) { fprintf(stderr, "Power must be >0\n"); exit(-1); }
        }
        
        if(strcmp(av[pNum], "--r") == 0){
            *strand = 'r';
        }
        if(strcmp(av[pNum], "-l") == 0){
            *min_len = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) <= 0) { fprintf(stderr, "Min len size must be >0\n"); exit(-1); }
        }
        if(strcmp(av[pNum], "-s") == 0){
            *min_sim = atof(av[pNum+1])/100;
            if(*min_sim <= 0 || *min_sim >= 1) { fprintf(stderr, "Min sim must be >0 and <100\n"); exit(-1); }
        }
        if(strcmp(av[pNum], "-e") == 0){
            *expected_value = (long double) atof(av[pNum+1]);
            if(atof(av[pNum+1]) < 0) { fprintf(stderr, "Min expected value must be >=0\n"); exit(-1); }
        }




        pNum++;

    }   
    

    if(*query==NULL || *ref==NULL || *frags==NULL){ fprintf(stderr, "Include query and reference sequence and output frags\n"); exit(-1); }
    
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


char * read_seq(FILE * f, uint64_t * l, uint64_t ** seq_index_numbers, uint64_t * n_seqs) {
    char c;
    uint64_t lon = 0, k = 0;
    uint64_t SIZE = 0;
    uint64_t i = 0, r = 0, curr_idx = 0;
    char * seq = NULL;
    char * buffer = NULL;
    uint64_t * seq_index = NULL;


    //Memory

    fseek(f, 0, SEEK_END);
    SIZE = ftell(f);
    fseek(f, 0, SEEK_SET);

    if ((seq = (char *) calloc(SIZE, sizeof(char))) == NULL) {
        fprintf(stderr, "Could not allocate sequence\n"); exit(-1);
    }
    if ((buffer = (char *) calloc(READBUF, sizeof(char))) == NULL) {
        fprintf(stderr, "Could not allocate buffer\n"); exit(-1);
    }
    if ((seq_index = (uint64_t *) calloc(SIZE, sizeof(uint64_t))) == NULL) {
        fprintf(stderr, "Could not allocate sequence index\n"); exit(-1);
    }

    i = READBUF + 1;

    
    while ((c = buffered_fgetc(buffer, &i, &r, f)) != '>' && (!feof(f) || (feof(f) &&  i < r ) )); //start seq

    if (feof(f) && i >= r)
        return 0;
    while ((c = buffered_fgetc(buffer, &i, &r, f)) == ' ');

    while (k < 50 && c != '\n' && c != ' ') {
        if (feof(f) && i >= r)
            return 0;

        //sX->ident[k++] = c;
        c = buffered_fgetc(buffer, &i, &r, f);
    }

    //sX->ident[k] = 0; //end of data.
    while (c != '\n')
        c = buffered_fgetc(buffer, &i, &r, f);
    c = buffered_fgetc(buffer, &i, &r, f);

    //start list with sX2
    while (!feof(f) || (feof(f) &&  i < r ) ) {
        c = toupper(c);
        if (c == '>') {            
            seq_index[lon] = curr_idx;
            seq[lon++] = '*';
            ++curr_idx;
            while (c != '\n') {
                if ((feof(f) &&  i >= r ))
                    return 0;
                c = buffered_fgetc(buffer, &i, &r, f);
            }
            //break;
        }
        if (isupper(c) && c != '\n'){
            seq_index[lon] = curr_idx;
            seq[lon++] = c;
        }
        if (c == '*') {
            seq_index[lon] = curr_idx;
            seq[lon++] = c;
        }
        c = buffered_fgetc(buffer, &i, &r, f);
    }

    free(buffer);
    seq[lon] = 0x00;

    *l = lon;

    *seq_index_numbers = seq_index;

    *n_seqs = curr_idx;

    return seq;
}

void quicksort_frags(struct FragFile * array, int64_t x, int64_t y) {

    struct FragFile pivot, aux;
    int64_t x1, y1;

    pivot = array[(x+y)/2];

    x1 = x;
    y1 = y;

    do { 
        while (((int64_t)pivot.xStart-(int64_t)pivot.yStart) > ((int64_t)array[x1].xStart-(int64_t)array[x1].yStart)) x1++;
        while (((int64_t)pivot.xStart-(int64_t)pivot.yStart) < ((int64_t)array[y1].xStart-(int64_t)array[y1].yStart)) y1--;
        
        if (x1 < y1) { 
            aux = array[x1];
            array[x1] = array[y1]; 
            array[y1] = aux;
            x1++;
            y1--;
        }else if (x1 == y1) x1++;
    } while (x1 <= y1);

    if (x<y1) quicksort_frags(array, x, y1);
    if (x1<y) quicksort_frags(array, x1, y);
}

void filter_frags(struct FragFile * array, uint64_t n_frags, FILE * outfrags, uint64_t len_query, uint64_t len_ref){
    uint64_t i, j = 0;
    if(n_frags == 0) return;
    for(i=1; i<n_frags; i++){

        // Same diagonal and next frag start is included in previous frag
        if( (((int64_t)array[j].xStart - (int64_t)array[j].yStart) == ((int64_t)array[i].xStart - (int64_t)array[i].yStart)) 
        && array[i].xStart >= array[j].xStart && array[i].xStart <= array[j].xEnd){
                if(array[i].xEnd > array[j].xEnd){
                    array[j].ident += (array[i].length * array[i].ident) / (array[j].xEnd - array[j].xStart);
                    array[j].xEnd = array[i].xEnd;
                    array[j].yEnd = array[i].yEnd;
                    array[j].length = array[j].xEnd - array[j].xStart;

                    array[j].score = array[j].ident * POINT + (array[j].length - array[j].ident) * (-POINT);
                    array[j].similarity = 100*(float) array[i-1].score / (float) (array[j].length * POINT);
                    array[j].evalue = 0.333 * len_query * len_ref * expl(-0.275 * array[j].score); //0.333 is karlin and 0.275 is lambda
                } // Otherwise it is just fully included
                
                continue;
        }

        // Write the frag
        fwrite(&array[j], sizeof(struct FragFile), 1, outfrags);
        //fprintf(stdout, "d=%"PRId64" sx=%"PRIu64", sy=%"PRIu64", x=%"PRIu64", y=%"PRIu64"\n", array[j].diag, array[j].seqX, array[j].seqY, array[j].xStart, array[j].yStart);
        j = i;

    }
}