
// Standard utilities and common systems includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "structs.h"


#define BUFFER_SIZE 2048
#define MAX_KERNEL_SIZE BUFFER_SIZE*100
#define CORES_PER_COMPUTE_UNIT 32
#define DIMENSION 1000

char * get_dirname(char * path);
char * get_basename(char * path);
void init_args(int argc, char ** av, FILE ** query, cl_uint * selected_device, ulong * words_per_work_item, ulong * kmer_size, char * out);

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
    cl_ulong device_RAM, global_device_RAM;
    cl_ulong words_per_work_item, kmer_size;
    size_t work_group_size[3], work_group_size_global;
    cl_int ret;
    char outname[1024]; outname[0] = '\0';
    char * path_kernels = get_dirname(argv[0]);
    fprintf(stdout, "[INFO] Working on directory %s\n", path_kernels);

    FILE * query = NULL, * out = NULL;

    init_args(argc, argv, &query, &selected_device, &words_per_work_item, &kmer_size, outname);

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
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &global_device_RAM, NULL))){ fprintf(stderr, "Failed to get device global memory %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tGlobal mem   : %"PRIu64" (%"PRIu64" MB)\n", (uint64_t) global_device_RAM, (uint64_t) global_device_RAM / (1024*1024));
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &device_RAM, NULL))){ fprintf(stderr, "Failed to get device local memory %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tLocal mem    : %"PRIu64" (%"PRIu64" KB)\n", (uint64_t) device_RAM, (uint64_t) device_RAM / (1024));
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &compute_units, NULL))){ fprintf(stderr, "Failed to get device local memory %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tCompute units: %"PRIu64"\n", (uint64_t) compute_units);
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size_global, NULL))){ fprintf(stderr, "Failed to get device global work items size %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tMax work group size: %zu\n", work_group_size_global);
        if(CL_SUCCESS != (ret = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, 3*sizeof(size_t), &work_group_size, NULL))){ fprintf(stderr, "Failed to get device work items size %d\n", i); exit(-1); }
        fprintf(stdout, "\t\tWork size items: (%zu, %zu, %zu)\n", work_group_size[0], work_group_size[1], work_group_size[2]);
    }
    selected_device = 3;
    fprintf(stdout, "[INFO] Using device %d\n", selected_device);
    
    // Create an OpenCL context
    cl_context context = clCreateContext(NULL, 1, &devices[selected_device], NULL, NULL, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not create context. Error: %d\n", ret); exit(-1); }

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, devices[selected_device], 0, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not create command queue. Error: %d\n", ret); exit(-1); }

    ////////////////////////////////////////////////////////////////////////////////
    // Calculate how many words we can load at once
    ////////////////////////////////////////////////////////////////////////////////

    ulong N_at_once = (global_device_RAM-100*1024*1024) / sizeof(wordGPU); // Minus 100 MB
    
    // Find closest power of two to use
    fprintf(stdout, "[INFO] Maximum words to load at once: %"PRIu64"\n", (uint64_t) N_at_once);
    N_at_once = (ulong) floor(log2((float) N_at_once));
    N_at_once = (ulong) pow(2, (float) N_at_once);
    fprintf(stdout, "[INFO] Reading frames of %"PRIu64" words, i.e. %"PRIu64" MB.\n", (uint64_t) N_at_once, (uint64_t) ((N_at_once * sizeof(wordGPU)) / (1024*1024)));



    // Load kernel
    FILE * read_kernel; 
    char kernel_temp_path[BUFFER_SIZE];
    kernel_temp_path[0] = '\0';
    strcat(kernel_temp_path, path_kernels);
    
    strcat(kernel_temp_path, "/kernels/global_sorting.cl"); 
    read_kernel = fopen(kernel_temp_path, "r");
    if(!read_kernel){ fprintf(stderr, "Failed to load kernel (1).\n"); exit(-1); }
    char * source_str = (char *) malloc(MAX_KERNEL_SIZE);
    if(source_str == NULL) { fprintf(stderr, "Could not allocate kernel\n"); exit(-1); }
    size_t source_size = fread(source_str, 1, MAX_KERNEL_SIZE, read_kernel);
    fclose(read_kernel);

    // Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, (const char **) &source_str, (const size_t *) &source_size, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating program (1): %d\n", ret); exit(-1); }


    // Build the program
    ret = clBuildProgram(program, 1, &devices[selected_device], NULL, NULL, NULL);
    if(ret != CL_SUCCESS){ 
        fprintf(stderr, "Error building program (1): %d\n", ret); 
        if (ret == CL_BUILD_PROGRAM_FAILURE) {
            // Determine the size of the log
            size_t log_size;
            clGetProgramBuildInfo(program, devices[selected_device], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
            // Allocate memory for the log
            char *log = (char *) malloc(log_size);
            // Get the log
            clGetProgramBuildInfo(program, devices[selected_device], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
            // Print the log
            fprintf(stdout, "%s\n", log);
        }
        exit(-1); 
    }

    
    // Create the OpenCL kernel
    cl_kernel kernel = clCreateKernel(program, "kernel_sort", &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating kernel (1): %d\n", ret); exit(-1); }
    
    // Set working size
    size_t local_item_size = 256; // Number of work items in a work group
    size_t global_item_size = N_at_once/words_per_work_item + (local_item_size - (N_at_once/words_per_work_item) % local_item_size);

    fprintf(stdout, "[INFO] Work group size: %"PRIu64", work-items: %"PRIu64"\n", local_item_size, global_item_size);

    // Parameters for the kernels
    parameter_sort p;
    p.N = N_at_once;
    p.kmer_size = kmer_size;
    p.comparators_per_wi = words_per_work_item;

    // Allocate space on host for reading the words
    wordGPU * wgpu = (wordGPU *) malloc(N_at_once * sizeof(wordGPU));
    if(wgpu == NULL){ fprintf(stderr, "Could not allocate memory for host words"); exit(-1); }

    // Allocate space on host for writing the sorted words
    wordGPU * wgpu_sorted = (wordGPU *) malloc(N_at_once * sizeof(wordGPU));
    if(wgpu_sorted == NULL){ fprintf(stderr, "Could not allocate memory for host words"); exit(-1); }

    // This is the actual control loop that launches a kernel for each pass of the bitonic algorithm

    long logn = (ulong)log2((float)(N_at_once));
    ulong read_items;
    long sorted_frames = 0;

    while(!feof(query)){

        read_items = fread(wgpu, sizeof(wordGPU), N_at_once, query);
        if(read_items == 0){ fprintf(stderr, "Read 0 items\n"); exit(-1); }
        if(read_items < N_at_once){
            // In this case we have to make it a power of two (convert to bad words)
            ulong idx;
            for(idx=read_items; idx<N_at_once; idx++) wgpu[idx].pos = 0xFFFFFFFFFFFFFFFF;
        }

        // Allocate words
        cl_mem cl_words = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N_at_once * sizeof(wordGPU), wgpu, &ret);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for ranndom numbers in device. Error: %d\n", ret); exit(-1); }

        ////////////////////////////////////////////////////////////////////////////////
        // Bitonic algorithm control loop
        ////////////////////////////////////////////////////////////////////////////////


        long step, stage;
        for(step=1; step<=logn; step++){
            for(stage=step; stage>0; stage--){

                p.step = (ulong) step;
                p.stage = (ulong) stage;

                // Here goes kernel execution with new params

                cl_mem params = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(parameter_sort), &p, &ret);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for ranndom numbers in device. Error: %d\n", ret); exit(-1); }

                // Set the arguments of the kernel
                ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&cl_words);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (1): %d\n", ret); exit(-1); }

                ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&params);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (2): %d\n", ret); exit(-1); }

                ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, 
                    &global_item_size, &local_item_size, 0, NULL, NULL);
                if(ret != CL_SUCCESS){ fprintf(stderr, "Error enqueuing kernel (1): %d\n", ret); exit(-1); }

                // Wait for kernel to finish
                ret = clFlush(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
                ret = clFinish(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event: %d\n", ret); exit(-1); }

                ret = clReleaseMemObject(params); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (6)\n"); exit(-1); }

                fprintf(stdout, "[INFO] Enqueued kernel for %"PRIu64", %"PRIu64" on frame %"PRIu64"\n", (uint64_t) step, (uint64_t) stage, (uint64_t) sorted_frames);

            }
        }

        ret = clEnqueueReadBuffer(command_queue, cl_words, CL_TRUE, 0, N_at_once*sizeof(wordGPU), wgpu_sorted, 0, NULL, NULL);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not read from buffer: %d\n", ret); exit(-1); }

        char full_name[1024];
        full_name[0] = '\0';
        sprintf(full_name, "%s-%"PRIu64".sorted", outname, (uint64_t) sorted_frames);
        sorted_frames++;

        // Writing frames

        out = fopen(full_name, "wb");
        if(out == NULL){ fprintf(stderr, "Could not open output file %s\n", full_name); exit(-1); }
        fwrite(wgpu_sorted, sizeof(wordGPU), N_at_once, out);
        fclose(out);

        ret = clReleaseMemObject(cl_words); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (Words)\n"); exit(-1); }

    }

    fprintf(stdout, "[INFO] Completed execution of all kernels\n");

    /*
    ret = clEnqueueReadBuffer(command_queue, cl_randomnumbers, CL_TRUE, 0, N*sizeof(int), sortednumbers, 0, NULL, NULL);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not read from buffer: %d\n", ret); exit(-1); }
    */
    


    free(wgpu);
    free(wgpu_sorted);
    
    ret = clReleaseKernel(kernel); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (3)\n"); exit(-1); }
    ret = clReleaseProgram(program); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (4)\n"); exit(-1); }

    fprintf(stdout, "[INFO] Completed processing query splits\n");
    
    return 0;
}


void init_args(int argc, char ** av, FILE ** query, cl_uint * selected_device, ulong * words_per_work_item, ulong * kmer_size, char * out){
    
    int pNum = 0;
    char * p1;
    char outname[2048];
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           words -wfile [file]\n");
            fprintf(stdout, "OPTIONAL:\n");
            
            fprintf(stdout, "           -dev        [Integer: d>=0] Selects the device to be used\n");
            fprintf(stdout, "           -kmer       [Integer: k>=1] Size of K-mer to be used\n");
            fprintf(stdout, "           -wpi        [Integer: k>=1] Number of words per work item to be used\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }

        if(strcmp(av[pNum], "-wfile") == 0){
            *query = fopen(av[pNum+1], "rb");
            if(*query==NULL){ fprintf(stderr, "Could not open input words file\n"); exit(-1); }
            p1 = get_basename(av[pNum+1]);
        }

        if(strcmp(av[pNum], "-dev") == 0){
            *selected_device = (cl_uint) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 0) { fprintf(stderr, "Device must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-wpi") == 0){
            *words_per_work_item = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 1) { fprintf(stderr, "words per work item must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-kmer") == 0){
            *kmer_size = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) <= 0) { fprintf(stderr, "Kmer size must be >0\n"); exit(-1); }
        }


        pNum++;

    }   
    
    if(*query==NULL){ fprintf(stderr, "You have to include a query sequence!\n"); exit(-1); }
    strcat(outname, p1);
    strcpy(out, outname);
    
    free(p1);  
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