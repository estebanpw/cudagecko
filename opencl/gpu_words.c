/*

Input fasta sequence must not contain \n's!

*/


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

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

void init_args(int argc, char ** av, FILE ** query, cl_uint * selected_device, ulong * kmers_per_work_item, ulong * kmer_size, FILE ** out, ulong * generate_reverse);
char * get_dirname(char * path);
char * get_basename(char * path);

int main(int argc, char ** argv)
{
    cl_uint selected_device = 0;
    ulong kmer_size = 32;
    ulong generate_reverse = 0;
    ulong kmers_per_work_item = 32;
    FILE * query = NULL, * out = NULL;
    init_args(argc, argv, &query, &selected_device, &kmers_per_work_item, &kmer_size, &out, &generate_reverse);

    char * path_kernels = get_dirname(argv[0]);
    fprintf(stdout, "[INFO] Working on directory %s\n", path_kernels);

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
    size_t work_group_size[3], work_group_size_global;
    cl_int ret;
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

    fprintf(stdout, "[INFO] Using device %d\n", selected_device);
    
    // Create an OpenCL context
    cl_context context = clCreateContext(NULL, 1, &devices[selected_device], NULL, NULL, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not create context. Error: %d\n", ret); exit(-1); }

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, devices[selected_device], 0, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not create command queue. Error: %d\n", ret); exit(-1); }

    

    ////////////////////////////////////////////////////////////////////////////////
    // Make index dictionary
    ////////////////////////////////////////////////////////////////////////////////

    ulong total_ram_asked = 0;

    // Check size of input sequence
    fseek(query, 0, SEEK_END);
    ulong query_len_bytes = (ulong) ftell(query);
    rewind(query);

    // Calculate how much ram we can use for every chunk
    ulong words_to_store = query_len_bytes - kmer_size + 1;
    if(query_len_bytes*sizeof(wordGPU) >= global_device_RAM){
        
        words_to_store = (-(1000*1024*1024)+global_device_RAM)/(sizeof(wordGPU)+1); // This equation finds the maximum chunk size to fit the words structure in ram as well as the sequence chunk (in MB) while leaving some space (negative number)
        fprintf(stdout, "[INFO] Attempting to store %"PRIu64" words\n", words_to_store);
    }

    // Initialize words list
    /*
    Hash_item empty_hash_item;
    memset(&empty_hash_item, 0x0, sizeof(Hash_item));
    empty_hash_item.pos_in_y = 0xFFFFFFFFFFFFFFFF;
    ret = clEnqueueFillBuffer(command_queue, hash_table_mem, (const void *) &empty_hash_item, sizeof(Hash_item), 0, hash_table_size * sizeof(Hash_item), 0, NULL, NULL); 
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not initialize hash table. Error: %d\n", ret); exit(-1); }
    */

    // Allocate memory in host for sequence chunk
    char * query_mem_host = (char *) malloc(words_to_store * sizeof(char));
    if(query_mem_host == NULL){ fprintf(stderr, "Could not allocate host memory (i.e. %"PRIu64" bytes)for query sequence\n", words_to_store); exit(-1); }

    // Load kernel
    FILE * read_kernel; 
    char kernel_temp_path[BUFFER_SIZE];
    kernel_temp_path[0] = '\0';
    strcat(kernel_temp_path, path_kernels);
    switch(generate_reverse){
        case 0: { strcat(kernel_temp_path, "/kernels/kernel_words_f.cl") ; read_kernel = fopen(kernel_temp_path, "r"); fprintf(stdout, "[INFO] Generating only forward strand.\n"); } 
        break;
        case 1: { strcat(kernel_temp_path, "/kernels/kernel_words_r.cl") ; read_kernel = fopen(kernel_temp_path, "r"); fprintf(stdout, "[INFO] Generating forward and reverse strand.\n"); }
        break;
        
        default: { fprintf(stderr, "Bad choice of strand parameter\n"); exit(-1); }
        break;
    }
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
    cl_kernel kernel = clCreateKernel(program, "kernel_words", &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating kernel (1): %d\n", ret); exit(-1); }

    // Set working size
    size_t local_item_size = CORES_PER_COMPUTE_UNIT * 8; // Number of work items in a work group
    size_t global_item_size;

    // Allocate host memory to retrieve words
    wordGPU * words = (wordGPU *) calloc(words_to_store, sizeof(wordGPU));
    if(words == NULL){ fprintf(stderr, "Error allocating host words\n"); exit(-1); }

    // Read the input query in chunks
    int split = 0;
    uint64_t items_read = 0;
    query_len_bytes = 0;
    while(!feof(query)){

        // Load sequence chunk into ram
        items_read = fread(query_mem_host, sizeof(char), words_to_store, &query[0]);
        query_len_bytes += items_read;


        // Allocate words list
        cl_mem words_list = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, words_to_store * sizeof(wordGPU), words, &ret);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for word list in device. Error: %d\n", ret); exit(-1); }

        total_ram_asked = words_to_store * sizeof(wordGPU);

        // Allocate memory in device for sequence chunk
        cl_mem query_mem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, items_read * sizeof(char), query_mem_host, &ret);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory (i.e. %"PRIu64" bytes) for query sequence in device. Error: %d\n", items_read, ret); exit(-1); }

        total_ram_asked += items_read * sizeof(char);

        // Set global working sizes
        fprintf(stdout, "[INFO] Split #%d: %"PRIu64"\n", split, items_read);
    
        global_item_size = ((items_read - kmer_size + 1)) / kmers_per_work_item;
        global_item_size = global_item_size - (global_item_size % local_item_size); // Make it evenly divisable

        fprintf(stdout, "[INFO] Work items: %"PRIu64". Work groups: %"PRIu64". Total K-mers to be computed %"PRIu64"\n", (uint64_t) global_item_size, (uint64_t)(global_item_size/local_item_size), global_item_size * kmers_per_work_item);

        // Load parameters
        parameter_word params = {kmer_size, items_read, (ulong) global_item_size, (ulong) kmers_per_work_item, query_len_bytes - items_read};    
        cl_mem params_mem = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(parameter_word), &params, &ret);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for kmer sizes variable in device. Error: %d\n", ret); exit(-1); }

        total_ram_asked += sizeof(parameter_word);

        // Set the arguments of the kernel
        
        ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&words_list);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (1): %d\n", ret); exit(-1); }
        ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&params_mem);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (2): %d\n", ret); exit(-1); }
        ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&query_mem);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (3): %d\n", ret); exit(-1); }
        

        fprintf(stdout, "[INFO] Executing the kernel on split %d, global RAM usage: %"PRIu64" MB\n", split++, total_ram_asked/(1024*1024));

        // Execute the OpenCL kernel on the lists
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, 
                &global_item_size, &local_item_size, 0, NULL, NULL);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Error enqueuing kernel (1): %d\n", ret); exit(-1); }

        // Wait for kernel to finish
        ret = clFlush(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
        ret = clFinish(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event: %d\n", ret); exit(-1); }

        // Transfer words
        ret = clEnqueueReadBuffer(command_queue, words_list, CL_TRUE, 0, words_to_store*sizeof(wordGPU), words, 0, NULL, NULL);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not read from buffer: %d\n", ret); exit(-1); }

        // Before deallocating, write to disk the words
        for(i=0; i<words_to_store; i++){
            if(words[i].pos != 0xFFFFFFFFFFFFFFFF) fwrite(&words[i], sizeof(wordGPU), 1, out);
        }

        // Deallocation & cleanup for next round
        ret = clReleaseMemObject(words_list); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (4)\n"); exit(-1); }
        ret = clReleaseMemObject(query_mem); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (5)\n"); exit(-1); }
        ret = clReleaseMemObject(params_mem); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (6)\n"); exit(-1); }
            
    }

    fclose(query);
    fclose(out);
    free(query_mem_host);


    // Wait for kernel to finish
    ret = clFlush(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
    ret = clFinish(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event: %d\n", ret); exit(-1); }
    ret = clReleaseKernel(kernel); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (3)\n"); exit(-1); }
    ret = clReleaseProgram(program); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (4)\n"); exit(-1); }

    fprintf(stdout, "[INFO] Completed processing query splits\n");

    

    return 0;
}


void init_args(int argc, char ** av, FILE ** query, cl_uint * selected_device, ulong * kmers_per_work_item, ulong * kmer_size, FILE ** out, ulong * generate_reverse){
    
    int pNum = 0;
    char * p1;
    char outname[2048];
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           words -query [file]\n");
            fprintf(stdout, "OPTIONAL:\n");
            
            fprintf(stdout, "           -dev        [Integer: d>=0] Selects the device to be used\n");
            fprintf(stdout, "           -kmer       [Integer: k>=1] Size of K-mer to be used\n");
            fprintf(stdout, "           -kwi        [Integer: k>=1] Number of kmers per work item to be used\n");
            fprintf(stdout, "           --r         Generate the reverse strand\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }

        if(strcmp(av[pNum], "-query") == 0){
            *query = fopen(av[pNum+1], "rt");
            if(*query==NULL){ fprintf(stderr, "Could not open query file\n"); exit(-1); }
            p1 = get_basename(av[pNum+1]);
        }

        if(strcmp(av[pNum], "-dev") == 0){
            *selected_device = (cl_uint) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 0) { fprintf(stderr, "Device must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-kwi") == 0){
            *kmers_per_work_item = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 1) { fprintf(stderr, "Kmers per work item size must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-kmer") == 0){
            *kmer_size = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) <= 0) { fprintf(stderr, "Kmer size must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "--r") == 0){
            *generate_reverse = 1;
        }


        pNum++;

    }   
    
    if(*query==NULL){ fprintf(stderr, "You have to include a query sequence!\n"); exit(-1); }
    strcat(outname, p1);
    if(*generate_reverse == 1) strcat(outname, "-words-r.bin"); else strcat(outname, "-words-f.bin");
    *out = fopen(outname, "wb");
    if(*out == NULL){ fprintf(stderr, "Could not open output file\n"); exit(-1); }
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