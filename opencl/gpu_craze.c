
// Standard utilities and common systems includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <math.h>
#include "commonFunctions.h"
#include "structs.h"


#define BUFFER_SIZE 2048
#define MAX_KERNEL_SIZE BUFFER_SIZE*100
#define CORES_PER_COMPUTE_UNIT 32
#define DIMENSION 1000

char * get_dirname(char * path);
char * get_basename(char * path);
void init_args(int argc, char ** av, FILE ** query, FILE ** ref, FILE ** hits, FILE ** frags, cl_uint * selected_device, ulong * kmer_size, ulong * min_len, float * min_sim, long double * expected_value, unsigned char * reverse, FILE ** info, ulong * n_seqs_y);
char * read_seq(FILE * f, uint64_t * l);
void print_hit(hitGPU * h);
void print_fragment(struct FragFile * f);

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
    cl_ulong kmer_size = 32, min_len = 30, n_seqs_y = 0;
    float min_sim = 30.0;
    long double min_expected_value = 0.01;
    size_t work_group_size[3], work_group_size_global;
    cl_int ret;
    unsigned char reverse = 'f';
    char * path_kernels = get_dirname(argv[0]);
    fprintf(stdout, "[INFO] Working on directory %s\n", path_kernels);

    FILE * query = NULL, * ref = NULL, * hits, * frags = NULL, * info = NULL;

    //init_args(argc, argv, &query, &ref, &hits, &frags, &selected_device, &kmer_size, &min_len, &min_sim, &min_expected_value, &reverse, &info, &n_seqs_y);
    fprintf(stdout, "Working on parameters: k:%"PRIu64", l:%"PRIu64", s:%f, e:%Le\n", kmer_size, min_len, min_sim, min_expected_value);

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
    selected_device = 3;
    fprintf(stdout, "[INFO] Using device %d\n", selected_device);
    
    // Create an OpenCL context
    cl_context context = clCreateContext(NULL, 1, &devices[selected_device], NULL, NULL, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not create context. Error: %d\n", ret); exit(-1); }

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, devices[selected_device], 0, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not create command queue. Error: %d\n", ret); exit(-1); }

    ////////////////////////////////////////////////////////////////////////////////
    // Load sequences into GPU ram
    ////////////////////////////////////////////////////////////////////////////////

    // Check sequence and hits sizes
    uint64_t n_huge = 10 * 1000 * 1000;
    char * seq_x = (char *) malloc(n_huge);
    char * seq_y = (char *) malloc(n_huge);
    if(seq_x == NULL || seq_y == NULL){ fprintf(stderr, "Could not allocate memory for sequences\n"); exit(-1); }

    // Load into GPU
    cl_mem cl_seq_x = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, n_huge * sizeof(char), seq_x, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for seq_x in device. Error: %d\n", ret); exit(-1); }

    cl_mem cl_seq_y = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, n_huge * sizeof(char), seq_y, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for seq_y in device. Error: %d\n", ret); exit(-1); }

    // No longer needed, its in GPU
    free(seq_x);
    free(seq_y);

    ////////////////////////////////////////////////////////////////////////////////
    // Calculate how many hits and frags we can load at once
    ////////////////////////////////////////////////////////////////////////////////


    // Load kernel
    FILE * read_kernel; 
    char kernel_temp_path[BUFFER_SIZE];
    kernel_temp_path[0] = '\0';
    strcat(kernel_temp_path, path_kernels);
    
    strcat(kernel_temp_path, "/kernels/kernel_craze.cl"); 
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
    cl_kernel kernel = clCreateKernel(program, "kernel_craze", &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Error creating kernel (1): %d\n", ret); exit(-1); }
    
    // Set working size
    size_t local_item_size = 256; // Number of work items in a work group
    uint64_t two_D_kernel_size = n_huge * n_huge;
    size_t global_item_size = two_D_kernel_size + (local_item_size - two_D_kernel_size % local_item_size);

    fprintf(stdout, "[INFO] Work group size: %"PRIu64", work-items: %"PRIu64"\n", local_item_size, global_item_size);

    // Allocate frags on GPU
    uint64_t in_once = 1000 * 1000 * 20, computed = 0; // 640 MB
    // Allocate space on host for writing the output
    reduced_frag * frag_write = (reduced_frag *) calloc(in_once, sizeof(reduced_frag));
    if(frag_write == NULL){ fprintf(stderr, "Could not allocate memory for host frags"); exit(-1); }
    cl_mem cl_frags = clCreateBuffer(context, CL_MEM_READ_WRITE, in_once * sizeof(reduced_frag), NULL, &ret);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Could not allocate memory for frags in device. Error: %d\n", ret); exit(-1); }

    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&cl_seq_x);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (1): %d\n", ret); exit(-1); }

    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&cl_seq_y);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (2): %d\n", ret); exit(-1); }

    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&cl_frags);
    if(ret != CL_SUCCESS){ fprintf(stderr, "Bad setting of param (4): %d\n", ret); exit(-1); }

    while(computed < two_D_kernel_size){
       
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, 
            &global_item_size, &local_item_size, 0, NULL, NULL);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Error enqueuing kernel (1): %d\n", ret); exit(-1); }

        // Wait for kernel to finish
        ret = clFlush(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad flush of event: %d\n", ret); exit(-1); }
        ret = clFinish(command_queue); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad finish of event: %d\n", ret); exit(-1); }

        computed += in_once;
        // Simulate read
        ret = clEnqueueReadBuffer(command_queue, cl_frags, CL_TRUE, 0, in_once*sizeof(reduced_frag), frag_write, 0, NULL, NULL);
        if(ret != CL_SUCCESS){ fprintf(stderr, "Could not read from buffer: %d\n", ret); exit(-1); }

        //See how long it takes
        if(computed % 100000000000 == 0) printf("Computed a thousand\n");

    }


    
    ret = clReleaseKernel(kernel); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (3)\n"); exit(-1); }
    ret = clReleaseProgram(program); if(ret != CL_SUCCESS){ fprintf(stderr, "Bad free (4)\n"); exit(-1); }

    fprintf(stdout, "[INFO] Completed processing query splits\n");

    // Write rest of info regarding sequences
    
    return 0;
}


void init_args(int argc, char ** av, FILE ** query, FILE ** ref, FILE ** hits, FILE ** frags, cl_uint * selected_device, ulong * kmer_size, ulong * min_len, float * min_sim, long double * expected_value, unsigned char * reverse, FILE ** info, ulong * n_seqs_y){
    
    int pNum = 0;
    char * p1, * p2;
    char outname[2048];
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           gpu_frags -q [file] -r [file] -hits [file]\n");
            fprintf(stdout, "OPTIONAL:\n");
            
            fprintf(stdout, "           -dev        [Integer: d>=0] Selects the device to be used\n");
            fprintf(stdout, "           -kmer       [Integer: k>=1] Size of K-mer to be used\n");
            fprintf(stdout, "           -l          [Integer: k>=1] Minimum length\n");
            fprintf(stdout, "           -s          [Float : 0<k<100] Minimum similarity\n");
            fprintf(stdout, "           -e          [Float  : 0<=e] Minimum attained e-value\n");
            fprintf(stdout, "           --r         It is a reverse comparison\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }

        if(strcmp(av[pNum], "-q") == 0){
            *query = fopen(av[pNum+1], "rt");
            if(*query==NULL){ fprintf(stderr, "Could not open input query sequence file\n"); exit(-1); }
            p1 = get_basename(av[pNum+1]);
        }
        if(strcmp(av[pNum], "-r") == 0){
            *ref = fopen(av[pNum+1], "rt");
            if(*ref==NULL){ fprintf(stderr, "Could not open input ref sequence file\n"); exit(-1); }
            p2 = get_basename(av[pNum+1]);
        }

        if(strcmp(av[pNum], "-hits") == 0){
            *hits = fopen(av[pNum+1], "rb");
            if(*hits==NULL){ fprintf(stderr, "Could not open input hits file\n"); exit(-1); }

            char temp_name[2048]; temp_name[0] = '\0';
            strcat(temp_name, av[pNum+1]);
            strcat(temp_name, ".nseqs");

            FILE * read_seqs_y = fopen(temp_name, "rt");
            if(read_seqs_y == NULL){ fprintf(stderr, "Could not open %s\n", temp_name); exit(-1); }
            if(1 != fscanf(read_seqs_y, "%"PRIu64, n_seqs_y)){ fprintf(stderr, "Could not read number of y-seqs\n"); exit(-1); }
            fprintf(stdout, "Read y-seqs %"PRIu64" successfully from %s\n", *n_seqs_y, temp_name);
            fclose(read_seqs_y);
        }

        if(strcmp(av[pNum], "--r") == 0){
            *reverse = 'r';
        }

        if(strcmp(av[pNum], "-dev") == 0){
            *selected_device = (cl_uint) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 0) { fprintf(stderr, "Device must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-kmer") == 0){
            *kmer_size = (ulong) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) <= 0) { fprintf(stderr, "Kmer size must be >0\n"); exit(-1); }
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
    
    if(*query != NULL && *ref != NULL){
        strcat(outname, p1);
        strcat(outname, "-");
        strcat(outname, p2);
        if(*reverse == 'f') strcat(outname, "-f"); else strcat(outname, "-r");
        strcat(outname, ".frags");
        fprintf(stdout, "Writing frags to %s.\n", outname);
        *frags = fopen(outname, "wb");
        strcat(outname, ".INF");
        *info = fopen(outname, "wt");
        if(*info == NULL) { fprintf(stderr, "Could not open info file\n"); exit(-1); }
        fprintf(*info, "All by-Identity Ungapped Fragments (Hits based approach)\n");
        fprintf(*info, "[Abr.98/Apr.2010/Dec.2011/Jun.2018 -- BITLAB\n");
        fprintf(*info, "SeqX filename        : %s\n", p1);
        fprintf(*info, "SeqY filename        : %s\n", p2);
        fprintf(*info, "SeqX name            : %s\n", p1);
	    fprintf(*info, "SeqY name            : %s\n", p2);
    }

    if(*query==NULL || *ref==NULL || *hits==NULL || *frags==NULL){ fprintf(stderr, "You have to include a query, reference, hits and frags output!\n"); exit(-1); }
    
    free(p1);
    free(p2);

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



char * read_seq(FILE * f, uint64_t * l) {
    char c;
    uint64_t lon = 0, k = 0;
    uint64_t SIZE = 0;
    uint64_t i = 0, r = 0;
    char * seq = NULL;
    char * buffer = NULL;


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
            seq[lon++] = '*';
            while (c != '\n') {
                if ((feof(f) &&  i >= r ))
                    return 0;
                c = buffered_fgetc(buffer, &i, &r, f);
            }
            //break;
        }
        if (isupper(c))
            seq[lon++] = c;
        if (c == '*') {
            seq[lon++] = c;
        }
        c = buffered_fgetc(buffer, &i, &r, f);
    }

    free(buffer);
    seq[lon] = 0x00;

    *l = lon;

    return seq;
}

void print_hit(hitGPU * h){
    fprintf(stdout, "HIT::pos(x,y) (%"PRIu64", %"PRIu64") from seq(x,y) (%"PRIu64", %"PRIu64")", h->pos_x, h->pos_y, h->seq_x, h->seq_y);
}

void print_fragment(struct FragFile * f){
    
    fprintf(stdout, "FRAG::(%"PRIu64", %"PRIu64") to (%"PRIu64", %"PRIu64"): [%"PRIu64"]-[%"PRIu64"] %c LEN[%"PRIu64"] IDE[%"PRIu64"] SIM:%f\n", f->xStart, f->yStart, f->xEnd, f->yEnd, f->seqX, f->seqY, f->strand, f->length, f->ident, f->similarity);
}