
// Standard utilities and common systems includes
#include "kernels.cuh"
#include "cub/cub.cuh"

#define BUFFER_SIZE 2048
#define CORES_PER_COMPUTE_UNIT 32
#define KMER_SIZE 32
//#define DIMENSION 1000


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

void init_args(int argc, char ** av, FILE ** query, unsigned * selected_device, FILE ** ref, FILE ** out, unsigned * write);
void perfect_hash_to_word(char * word, uint64_t hash, uint64_t k);
void print_kmers_to_file(Word * table_mem, uint64_t table_size, FILE * fout);
char * get_dirname(char * path);
char * get_basename(char * path);
uint64_t load_seq(FILE * f, char * seq);
uint64_t get_seq_len(FILE * f);
void terror(const char * s){ fprintf(stderr, "\n%s\n", s); exit(-1); }

int main(int argc, char ** argv)
{
    uint64_t i;
    unsigned selected_device = 0, write = 0;
    FILE * query = NULL, * ref = NULL, * out = NULL;
    init_args(argc, argv, &query, &selected_device, &ref, &out, &write);

    ////////////////////////////////////////////////////////////////////////////////
    // Get info of devices
    ////////////////////////////////////////////////////////////////////////////////

    int ret_num_devices;
    unsigned compute_units;
    uint64_t local_device_RAM, global_device_RAM;
    int work_group_dimensions[3], work_group_size_local;
    int ret;
    
    // Query how many devices there are
    if(cudaSuccess != (ret = cudaGetDeviceCount(&ret_num_devices))){ fprintf(stderr, "Failed to query number of devices\n"); exit(-1); }

    cudaDeviceProp device;

    for(i=0; i<ret_num_devices; i++){
        if( cudaSuccess != (ret = cudaGetDeviceProperties(&device, i))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }

        fprintf(stdout, "\tDevice [%"PRIu64"]: %s\n", i, device.name);
        global_device_RAM = device.totalGlobalMem;
        fprintf(stdout, "\t\tGlobal mem   : %"PRIu64" (%"PRIu64" MB)\n", (uint64_t) global_device_RAM, (uint64_t) global_device_RAM / (1024*1024));
        local_device_RAM = device.sharedMemPerBlock;
        fprintf(stdout, "\t\tLocal mem    : %"PRIu64" (%"PRIu64" KB)\n", (uint64_t) local_device_RAM, (uint64_t) local_device_RAM / (1024));
        compute_units = device.multiProcessorCount;
        fprintf(stdout, "\t\tCompute units: %"PRIu64"\n", (uint64_t) compute_units);
        work_group_size_local = device.maxThreadsPerBlock;
        fprintf(stdout, "\t\tMax work group size: %d\n", work_group_size_local);
        memcpy(work_group_dimensions, device.maxThreadsDim, 3*sizeof(int));
        fprintf(stdout, "\t\tWork size dimensions: (%d, %d, %d)\n", work_group_dimensions[0], work_group_dimensions[1], work_group_dimensions[2]);
        fprintf(stdout, "\t\tWarp size: %d\n", device.warpSize);
        fprintf(stdout, "\t\tGrid dimensions: (%d, %d, %d)\n", device.maxGridSize[0], device.maxGridSize[1], device.maxGridSize[2]);
    }
    //selected_device = 3; // REMOVE --- ONLY FOR TESTING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    if( cudaSuccess != (ret = cudaSetDevice(selected_device))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }
    fprintf(stdout, "[INFO] Using device %d\n", selected_device);

    if( cudaSuccess != (ret = cudaGetDeviceProperties(&device, selected_device))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }
    global_device_RAM = device.totalGlobalMem;

    
    
    ////////////////////////////////////////////////////////////////////////////////
    // Make index dictionary
    ////////////////////////////////////////////////////////////////////////////////

    

    // Calculate how much ram we can use for every chunk
    uint64_t effective_global_ram =  (global_device_RAM - 100*1024*1024); //Minus 100 MBs
    uint64_t ram_to_be_used = (effective_global_ram) / (sizeof(Word) + sizeof(char)); //
    uint64_t words_at_once = ram_to_be_used;


    // Allocate words table
    Word * dictionary_dev = NULL;
    ret = cudaMalloc(&dictionary_dev, words_at_once * sizeof(Word));
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device. Error: %d\n", ret); exit(-1); }
    fprintf(stdout, "[INFO] Allocated %"PRIu64" bytes for hash %"PRIu64" entries\n", words_at_once * sizeof(Word), words_at_once);

    // Initialize table
    ret = cudaMemset(dictionary_dev, 0x0, words_at_once * sizeof(Word));
    if(ret != cudaSuccess){ fprintf(stderr, "Could not initialize words table. Error: %d\n", ret); exit(-1); }

    
    // Set working size
    size_t threads_number = 32;
    size_t number_of_blocks;
    //cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte); // NOTICE: MAXWELL ignores this--

    

    // Inspect shared memory configuration
    cudaSharedMemConfig shared_mem_conf;
    ret = cudaDeviceGetSharedMemConfig(&shared_mem_conf);
    if(ret != cudaSuccess){ fprintf(stdout, "[WARNING] Could not get shared memory configuration. Error: %d\n", ret); }
    else { fprintf(stdout, "[INFO] Shared memory configuration is: %s\n", (shared_mem_conf == cudaSharedMemBankSizeFourByte) ? ("4 bytes") : ("8 bytes")); }

    // Load DNA sequences
    uint64_t query_len = get_seq_len(query);
    uint64_t ref_len = get_seq_len(ref);

    char * query_seq_host = (char *) malloc(query_len * sizeof(char));
    char * ref_seq_host = (char *) malloc(ref_len * sizeof(char));

    if(query_seq_host == NULL || ref_seq_host == NULL) terror("Could not allocate memory for sequences in host");

    load_seq(query, query_seq_host);
    load_seq(ref, ref_seq_host);
     
    
    clock_t begin;

    // Pointer to device memory allocating the query sequence
    char * seq_dev_mem = NULL;

    // Allocate memory in device for sequence chunk
    ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char));
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for query sequence in device (Attempted %"PRIu64" bytes). Error: %d\n", words_at_once * sizeof(char), ret); exit(-1); }


    // Read the input query in chunks
    int split = 0;
    uint64_t pos_in_query = 0, pos_in_ref = 0;
    while(pos_in_query < query_len){

        uint64_t items_read = MIN(query_len - pos_in_query, words_at_once);
        
        // Load sequence chunk into ram
        ret = cudaMemcpy(seq_dev_mem, &query_seq_host[pos_in_query], items_read, cudaMemcpyHostToDevice);
        if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device. Error: %d\n", ret); exit(-1); }

        // Run kmers
        number_of_blocks = (((items_read - KMER_SIZE + 1)) / (threads_number*4)); printf("[1] Processing: %"PRIu64"\r", number_of_blocks*threads_number*4); // Blocks
        ret = cudaMemset(dictionary_dev, 0x0, words_at_once * sizeof(Word));
        kernel_register_fast_hash_rotational<<<number_of_blocks, threads_number>>>(dictionary_dev, seq_dev_mem, pos_in_query);
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on query. Error: %d\n", ret); exit(-1); }

        // Sort kmers

        /*

        // Declare, allocate, and initialize device-accessible pointers for sorting data
        uint64_t num_items;          // e.g., 7
        uint64_t * d_key_buf;         // e.g., [8, 6, 7, 5, 3, 0, 9]
        uint64_t * d_key_alt_buf;     // e.g., [        ...        ]
        uint64_t * d_value_buf;       // e.g., [0, 1, 2, 3, 4, 5, 6]
        uint64_t * d_value_alt_buf;   // e.g., [        ...        ]

        cub::DoubleBuffer<uint64_t> d_keys(d_key_buf, d_key_alt_buf);
        cub::DoubleBuffer<uint64_t> d_values(d_value_buf, d_value_alt_buf);
        void * d_temp_storage = NULL;
        size_t temp_storage_bytes = 0;
        cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, num_items);
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "CUB sorting failed on query. Error: %d\n", ret); exit(-1); }
        */

        // Download kmers
        Word * dict_x = (Word *) malloc(items_read*sizeof(Word));
        ret = cudaMemcpy(dict_x, dictionary_dev, items_read*sizeof(Word), cudaMemcpyDeviceToHost);
        if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers. Error: %d\n", ret); exit(-1); }

        pos_in_query += words_at_once;

        if(ret != cudaSuccess){ fprintf(stderr, "Bad finish of indexing kernel: %d : %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

        /*

        while(pos_in_ref < ref_len){

            items_read = MIN(ref_len - pos_in_ref, words_at_once);

            // Load sequence chunk into ram
            ret = cudaMemcpy(seq_dev_mem, &query_seq_host[pos_in_query], items_read, cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device. Error: %d\n", ret); exit(-1); }

            // Run kmers
            number_of_blocks = (((items_read - KMER_SIZE + 1)) / (threads_number*4)); printf("[2] Processing: %"PRIu64"\r", number_of_blocks*threads_number*4); // Blocks
            ret = cudaMemset(dictionary_dev, 0x0, words_at_once * sizeof(Word));
            kernel_register_fast_hash_rotational<<<number_of_blocks, threads_number>>>(dictionary_dev, seq_dev_mem, pos_in_query);
            ret = cudaDeviceSynchronize();

            // Download kmers
            Word * dict_y = (Word *) malloc(items_read*sizeof(Word));
            ret = cudaMemcpy(dict_y, dictionary_dev, items_read*sizeof(Word), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers. Error: %d\n", ret); exit(-1); }

            pos_in_ref += words_at_once;
            
        }

        */

        // Free dict_x since it will be created again
        free(dict_x);
        

        
        // Deallocation & cleanup for next round

        ret = cudaFree(query_mem);
        if(ret != cudaSuccess){ fprintf(stderr, "Bad free of query memory in indexing: %d\n", ret); exit(-1); }

        ++split;
        fprintf(stdout, "[INFO] Split %d\n", split);
    }

    if(write == 1)
    {
        print_kmers_to_file(dictionary_dev, query_len, out);
        //print_kmers_to_file_paused(table_mem, query_len_bytes);
    }

    /*
    uint64_t * debug = (uint64_t *) malloc(table_size*sizeof(uint64_t));
    ret = cudaMemcpy(debug, table_mem, table_size*sizeof(uint64_t), cudaMemcpyDeviceToHost); 
    if(ret != cudaSuccess){ fprintf(stderr, "DEBUG. Error: %d\n", ret); exit(-1); }
    
    for(i=0;i<12;i++){
        fprintf(stdout, "#%"PRIu64": %"PRIu64"\n", i, debug[i]);
    } 
    */
    
    
    
    
    

    

    fclose(query);
    fclose(ref);
    free(query_seq_host);
    free(ref_seq_host);

    ret = cudaFree(dictionary_q_f);
    if(ret != cudaSuccess){ fprintf(stderr, "Bad free of k-mer table: %d\n", ret); exit(-1); }
    

    return 0;
}


void init_args(int argc, char ** av, FILE ** query, unsigned * selected_device, FILE ** ref, FILE ** out, unsigned * write){
    
    int pNum = 0;
    char * p1 = NULL, * p2 = NULL;
    char outname[2048]; outname[0] = '\0';
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           GECKOCUDA -query [file] -ref [file] -device [device]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           --write     Enables writing output\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }
        if(strcmp(av[pNum], "--write") == 0){
            *write = 1;
        }

        if(strcmp(av[pNum], "-query") == 0){
            *query = fopen(av[pNum+1], "rt");
            if(*query==NULL){ fprintf(stderr, "Could not open query file\n"); exit(-1); }
            p1 = get_basename(av[pNum+1]);
        }
        
        if(strcmp(av[pNum], "-ref") == 0){
            *ref = fopen(av[pNum+1], "rt");
            if(*ref==NULL){ fprintf(stderr, "Could not open reference file\n"); exit(-1); }
            p2 = get_basename(av[pNum+1]);
        }

        if(strcmp(av[pNum], "-dev") == 0){
            *selected_device = (unsigned) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 0) { fprintf(stderr, "Device must be >0\n"); exit(-1); }
        }

        pNum++;

    }   
    
    if(*query==NULL || *ref==NULL){ fprintf(stderr, "You have to include a query and a reference sequence!\n"); exit(-1); }
    strcat(outname, p1);
    strcat(outname, "-");
    strcat(outname, p2);
    strcat(outname, ".kmers");
    *out = fopen(outname, "wt");
    if(*out == NULL){ fprintf(stderr, "Could not open output file\n"); exit(-1); }
    if(p1 != NULL) free(p1);
    if(p2 != NULL) free(p2);   
}

void perfect_hash_to_word(char * word, uint64_t hash, uint64_t k){
    /*
    int64_t jIdx = (int64_t) (k-1), upIdx = 31;
    uint64_t v;
    while(jIdx >= 0){
        v = (uint64_t) floor(hash / (pow(4, jIdx)));
        if(v == 0){ word[upIdx--] = (char) 'A'; hash -= ((uint64_t) pow(4, (uint64_t) jIdx) * 0); }
        if(v == 1){ word[upIdx--] = (char) 'C'; hash -= ((uint64_t) pow(4, (uint64_t) jIdx) * 1); }
        if(v == 2){ word[upIdx--] = (char) 'G'; hash -= ((uint64_t) pow(4, (uint64_t) jIdx) * 2); }
        if(v == 3){ word[upIdx--] = (char) 'T'; hash -= ((uint64_t) pow(4, (uint64_t) jIdx) * 3); }
        
        if(jIdx == 0) break;
        --jIdx;
    }
    */
    
    int64_t jIdx = (int64_t) (k-1), upIdx = 0;
    uint64_t v;
    while(jIdx >= 0){
        v = (uint64_t) floor(hash / (pow(4, jIdx)));
        if(v == 0){ word[upIdx++] = (char) 'A'; hash -= ((uint64_t) pow(4, (uint64_t) jIdx) * 0); }
        if(v == 1){ word[upIdx++] = (char) 'C'; hash -= ((uint64_t) pow(4, (uint64_t) jIdx) * 1); }
        if(v == 2){ word[upIdx++] = (char) 'G'; hash -= ((uint64_t) pow(4, (uint64_t) jIdx) * 2); }
        if(v == 3){ word[upIdx++] = (char) 'T'; hash -= ((uint64_t) pow(4, (uint64_t) jIdx) * 3); }
        
        if(jIdx == 0) break;
        --jIdx;
    }
}

void print_kmers_to_file(Word * table_mem, uint64_t table_size, FILE * fout){
    
    Word * debug = (Word *) malloc(table_size*sizeof(Word));
    int ret = cudaMemcpy(debug, table_mem, table_size*sizeof(Word), cudaMemcpyDeviceToHost);
    if(ret != cudaSuccess){ fprintf(stderr, "DEBUG. Error: %d\n", ret); exit(-1); }
    
    uint64_t i;
    char word[KMER_SIZE+1];
    for(i=0;i<table_size;i++){
        perfect_hash_to_word(word, debug[i].hash, KMER_SIZE);
        word[KMER_SIZE] = '\0';
        fprintf(fout, "#%"PRIu64", %s\n", i, word);
        fprintf(fout, "#%"PRIu64", %"PRIu64" at %" PRIu64"\n", i, debug[i].hash, debug[i].pos);
    } 
    fclose(fout);
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

uint64_t get_seq_len(FILE * f) {
    char c = '\0';
    uint64_t l = 0;

    while(!feof(f)){
        c = getc(f);

        if(c == '>'){
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            ++l;
        }
    }


    rewind(f);
    return l;
}

uint64_t load_seq(FILE * f, char * seq) {
    char c = '\0';
    uint64_t l = 0;

    while(!feof(f)){
        c = getc(f);

        if(c == '>'){
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            seq[l] = c;
            ++l;
        }
    }


    rewind(f);
    return l;
}