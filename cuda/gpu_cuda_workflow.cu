
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

uint64_t generate_hits(uint64_t words_at_once, uint64_t * h_pos1, uint64_t * h_pos2, uint64_t * keys_x, uint64_t * keys_y, uint64_t * values_x, uint64_t * values_y);
void read_kmers(uint64_t query_l, char * seq_x, uint64_t * keys_x, uint64_t * values_x);
void init_args(int argc, char ** av, FILE ** query, unsigned * selected_device, FILE ** ref, FILE ** out, unsigned * write);
void perfect_hash_to_word(char * word, uint64_t hash, uint64_t k);
void print_kmers_to_file(uint64_t * keys, uint64_t * values, uint64_t table_size, FILE * fout);
char * get_dirname(char * path);
char * get_basename(char * path);
uint64_t load_seq(FILE * f, char * seq);
uint64_t get_seq_len(FILE * f);
void terror(const char * s){ fprintf(stderr, "\n%s\n", s); exit(-1); }
void Qsort(uint64_t * keys, uint64_t * values, int64_t x, int64_t y);

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
    uint64_t ram_to_be_used = (effective_global_ram) / (2 * (sizeof(Word) + sizeof(char))); //
    uint64_t words_at_once = ram_to_be_used;


    // Allocate words table
    uint64_t * keys, * values, * keys_buf, * values_buf;
    ret = cudaMalloc(&keys, words_at_once * sizeof(uint64_t));
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (1). Error: %d\n", ret); exit(-1); }
    ret = cudaMalloc(&values, words_at_once * sizeof(uint64_t));
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (2). Error: %d\n", ret); exit(-1); }
    ret = cudaMalloc(&keys_buf, words_at_once * sizeof(uint64_t));
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (3). Error: %d\n", ret); exit(-1); }
    ret = cudaMalloc(&values_buf, words_at_once * sizeof(uint64_t));
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (4). Error: %d\n", ret); exit(-1); }
    fprintf(stdout, "[INFO] Allocated %"PRIu64" bytes for hash (you can have %"PRIu64" entries) and their buffers\n", words_at_once * 2 * sizeof(Word), words_at_once);

    // Initialize table
    ret = cudaMemset(keys, 0x0, words_at_once * sizeof(uint64_t));
    ret = cudaMemset(values, 0x0, words_at_once * sizeof(uint64_t));
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

    fprintf(stdout, "[INFO] Loading query\n");
    load_seq(query, query_seq_host);
    fprintf(stdout, "[INFO] Loading reference\n");
    load_seq(ref, ref_seq_host);

    
     
    
    clock_t begin;

    // Pointer to device memory allocating the query sequence
    char * seq_dev_mem = NULL;

    // Allocate memory in device for sequence chunk
    ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char));
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for query sequence in device (Attempted %"PRIu64" bytes). Error: %d\n", words_at_once * sizeof(char), ret); exit(-1); }

    // Allocate memory in host to download kmers and store hits
    uint64_t * dict_x_keys, * dict_x_values, * dict_y_keys, * dict_y_values, * h_pos1, * h_pos2;
    dict_x_keys = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    dict_x_values = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    if(dict_x_keys == NULL || dict_x_values == NULL) { fprintf(stderr, "Allocating for kmer download in query. Error: %d\n", ret); exit(-1); }
    dict_y_keys = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    dict_y_values = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    if(dict_y_keys == NULL || dict_y_values == NULL) { fprintf(stderr, "Allocating for kmer download in ref. Error: %d\n", ret); exit(-1); }
    h_pos1 = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    h_pos2 = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    if(h_pos1 == NULL || h_pos2 == NULL) { fprintf(stderr, "Allocating for hits download. Error: %d\n", ret); exit(-1); }


    ////////////////////////////////////////////////////////////////////////////////
    // Read the query and reference in blocks
    ////////////////////////////////////////////////////////////////////////////////

    FILE * debug = fopen("yo", "wt");
    int split = 0;
    uint64_t pos_in_query = 0, pos_in_ref = 0;
    while(pos_in_query < query_len){

        fprintf(stdout, "[EXECUTING] Running split %d -> (%d%%)\n", split, (int)((100*pos_in_query)/query_len));
        uint64_t items_read = MIN(query_len - pos_in_query, words_at_once);


        ////////////////////////////////////////////////////////////////////////////////
        // Run kmers for query
        ////////////////////////////////////////////////////////////////////////////////
        
        // Load sequence chunk into ram
        ret = cudaMemcpy(seq_dev_mem, &query_seq_host[pos_in_query], items_read, cudaMemcpyHostToDevice);
        if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device. Error: %d\n", ret); exit(-1); }

        // Run kmers
        number_of_blocks = (((items_read - KMER_SIZE + 1)) / (threads_number*4));
        ret = cudaMemset(keys, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
        ret = cudaMemset(values, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
        kernel_register_fast_hash_rotational<<<number_of_blocks, threads_number>>>(keys, values, seq_dev_mem, pos_in_query);
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on query. Error: %d\n", ret); exit(-1); }

        ////////////////////////////////////////////////////////////////////////////////
        // Sort the query kmers
        ////////////////////////////////////////////////////////////////////////////////

        cub::DoubleBuffer<uint64_t> d_keys(keys, keys_buf);
        cub::DoubleBuffer<uint64_t> d_values(values, values_buf);
        void * d_temp_storage = NULL;
        size_t temp_storage_bytes = 0;
        cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, items_read);
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "Bad pre-sorting (1). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

        // Allocate temporary storage
        ret = cudaMalloc(&d_temp_storage, temp_storage_bytes);
        if(ret != cudaSuccess){ fprintf(stderr, "Bad allocating of temp storage for words sorting (1). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

        // Remove this for  debug
        //read_kmers(query_len, query_seq_host, dict_x_keys, dict_x_values);
        //ret = cudaMemcpy(keys, dict_x_keys, items_read*sizeof(uint64_t), cudaMemcpyHostToDevice);
        //ret = cudaMemcpy(values, dict_x_values, items_read*sizeof(uint64_t), cudaMemcpyHostToDevice);
        
        cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, items_read);
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "CUB sorting failed on query. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
        

        // Download kmers        
        ret = cudaMemcpy(dict_x_keys, keys, items_read*sizeof(uint64_t), cudaMemcpyDeviceToHost);
        if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (1). Error: %d\n", ret); exit(-1); }
        ret = cudaMemcpy(dict_x_values, values, items_read*sizeof(uint64_t), cudaMemcpyDeviceToHost);
        if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (2). Error: %d\n", ret); exit(-1); }

        ret = cudaFree(d_temp_storage);
        if(ret != cudaSuccess){ fprintf(stderr, "Bad free of temp storage (1): %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

        pos_in_query += words_at_once;



        ////////////////////////////////////////////////////////////////////////////////
        // Run the reference blocks
        ////////////////////////////////////////////////////////////////////////////////

        while(pos_in_ref < ref_len){

            items_read = MIN(ref_len - pos_in_ref, words_at_once);

            // Load sequence chunk into ram
            ret = cudaMemcpy(seq_dev_mem, &ref_seq_host[pos_in_ref], items_read, cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device. Error: %d\n", ret); exit(-1); }

            // Run kmers
            number_of_blocks = (((items_read - KMER_SIZE + 1)) / (threads_number*4)); 
            ret = cudaMemset(keys, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(values, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            kernel_register_fast_hash_rotational<<<number_of_blocks, threads_number>>>(keys, values, seq_dev_mem, pos_in_ref);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on ref. Error: %d\n", ret); exit(-1); }

            // remove all this             
            
            //read_kmers(ref_len, ref_seq_host, dict_y_keys, dict_y_values);
            //ret = cudaMemcpy(keys, dict_y_keys, items_read*sizeof(uint64_t), cudaMemcpyHostToDevice);
            //ret = cudaMemcpy(values, dict_y_values, items_read*sizeof(uint64_t), cudaMemcpyHostToDevice);


            ////////////////////////////////////////////////////////////////////////////////
            // Sort reference kmers
            ////////////////////////////////////////////////////////////////////////////////

            ret = cudaMemset(keys_buf, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(values_buf, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            cub::DoubleBuffer<uint64_t> d_keys_ref(keys, keys_buf);
            cub::DoubleBuffer<uint64_t> d_values_ref(values, values_buf);
            d_temp_storage = NULL;
            temp_storage_bytes = 0;
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys_ref, d_values_ref, items_read);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Bad pre-sorting (2). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            // Allocate temporary storage
            ret = cudaMalloc(&d_temp_storage, temp_storage_bytes);
            if(ret != cudaSuccess){ fprintf(stderr, "Bad allocating of temp storage for words sorting (2). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys_ref, d_values_ref, items_read);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "CUB sorting failed on ref. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            


            // Download kmers
            ret = cudaMemcpy(dict_y_keys, keys, items_read*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (3). Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(dict_y_values, values, items_read*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (4). Error: %d\n", ret); exit(-1); }

            ret = cudaFree(d_temp_storage);
            if(ret != cudaSuccess){ fprintf(stderr, "Bad free of temp storage (2): %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            pos_in_ref += words_at_once;


            ////////////////////////////////////////////////////////////////////////////////
            // Generate hits for the current split
            ////////////////////////////////////////////////////////////////////////////////

            // SUbstitute kmers for debugggggg

            //read_kmers(query_len, query_seq_host, dict_x_keys, dict_x_values);
            //Qsort(dict_x_keys, dict_x_values, 0, (int64_t) query_len);
            //for(i=0; i<words_at_once; i++) printf("%" PRIu64" %"PRIu64"\n", dict_x_keys[i], dict_x_values[i]);
            //read_kmers(ref_len, ref_seq_host, dict_y_keys, dict_y_values);
            //Qsort(dict_y_keys, dict_y_values, 0, (int64_t) ref_len);
            //for(i=0; i<words_at_once; i++) printf("%" PRIu64" %"PRIu64"\n", dict_x_keys[i], dict_x_values[i]);

            uint64_t n_hits_found = generate_hits(words_at_once, h_pos1, h_pos2, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values);
            
            printf("generated %"PRIu64"\n", n_hits_found);

            fprintf(debug, "All by-Identity Ungapped Fragments (Hits based approach)\n");
            fprintf(debug, "[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>\n");
            fprintf(debug, "SeqX filename : undef\n");
            fprintf(debug, "SeqY filename : undef\n");
            fprintf(debug, "SeqX name : undef\n");
            fprintf(debug, "SeqY name : undef\n");
            fprintf(debug, "SeqX length : %"PRIu64"\n", query_len);
            fprintf(debug, "SeqY length : %"PRIu64"\n", ref_len);
            fprintf(debug, "Min.fragment.length : undef\n");
            fprintf(debug, "Min.Identity : undef\n");
            fprintf(debug, "Tot Hits (seeds) : undef\n");
            fprintf(debug, "Tot Hits (seeds) used: undef\n");
            fprintf(debug, "Total fragments : undef\n");
            fprintf(debug, "========================================================\n");
            fprintf(debug, "Total CSB: 0\n");
            fprintf(debug, "========================================================\n");
            fprintf(debug, "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");

            for(i=0; i<n_hits_found; i++){
                fprintf(debug, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,32,32,32,1.0,1.0,0,0\n", h_pos1[i], h_pos2[i], h_pos1[i]+32, h_pos2[i]+32);
            }
            
        }

        // Restart the reference for every block in query
        pos_in_ref = 0;




        ++split;
    }

    fclose(debug);

    if(write == 1)
    {
        print_kmers_to_file(dict_x_keys, dict_x_values, query_len, out);
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
    
    
    
    
    

    fprintf(stdout, "[INFO] Completed\n");

    fclose(query);
    fclose(ref);
    free(query_seq_host);
    free(ref_seq_host);
    

    return 0;
}

uint64_t generate_hits(uint64_t words_at_once, uint64_t * h_pos1, uint64_t * h_pos2, 
    uint64_t * keys_x, uint64_t * keys_y, uint64_t * values_x, uint64_t * values_y){

    uint64_t id_x = 0, id_y = 0, n_hits_found = 0;
    
    while(id_x < words_at_once || id_y < words_at_once) {

        // Compare
        if(id_x == words_at_once || id_y == words_at_once){ printf("breaking\n");  break; }
        
        
        if(keys_x[id_x] == keys_y[id_y] && values_x[id_x] != 0xFFFFFFFFFFFFFFFF && values_y[id_y] != 0xFFFFFFFFFFFFFFFF) {
            // This is a hit
            //printf("Made hit: ");
            h_pos1[n_hits_found] = values_x[id_x];
            h_pos2[n_hits_found] = values_y[id_y];

            //printf("Matching hash %"PRIu64" with %"PRIu64" @ (%"PRIu64", %"PRIu64")\n", keys_x[id_x], keys_y[id_y], values_x[id_x], values_y[id_y]);

            ++n_hits_found;
            if(n_hits_found == words_at_once){ fprintf(stderr, "Reached maximum limit of hits\n"); }

            ++id_y;

            //printf("next comp is %"PRIu64" with %"PRIu64"\n", keys_x[id_x], keys_y[id_y]);
        }
        else if(keys_x[id_x] < keys_y[id_y]) ++id_x;
        else ++id_y;
    }

    //printf("Generated %"PRIu64" hits \n", n_hits_found);
    return n_hits_found;

}

void read_kmers(uint64_t query_l, char * seq_x, uint64_t * keys_x, uint64_t * values_x){


    memset(keys_x, 0x0, query_l * sizeof(uint64_t));
    memset(values_x, 0xFFFFFFFF, query_l * sizeof(uint64_t));

    uint64_t current_len = 0;
    char c;
    uint64_t word_size = 0, hash = 0;

    while( current_len < query_l ) {

        c = seq_x[current_len++];
            
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            
            //curr_kmer[word_size] = c;
            if(word_size < 32) ++word_size;

            if(c == 'A') { hash = (hash << 2) + 0; }
            if(c == 'C') { hash = (hash << 2) + 1; }
            if(c == 'G') { hash = (hash << 2) + 2; }
            if(c == 'T') { hash = (hash << 2) + 3; }

        }else{ //It can be anything (including N, Y, X ...)

            if(c != '\n' && c != '>') {
                
                hash = 0;
                word_size = 0;
                ++current_len;

            } 
        }

        if(word_size == 32){
            
            keys_x[current_len] = hash;
            values_x[current_len] = current_len;

        }
    }
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
    strcat(outname, ".csv");
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

void print_kmers_to_file(uint64_t * keys, uint64_t * values, uint64_t table_size, FILE * fout){
    
        
    uint64_t i;
    char word[KMER_SIZE+1];
    for(i=0;i<table_size;i++){
        perfect_hash_to_word(word, keys[i], KMER_SIZE);
        word[KMER_SIZE] = '\0';
        fprintf(fout, "#%"PRIu64", %s\n", i, word);
        fprintf(fout, "#%"PRIu64", %"PRIu64" at %" PRIu64"\n", i, keys[i], values[i]);
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

void Qsort(uint64_t * keys, uint64_t * values, int64_t x, int64_t y) {
    
    uint64_t pivote_k, aux_k, aux_val;

    int64_t x1, y1;

    pivote_k = keys[(x+y)/2];

    x1 = x;
    y1 = y;
    
    do { 
        while (pivote_k > keys[x1]) x1++;
        while (pivote_k < keys[y1]) y1--;
        if (x1 < y1) {
            aux_k = keys[x1]; aux_val = values[x1];
            keys[x1] = keys[y1]; values[x1] = values[y1];
            keys[y1] = aux_k; values[y1] = aux_val;

            x1++;
            y1--;
        }
        else if (x1 == y1) x1++;
    } while (x1 <=y1);
    
    if (x<y1) Qsort(keys, values, x, y1);
    if (x1<y) Qsort(keys, values, x1, y);
} 