
// Standard utilities and common systems includes
#include "kernels.cuh"
#include "cpu_functions.c"
#include "cub/cub.cuh"

#include <cuda_profiler_api.h>

//#define DIMENSION 1000


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////


void print_header(FILE * out, uint32_t query_len, uint32_t ref_len);
uint64_t memory_allocation_chooser(uint64_t total_memory);

int main(int argc, char ** argv)
{
    uint32_t i, min_length = 64, max_frequency = 0;
    float factor = 0.125;
    int fast = 0; // sensitive is default
    unsigned selected_device = 0;
    FILE * query = NULL, * ref = NULL, * out = NULL;
    init_args(argc, argv, &query, &selected_device, &ref, &out, &min_length, &fast, &max_frequency, &factor);

    clock_t start = clock();
    /*Do something*/
    clock_t end = clock();

    ////////////////////////////////////////////////////////////////////////////////
    // Get info of devices
    ////////////////////////////////////////////////////////////////////////////////

    int ret_num_devices;
    //unsigned compute_units;
    uint64_t global_device_RAM;
    //int work_group_size_local;
    int ret;
    
    // Query how many devices there are
    if(cudaSuccess != (ret = cudaGetDeviceCount(&ret_num_devices))){ fprintf(stderr, "Failed to query number of devices\n"); exit(-1); }

    cudaDeviceProp device;

    for(i=0; i<ret_num_devices; i++){
        if( cudaSuccess != (ret = cudaGetDeviceProperties(&device, i))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }

        fprintf(stdout, "\tDevice [%"PRIu32"]: %s\n", i, device.name);
        global_device_RAM = device.totalGlobalMem;
        fprintf(stdout, "\t\tGlobal mem   : %"PRIu64" (%"PRIu64" MB)\n", global_device_RAM, global_device_RAM / (1024*1024));
        //compute_units = device.multiProcessorCount;
        //fprintf(stdout, "\t\tCompute units: %"PRIu64"\n", (uint64_t) compute_units);
        //work_group_size_local = device.maxThreadsPerBlock;
        //fprintf(stdout, "\t\tMax work group size: %d\n", work_group_size_local);
        //fprintf(stdout, "\t\tWork size dimensions: (%d, %d, %d)\n", work_group_dimensions[0], work_group_dimensions[1], work_group_dimensions[2]);
        //fprintf(stdout, "\t\tWarp size: %d\n", device.warpSize);
        //fprintf(stdout, "\t\tGrid dimensions: (%d, %d, %d)\n", device.maxGridSize[0], device.maxGridSize[1], device.maxGridSize[2]);
    }
    //selected_device = 3; // REMOVE --- ONLY FOR TESTING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    if( cudaSuccess != (ret = cudaSetDevice(selected_device))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }
    fprintf(stdout, "[INFO] Using device %d\n", selected_device);

    if( cudaSuccess != (ret = cudaGetDeviceProperties(&device, selected_device))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }
    global_device_RAM = device.totalGlobalMem;

    


    // Calculate how much ram we can use for every chunk
    uint64_t effective_global_ram =  (global_device_RAM - memory_allocation_chooser(global_device_RAM)); //Minus 100 to 300 MBs for other stuff

    // We will do the one-time alloc here
    // i.e. allocate a pool once and used it manually

    char * data_mem;
    ret = cudaMalloc(&data_mem, effective_global_ram * sizeof(char)); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pool memory in device. Error: %d\n", ret); exit(-1); }

    // One workload depends on number of words (words + sortwords + generate hits)
    // The other one depends on number of hits (sort hits + filterhits + frags)
    /*
    uint32_t factor = 24;
    uint32_t ram_to_be_used = (effective_global_ram) / (factor); // update this with max usage
    uint32_t words_at_once = ram_to_be_used;
    if(fast == 0) words_at_once = words_at_once/4;
    uint32_t max_hits = effective_global_ram - (words_at_once*8 + words_at_once*4);
    */
    if(fast == 1) factor = 0.5;
    uint64_t bytes_for_words = (factor * effective_global_ram); // 512 MB for words
    uint64_t words_at_once = bytes_for_words / (8+8+4+4); 
    uint64_t max_hits = (effective_global_ram - bytes_for_words) / (8+8+4+4); // The rest for allocating hits



    fprintf(stdout, "[INFO] You can have %"PRIu64" MB for words (i.e. %"PRIu64" words), and %"PRIu64" MB for hits (i.e. %"PRIu64" hits)\n", 
        bytes_for_words / (1024*1024), words_at_once, (effective_global_ram - bytes_for_words) / (1024*1024), max_hits);

    fprintf(stdout, "[INFO] Filtering at a minimum length of %"PRIu32" bps\n", min_length);
    if(fast == 1) 
        fprintf(stdout, "[INFO] Running on fast mode (some repetitive seeds will be skipped)\n");
    else
        fprintf(stdout, "[INFO] Running on sensitive mode (ALL seeds are computed [mf:%"PRIu32"])\n", max_frequency);

    

    
    // Set working size
    size_t threads_number = 32;
    size_t number_of_blocks;
    //cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte); // NOTICE: MAXWELL ignores this--

    

    // Inspect shared memory configuration
    //cudaSharedMemConfig shared_mem_conf;
    //ret = cudaDeviceGetSharedMemConfig(&shared_mem_conf);
    //if(ret != cudaSuccess){ fprintf(stdout, "[WARNING] Could not get shared memory configuration. Error: %d\n", ret); }
    //else { fprintf(stdout, "[INFO] Shared memory configuration is: %s\n", (shared_mem_conf == cudaSharedMemBankSizeFourByte) ? ("4 bytes") : ("8 bytes")); }

    // Load DNA sequences
    //uint32_t query_len = get_seq_len(query);
    //uint32_t ref_len = get_seq_len(ref);


    // Load faster
    fseek(query, 0L, SEEK_END);
    uint32_t coarse_query_len = (uint32_t) ftell(query);
    rewind(query);
    char * s_buffer = (char *) malloc(coarse_query_len * sizeof(char)); if(s_buffer == NULL) {fprintf(stderr, "Bad loading buffer (1)\n"); exit(-1);}
    char * pro_q_buffer = (char *) malloc(coarse_query_len * sizeof(char));  if(pro_q_buffer == NULL) {fprintf(stderr, "Bad loading buffer (2)\n"); exit(-1);}
    uint32_t read_bytes = (uint32_t) fread(s_buffer, 1, coarse_query_len, query); if(read_bytes < coarse_query_len) {fprintf(stderr, "Bad bytes reading (1)\n"); exit(-1);}
    uint32_t query_len = from_ram_load(s_buffer, pro_q_buffer, coarse_query_len);

    free(s_buffer);

    fseek(ref, 0L, SEEK_END);
    uint32_t coarse_ref_len = (uint32_t) ftell(ref);
    rewind(ref);
    s_buffer = (char *) malloc(coarse_ref_len * sizeof(char));  if(s_buffer == NULL) {fprintf(stderr, "Bad loading buffer (3)\n"); exit(-1);}
    char * pro_r_buffer = (char *) malloc(coarse_ref_len * sizeof(char));  if(pro_r_buffer == NULL) {fprintf(stderr, "Bad loading buffer (4)\n"); exit(-1);}
    read_bytes = (uint32_t) fread(s_buffer, 1, coarse_ref_len, ref); if(read_bytes < coarse_ref_len) {fprintf(stderr, "Bad bytes reading (2)\n"); exit(-1);}
    uint32_t ref_len = from_ram_load(s_buffer, pro_r_buffer, coarse_ref_len);

    free(s_buffer);


    fprintf(stdout, "[INFO] Qlen: %"PRIu32"; Rlen: %"PRIu32"\n", query_len, ref_len);

    // Check that sequence length complies
    if(MAX(query_len, ref_len) >= 2147483648){
        fprintf(stdout, "[WARNING] !!!!!!!!!!!!!!!!!!!!!!\n");
        fprintf(stdout, "[WARNING] PLEASE READ CAREFULLY\n");
        fprintf(stdout, "[WARNING] THE INPUT SEQUENCES ARE TOO LONG (MAX LEN 2147483648)\n");
        fprintf(stdout, "[WARNING] THE PROGRAM WILL CONTINUE TO WORK BUT MIGHT PRODUCE SOME ERRORS\n");
        fprintf(stdout, "[WARNING] THESE CAN APPEAR PARTICULARLY IN THE LIMITS OF THE SEQUENCE\n");
        fprintf(stdout, "[WARNING] CHECK THIS ISSUE ON A DOTPLOT\n");
        fprintf(stdout, "[WARNING] !!!!!!!!!!!!!!!!!!!!!!\n");
    }

    // Create variables to load up sequences using PINNED MEMORY to increase transfer rate
    //char * query_seq_host = (char *) malloc(query_len * sizeof(char));
    //char * ref_seq_host = (char *) malloc(ref_len * sizeof(char));
    //char * ref_rev_seq_host = (char *) malloc(ref_len * sizeof(char));

    // How about one big alloc (save ~3 seconds on mallocs)
    char * host_pinned_mem, * base_ptr_pinned;
    
    uint64_t pinned_bytes_on_host = words_at_once * (sizeof(uint64_t) * (2) + sizeof(uint32_t) * (2));
    pinned_bytes_on_host = pinned_bytes_on_host + max_hits * (sizeof(uint64_t) * (1) + sizeof(uint32_t) * (8));
    pinned_bytes_on_host = pinned_bytes_on_host + sizeof(char) * (query_len + 2*ref_len);
    pinned_bytes_on_host += 1024*1024; // Adding 1 MB for the extra padding in the realignments
    uint64_t pinned_address_checker = 0;


    fprintf(stdout, "[INFO] Allocating on host %"PRIu64" bytes (i.e. %"PRIu64" MBs)\n", pinned_bytes_on_host, pinned_bytes_on_host / (1024*1024));
    ret = cudaHostAlloc(&host_pinned_mem, pinned_bytes_on_host, cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for pool. Error: %d\n", ret); exit(-1); }
    

    char * query_seq_host, * ref_seq_host, * ref_rev_seq_host;
    base_ptr_pinned = (char *) &host_pinned_mem[0];
    query_seq_host = (char *) &host_pinned_mem[0];
    pinned_address_checker = realign_address(pinned_address_checker + query_len, 4);

    ref_seq_host = (char *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + ref_len, 4);

    ref_rev_seq_host = (char *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + ref_len, 4);
   
    printf("Reverse starts at %p\n", ref_rev_seq_host); 


    
    //ret = cudaHostAlloc(&query_seq_host, query_len * sizeof(char), cudaHostAllocMapped); 
    //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for query_seq host. Error: %d\n", ret); exit(-1); }
    //ret = cudaHostAlloc(&ref_seq_host, ref_len * sizeof(char), cudaHostAllocMapped); 
    //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for ref_seq host. Error: %d\n", ret); exit(-1); }
    //ret = cudaHostAlloc(&ref_rev_seq_host, ref_len * sizeof(char), cudaHostAllocMapped); 
    //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for reverse ref_seq host. Error: %d\n", ret); exit(-1); }
    

    

    //cudaHostAlloc((void**)&a,n*sizeof(a),cudaHostAllocDefault);

    if(query_seq_host == NULL || ref_seq_host == NULL || ref_rev_seq_host == NULL) terror("Could not allocate memory for sequences in host");

    ////////////////////////////////////////////////////////////////////////////////
    // Read sequences and reverse the reference
    ////////////////////////////////////////////////////////////////////////////////

    // Create streams to allow concurrent copy and execute

    uint32_t n_streams = 4;
    cudaStream_t streams[n_streams];

    for(i=0; i<n_streams; i++) cudaStreamCreate(&streams[i]);

    // Pointer to device memory allocating the query sequence, reference and reversed reference

    fprintf(stdout, "[INFO] Loading query\n");
    //load_seq(query, query_seq_host);
    memcpy(query_seq_host, pro_q_buffer, query_len);
    fprintf(stdout, "[INFO] Loading reference\n");

    start = clock();
    //load_seq(ref, ref_seq_host);
    memcpy(ref_seq_host, pro_r_buffer, ref_len);
    fprintf(stdout, "[INFO] Reversing reference\n");

    free(pro_q_buffer);
    free(pro_r_buffer);
    
    /*
    ret = cudaMalloc(&seq_dev_mem_reverse_aux, ref_len * sizeof(char)); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for reverse reference sequence in device (Attempted %"PRIu32" bytes) at reversing. Error: %d\n", (uint32_t) (ref_len * sizeof(char)), ret); exit(-1); }
    */

    // ## POINTER SECTION 0
    uint64_t address_checker = 0;
    
    char * ptr_seq_dev_mem_aux = &data_mem[0];
    address_checker = realign_address(address_checker + ref_len, 128);
    char * ptr_seq_dev_mem_reverse_aux = &data_mem[address_checker];
    address_checker = realign_address(address_checker + ref_len, 4);

    char * ptr_reverse_write[n_streams];

    uint32_t chunk_size = ref_len/n_streams + 4; // Always force one divisor more cause of decimal loss
    if(chunk_size % 128 != 0){ chunk_size += 128 - (chunk_size % 128); printf("ENTER\n"); } // Since sequence starts at 0, making the chunks multiple of 128 guarantees 100% GL efficiency


    ptr_reverse_write[0] = ptr_seq_dev_mem_reverse_aux;
    for(i=1; i<n_streams-1; i++){
        ptr_reverse_write[i] = ptr_reverse_write[i-1] + chunk_size;
    }
    ptr_reverse_write[n_streams-1] = ptr_reverse_write[n_streams-2] + chunk_size;
    uintptr_t where_ends = (uintptr_t) (ptr_reverse_write[n_streams-1] + MIN(chunk_size, ref_len - chunk_size*(n_streams-1)));
    if(where_ends % 128 != 0) ptr_reverse_write[n_streams-1] += 128 - (where_ends % 128);

    //printf("Starts %p Endings %p\n", ptr_reverse_write[0], ptr_reverse_write[0]+chunk_size);
    //printf("Starts %p Endings %p\n", ptr_reverse_write[1], ptr_reverse_write[1]+chunk_size);
    //printf("Starts %p Endings %p\n", ptr_reverse_write[2], ptr_reverse_write[2]+chunk_size);
    //printf("Starts %p Endings %p\n", ptr_reverse_write[3], ptr_reverse_write[3]+MIN(chunk_size, ref_len - chunk_size*(n_streams-1)));
   
 
    threads_number = 128;
    //threads_number = 32;
    if(ref_len > 1024){
    

        number_of_blocks = chunk_size/threads_number + threads_number; // Same
        if(number_of_blocks % threads_number != 0){ number_of_blocks += 1; printf("ENTER 2\n");}
        uint32_t offset; //, inverse_offset;

        
        for(i=0; i<n_streams; i++)
        {
            offset = chunk_size * i;
            //inverse_offset = (chunk_size * (i+1) < ref_len) ? ref_len - chunk_size * (i+1) : 0;

            ret = cudaMemcpyAsync(ptr_seq_dev_mem_aux + offset, ref_seq_host + offset, MIN(chunk_size, ref_len - offset), cudaMemcpyHostToDevice, streams[i]);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy reference sequence to device for reversing. Error: %d\n", ret); exit(-1); }

            //printf("On the other hand, the load align %p\n", ptr_seq_dev_mem_aux + offset);

            //cudaProfilerStart();
            kernel_reverse_complement<<<number_of_blocks, threads_number, 0, streams[i]>>>(ptr_seq_dev_mem_aux + offset, ptr_reverse_write[i], MIN(chunk_size, ref_len - offset));
            //cudaProfilerStop();

        }

        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "Could not compute reverse on reference. Error: %d\n", ret); exit(-1); }

        // Perform copy from the nstreams
        uint32_t t_copy_rev = 0;
        for(i=0; i<n_streams; i++){

            ret = cudaMemcpy(ref_rev_seq_host + t_copy_rev, ptr_reverse_write[n_streams - (i+1)], MIN(chunk_size, ref_len - chunk_size*(n_streams-1-i)), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy reference sequence to device for reversing. Error: %d\n", ret); exit(-1); }
            t_copy_rev += MIN(chunk_size, ref_len - chunk_size*(n_streams-1-i));
        }



    }else{


        number_of_blocks = (ref_len)/threads_number + 1;


        uintptr_t integer_ptr = (uintptr_t) ptr_seq_dev_mem_reverse_aux;
        uint64_t ptr_end = ((uint64_t) integer_ptr) + ref_len;
    
        if(ptr_end % 128 != 0) ptr_seq_dev_mem_reverse_aux += 128 - (ptr_end % 128);




        ret = cudaMemcpy(ptr_seq_dev_mem_aux, ref_seq_host, ref_len, cudaMemcpyHostToDevice);
        if(ret != cudaSuccess){ fprintf(stderr, "Could not copy reference sequence to device for reversing. Error: %d\n", ret); exit(-1); }

        //cudaProfilerStart();
        kernel_reverse_complement<<<number_of_blocks, threads_number>>>(ptr_seq_dev_mem_aux, ptr_seq_dev_mem_reverse_aux, ref_len);
        //cudaProfilerStop();
        
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "Could not compute reverse on reference. Error: %d\n", ret); exit(-1); }

        ret = cudaMemcpy(ref_rev_seq_host, ptr_seq_dev_mem_reverse_aux, ref_len, cudaMemcpyDeviceToHost);
        if(ret != cudaSuccess){ fprintf(stderr, "Could not copy reference sequence to device for reversing. Error: %d\n", ret); exit(-1); }
    }
    
    threads_number = 32;



    /*

    }
    */


    ret = cudaDeviceSynchronize();
    end = clock();

    //cudaFree(seq_dev_mem_aux);
    //cudaFree(seq_dev_mem_reverse_aux);

    //for(i=0; i<ref_len; i++){
    //    if(isupper(ref_rev_seq_host[i]) != ref_rev_seq_host[i]) {
    //        printf("Found first at %u and it is %.32s\n", i, &ref_rev_seq_host[i]); break;
    //    }
    //}


    // Print some info
#ifdef SHOWTIME
    fprintf(stdout, "[INFO] rev comp t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif
    fprintf(stdout, "[INFO] Showing start of reference sequence:\n");
    
    fprintf(stdout, "\t(Begin ref)%.64s\n", ref_seq_host);
    fprintf(stdout, "\t(Begin rev)%.64s\n", ref_rev_seq_host);
    fprintf(stdout, "\t(End   ref)%.64s\n", &ref_seq_host[ref_len-64]);
    fprintf(stdout, "\t(End   rev)%.64s\n", &ref_rev_seq_host[ref_len-64]);

     
    //fprintf(stdout, "\t(Full que)%.*s\n", query_len, query_seq_host);
    //fprintf(stdout, "\t(Full ref)%.*s\n", ref_len, ref_seq_host);
    //fprintf(stdout, "\t(Full rev)%.*s\n", ref_len, ref_rev_seq_host);


    // Write header to CSV
    print_header(out, query_len, ref_len);

    ////////////////////////////////////////////////////////////////////////////////
    // Allocation of pointers
    ////////////////////////////////////////////////////////////////////////////////


    // Allocate memory in host to download kmers and store hits
    
    uint64_t * dict_x_keys, * dict_y_keys; // Keys are hashes (64-b), values are positions (32)
    uint32_t * dict_x_values, * dict_y_values;

    pinned_address_checker = realign_address(pinned_address_checker, 8);
    dict_x_keys = (uint64_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + words_at_once * sizeof(uint64_t), 4);

    dict_x_values = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + words_at_once * sizeof(uint32_t), 8);

    dict_y_keys = (uint64_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + words_at_once * sizeof(uint64_t), 4);

    dict_y_values = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + words_at_once * sizeof(uint32_t), 4);
    // These depends on the number of words
    /*
    //dict_x_keys = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    ret = cudaHostAlloc(&dict_x_keys, words_at_once * sizeof(uint64_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for dict_x_keys. Error: %d\n", ret); exit(-1); }
    //dict_x_values = (uint32_t *) malloc(words_at_once*sizeof(uint32_t));
    ret = cudaHostAlloc(&dict_x_values, words_at_once * sizeof(uint32_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for dict_x_values. Error: %d\n", ret); exit(-1); }
    //if(dict_x_keys == NULL || dict_x_values == NULL) { fprintf(stderr, "Allocating for kmer download in query. Error: %d\n", ret); exit(-1); }
    //dict_y_keys = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    ret = cudaHostAlloc(&dict_y_keys, words_at_once * sizeof(uint64_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for dict_y_keys. Error: %d\n", ret); exit(-1); }
    //dict_y_values = (uint32_t *) malloc(words_at_once*sizeof(uint32_t));
    ret = cudaHostAlloc(&dict_y_values, words_at_once * sizeof(uint32_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for dict_y_values. Error: %d\n", ret); exit(-1); }
    //if(dict_y_keys == NULL || dict_y_values == NULL) { fprintf(stderr, "Allocating for kmer download in ref. Error: %d\n", ret); exit(-1); }
    */

    

    
    // These are now depending on the number of hits
    Hit * hits;
    uint32_t * filtered_hits_x, * filtered_hits_y;
    /*
    //Hit * hits = (Hit *) malloc(max_hits*sizeof(Hit));
    Hit * hits;
    ret = cudaHostAlloc(&hits, max_hits * sizeof(Hit), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for hits. Error: %d\n", ret); exit(-1); }

    //uint32_t * filtered_hits_x = (uint32_t *) malloc(max_hits*sizeof(uint32_t));
    uint32_t * filtered_hits_x;
    ret = cudaHostAlloc(&filtered_hits_x, max_hits * sizeof(uint32_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for filtered_hits_x. Error: %d\n", ret); exit(-1); }

    //uint32_t * filtered_hits_y = (uint32_t *) malloc(max_hits*sizeof(uint32_t));
    uint32_t * filtered_hits_y;
    ret = cudaHostAlloc(&filtered_hits_y, max_hits * sizeof(uint32_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for filtered_hits_y. Error: %d\n", ret); exit(-1); }
    */

    pinned_address_checker = realign_address(pinned_address_checker, 8);
    hits = (Hit *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(Hit), 4);

    filtered_hits_x = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 4);

    filtered_hits_y = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 4);


    uint32_t * host_left_offset,  * host_right_offset, * ascending_numbers, * indexing_numbers;
    uint64_t * diagonals;
    /*
    // These are for the device
    uint32_t * device_filt_hits_x, * device_filt_hits_y, * left_offset, * right_offset;

    //uint32_t * host_left_offset = (uint32_t *) malloc(max_hits*sizeof(uint32_t));
    uint32_t * host_left_offset;
    ret = cudaHostAlloc(&host_left_offset, max_hits * sizeof(uint32_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for host_left_offset. Error: %d\n", ret); exit(-1); }

    //uint32_t * host_right_offset = (uint32_t *) malloc(max_hits*sizeof(uint32_t));
    uint32_t * host_right_offset;
    ret = cudaHostAlloc(&host_right_offset, max_hits * sizeof(uint32_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for host_right_offset. Error: %d\n", ret); exit(-1); }

    //if(host_left_offset == NULL || host_right_offset == NULL) terror("Could not allocate host offsets");
    //if(hits == NULL || filtered_hits_x == NULL || filtered_hits_y == NULL) terror("Could not allocate hits");

    //uint64_t * diagonals = (uint64_t *) malloc(max_hits*sizeof(uint64_t));
    uint64_t * diagonals;
    ret = cudaHostAlloc(&diagonals, max_hits * sizeof(uint64_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for diagonals. Error: %d\n", ret); exit(-1); }

    //uint64_t * device_diagonals, * device_diagonals_buf;
    //uint32_t * device_hits, * device_hits_buf; // These will actually be just indices to redirect the hits sorting

    //uint32_t * ascending_numbers = (uint32_t *) malloc(max_hits*sizeof(uint32_t)); for(i=0; i<max_hits; i++) ascending_numbers[i] = i;
    uint32_t * ascending_numbers;
    ret = cudaHostAlloc(&ascending_numbers, max_hits * sizeof(uint32_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for ascending numbers. Error: %d\n", ret); exit(-1); }
    for(i=0; i<max_hits; i++) ascending_numbers[i] = i;

    //uint32_t * indexing_numbers = (uint32_t *) malloc(max_hits*sizeof(uint32_t));
    uint32_t * indexing_numbers;
    ret = cudaHostAlloc(&indexing_numbers, max_hits * sizeof(uint32_t), cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for indexing numbers. Error: %d\n", ret); exit(-1); }
    
    //if(hits == NULL) { fprintf(stderr, "Allocating for hits download. Error: %d\n", ret); exit(-1); }
    */

    pinned_address_checker = realign_address(pinned_address_checker, 4);
    host_left_offset = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 4);

    host_right_offset = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 8);

    diagonals = (uint64_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint64_t), 4);

    ascending_numbers = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 4);
    for(i=0; i<max_hits; i++) ascending_numbers[i] = i;

    indexing_numbers = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 4);


    
    //printf("ALOHAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA------------------------ remove thisssssssssssssssssssssssssssssssssssssssssss\n");
    //cudaFree(data_mem);


    ////////////////////////////////////////////////////////////////////////////////
    // Read the query and reference in blocks
    ////////////////////////////////////////////////////////////////////////////////


    int split = 0;
    uint32_t pos_in_query = 0, pos_in_ref = 0;
    while(pos_in_query < query_len){


        /*
        printf("ALOHAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA------------------------ remove thisssssssssssssssssssssssssssssssssssssssssss\n");
        // This would not be here, just testing
        ret = cudaMalloc(&data_mem, effective_global_ram * sizeof(char)); 
        if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pool memory in device. Error: %d\n", ret); exit(-1); }
        */

        address_checker = 0;

        // Allocate memory in device for sequence chunk
        // We have to this here since later on we will have to free all memory to load the hits
        //ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char));
        //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for query sequence in device (Attempted %"PRIu32" bytes). Error: %d\n", (uint32_t) (words_at_once * sizeof(char)), ret); exit(-1); }

        // Allocate words table
        //ret = cudaMalloc(&keys, words_at_once * sizeof(uint64_t));
        //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (1). Error: %d\n", ret); exit(-1); }
        //ret = cudaMalloc(&values, words_at_once * sizeof(uint32_t));
        //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (2). Error: %d\n", ret); exit(-1); }

        start = clock();

        // ## POINTER SECTION 1
        char * ptr_seq_dev_mem = &data_mem[0];
        char * base_ptr = ptr_seq_dev_mem;
        address_checker = realign_address(address_checker + words_at_once, 8);
        
        uint64_t * ptr_keys = (uint64_t *) (base_ptr + address_checker); // We have to realign because of the arbitrary length of the sequence chars
        address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 4);

        uint32_t * ptr_values = (uint32_t *) (base_ptr + address_checker); 
        address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 4);
        

        fprintf(stdout, "[EXECUTING] Running split %d -> (%d%%)[%u,%u]\n", split, (int)((100*(uint64_t)pos_in_query)/(uint64_t)query_len), pos_in_query, pos_in_ref);

        uint32_t items_read_x = MIN(query_len - pos_in_query, words_at_once);


        ////////////////////////////////////////////////////////////////////////////////
        // Run kmers for query
        ////////////////////////////////////////////////////////////////////////////////
        
        // Load sequence chunk into ram

        //ret = cudaMemcpy(seq_dev_mem, &query_seq_host[pos_in_query], items_read_x, cudaMemcpyHostToDevice);
        ret = cudaMemcpy(ptr_seq_dev_mem, &query_seq_host[pos_in_query], items_read_x, cudaMemcpyHostToDevice);
        if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device. Error: %d\n", ret); exit(-1); }

        // Run kmers
        ret = cudaMemset(ptr_keys, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
        ret = cudaMemset(ptr_values, 0xFFFFFFFF, words_at_once * sizeof(uint32_t));
        
                
        //number_of_blocks = (items_read_x - KMER_SIZE + 1)/threads_number;
        //kernel_index_global32<<<number_of_blocks, threads_number>>>(ptr_keys, ptr_values, ptr_seq_dev_mem, pos_in_query);


    
        number_of_blocks = (items_read_x - KMER_SIZE + 1)/(64);


        if(number_of_blocks != 0)
        {
            //cudaProfilerStart();
            kernel_index_global32<<<number_of_blocks, 64>>>(ptr_keys, ptr_values, ptr_seq_dev_mem, pos_in_query);
            //kernel_index_global32_advanced<<<number_of_blocks, 64>>>(ptr_keys, ptr_values, (uchar4 *) ptr_seq_dev_mem, pos_in_query);
            //cudaProfilerStop();
        
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on query. Error: %d\n", ret); exit(-1); }
        }

        //cudaFree(seq_dev_mem);
        //cudaFree(data_mem);

        // FOR DEBUG
        // Copy kmers to local
        
        /*
        uint64_t * kmers = (uint64_t *) malloc(words_at_once * sizeof(uint64_t));
        uint64_t * poses = (uint64_t *) malloc(words_at_once * sizeof(uint64_t));
        ret = cudaMemcpy(kmers, keys, items_read_x*sizeof(uint64_t), cudaMemcpyDeviceToHost);
        ret = cudaMemcpy(poses, values, items_read_x*sizeof(uint64_t), cudaMemcpyDeviceToHost);
        FILE * anything8 = fopen("kmers", "a");
        for(i=0; i<words_at_once; i++){
            fprintf(anything8, "%"PRIu64" %"PRIu64" %"PRIu64"\n", i, poses[i], kmers[i]);
        }
        free(kmers); free(poses);
        fclose(anything8);
        */
        
        end = clock();
#ifdef SHOWTIME
        fprintf(stdout, "[INFO] words Q t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

        ////////////////////////////////////////////////////////////////////////////////
        // Sort the query kmers
        ////////////////////////////////////////////////////////////////////////////////

        // Notice --------- Now that we are usiung pooled memory
        // I have not "freed" the part corresponding to the sequence (ptr_seq_dev_mem)
        // ANd thus next points build upon that
        // But thats no problem because it is a small fraction of memory

        //ret = cudaMalloc(&keys_buf, words_at_once * sizeof(uint64_t));
        //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (3). Error: %d\n", ret); exit(-1); }
        //ret = cudaMalloc(&values_buf, words_at_once * sizeof(uint32_t));
        //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (4). Error: %d\n", ret); exit(-1); }


        // ## POINTER SECTION 2

        start = clock();
        
        address_checker = realign_address(address_checker, 8);
        uint64_t * ptr_keys_buf = (uint64_t *) (base_ptr + address_checker);
        address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 4);

        uint32_t * ptr_values_buf = (uint32_t *) (base_ptr + address_checker); // Each alloc adds on top of the previous one
        address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 4);

        //cub::DoubleBuffer<uint64_t> d_keys(keys, keys_buf);
        //cub::DoubleBuffer<uint32_t> d_values(values, values_buf);

        cub::DoubleBuffer<uint64_t> d_keys(ptr_keys, ptr_keys_buf);
        cub::DoubleBuffer<uint32_t> d_values(ptr_values, ptr_values_buf);

        void * d_temp_storage = NULL;
        size_t temp_storage_bytes = 0;
        cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, items_read_x);
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "Bad pre-sorting (1). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

        // Allocate temporary storage
        ret = cudaMalloc(&d_temp_storage, temp_storage_bytes);
        if(ret != cudaSuccess){ fprintf(stderr, "Bad allocating of temp storage for words sorting (1). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

        // Remove this for  debug
        //read_kmers(query_len, query_seq_host, dict_x_keys, dict_x_values);
        //ret = cudaMemcpy(keys, dict_x_keys, items_read_x*sizeof(uint64_t), cudaMemcpyHostToDevice);
        //ret = cudaMemcpy(values, dict_x_values, items_read_x*sizeof(uint64_t), cudaMemcpyHostToDevice);
        
        cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, items_read_x);
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "CUB sorting failed on query. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
        

        // Download sorted kmers
        ret = cudaMemcpyAsync(dict_x_keys, ptr_keys_buf, items_read_x*sizeof(uint64_t), cudaMemcpyDeviceToHost, streams[0]);
        if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (1). Error: %d\n", ret); exit(-1); }
        ret = cudaMemcpyAsync(dict_x_values, ptr_values_buf, items_read_x*sizeof(uint32_t), cudaMemcpyDeviceToHost, streams[1]);
        if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (2). Error: %d\n", ret); exit(-1); }

        // Print hits for debug
        //for(i=0; i<items_read_x; i++){
        //    fprintf(out, "%"PRIu64"\n", dict_x_values[i]);
        //}

        ret = cudaFree(d_temp_storage);
        if(ret != cudaSuccess){ fprintf(stderr, "Bad free of temp storage (1): %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
        
        end = clock();
#ifdef SHOWTIME
        fprintf(stdout, "[INFO] sort words Q t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

        pos_in_query += words_at_once;

        //cudaFree(keys);
        //cudaFree(values);
        //cudaFree(keys_buf);
        //cudaFree(values_buf);

        ////////////////////////////////////////////////////////////////////////////////
        // Run the reference blocks
        ////////////////////////////////////////////////////////////////////////////////

        /*
        printf("ALOHAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA------------------------ remove thisssssssssssssssssssssssssssssssssssssssssss\n");
        cudaFree(data_mem);

        printf("ALOHAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA------------------------ remove thisssssssssssssssssssssssssssssssssssssssssss\n");
        // This would not be here, just testing
        ret = cudaMalloc(&data_mem, effective_global_ram * sizeof(char)); 
        if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pool memory in device. Error: %d\n", ret); exit(-1); }
        */

        while(pos_in_ref < ref_len){

            ////////////////////////////////////////////////////////////////////////////////
            // FORWARD strand in the reference
            ////////////////////////////////////////////////////////////////////////////////

            start = clock();

            uint32_t items_read_y = MIN(ref_len - pos_in_ref, words_at_once);

            //ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char));

            // Allocate words table
            //ret = cudaMalloc(&keys, words_at_once * sizeof(uint64_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (1). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&values, words_at_once * sizeof(uint32_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (2). Error: %d\n", ret); exit(-1); }

            // ## POINTER SECTION 3
            ptr_seq_dev_mem = &data_mem[0];
            base_ptr = ptr_seq_dev_mem;
            address_checker = realign_address(words_at_once, 8);

            ptr_keys = (uint64_t *) (base_ptr + address_checker); // We have to realign because of the arbitrary length of the sequence chars
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 4);

            ptr_values = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 4);
            

            // Load sequence chunk into ram
            ret = cudaMemcpy(ptr_seq_dev_mem, &ref_seq_host[pos_in_ref], items_read_y, cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device. Error: %d\n", ret); exit(-1); }

            // Run kmers
            ret = cudaMemset(ptr_keys, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(ptr_values, 0xFFFFFFFF, words_at_once * sizeof(uint32_t));
            //number_of_blocks = (((items_read_y - KMER_SIZE + 1)) / (threads_number*4)); 
            //kernel_register_fast_hash_rotational<<<number_of_blocks, threads_number>>>(keys, values, seq_dev_mem, pos_in_ref);
            
            number_of_blocks = ((items_read_y - KMER_SIZE + 1))/(64);
            if(number_of_blocks != 0)
            {

                kernel_index_global32<<<number_of_blocks, 64>>>(ptr_keys, ptr_values, ptr_seq_dev_mem, pos_in_ref);
            
                //number_of_blocks = ((items_read_y - KMER_SIZE + 1))/threads_number;
                //kernel_index_global32<<<number_of_blocks, threads_number>>>(ptr_keys, ptr_values, ptr_seq_dev_mem, pos_in_ref);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on ref. Error: %d\n", ret); exit(-1); }

            }

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] words R t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif


            //cudaFree(seq_dev_mem);

            

            ////////////////////////////////////////////////////////////////////////////////
            // Sort reference FORWARD kmers
            ////////////////////////////////////////////////////////////////////////////////

            //ret = cudaMalloc(&keys_buf, words_at_once * sizeof(uint64_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (3). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&values_buf, words_at_once * sizeof(uint32_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (4). Error: %d\n", ret); exit(-1); }

            // ## POINTER SECTION 4

            start = clock();
        
            address_checker = realign_address(address_checker, 8);
            ptr_keys_buf = (uint64_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 4);

            ptr_values_buf = (uint32_t *) (base_ptr + address_checker); // Each alloc adds on top of the previous one
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 4);

            ret = cudaMemset(ptr_keys_buf, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(ptr_values_buf, 0xFFFFFFFF, words_at_once * sizeof(uint32_t));

            cub::DoubleBuffer<uint64_t> d_keys_ref(ptr_keys, ptr_keys_buf);
            cub::DoubleBuffer<uint32_t> d_values_ref(ptr_values, ptr_values_buf);

            d_temp_storage = NULL;
            temp_storage_bytes = 0;

            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys_ref, d_values_ref, items_read_y);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Bad pre-sorting (2). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            // Allocate temporary storage
            ret = cudaMalloc(&d_temp_storage, temp_storage_bytes);
            if(ret != cudaSuccess){ fprintf(stderr, "Bad allocating of temp storage for words sorting (2). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys_ref, d_values_ref, items_read_y);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "CUB sorting failed on ref. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            


            // Download sorted reference kmers
            ret = cudaMemcpy(dict_y_keys, ptr_keys_buf, items_read_y*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (3). Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(dict_y_values, ptr_values_buf, items_read_y*sizeof(uint32_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (4). Error: %d\n", ret); exit(-1); }

            ret = cudaFree(d_temp_storage);
            if(ret != cudaSuccess){ fprintf(stderr, "Bad free of temp storage (2): %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            pos_in_ref += words_at_once;

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] sort words R t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif


            ////////////////////////////////////////////////////////////////////////////////
            // Generate FORWARD hits for the current split
            ////////////////////////////////////////////////////////////////////////////////

            
            //read_kmers(query_len, query_seq_host, dict_x_keys, dict_x_values);
            //Qsort(dict_x_keys, dict_x_values, 0, (int64_t) query_len);
            //for(i=0; i<words_at_once; i++) printf("%" PRIu64" %"PRIu64"\n", dict_x_keys[i], dict_x_values[i]);
            //read_kmers(ref_len, ref_seq_host, dict_y_keys, dict_y_values);
            //Qsort(dict_y_keys, dict_y_values, 0, (int64_t) ref_len);
            //for(i=0; i<words_at_once; i++) printf("%" PRIu64" %"PRIu64"\n", dict_y_keys[i], dict_y_values[i]);

            
            //cudaFree(keys);
            //cudaFree(values);
            //cudaFree(keys_buf);
            //cudaFree(values_buf);

            start = clock();

            uint32_t n_hits_found;
            if(fast == 1)
                n_hits_found = generate_hits_fast(max_hits, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len);
            else
                n_hits_found = generate_hits_sensitive(max_hits, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len, max_frequency);
            

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] hits Q-R t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            fprintf(stdout, "[INFO] Generated %"PRIu32" hits on split %d -> (%d%%)[%u,%u]{%u,%u}\n", n_hits_found, split, (int)((100*MIN((uint64_t)pos_in_ref, (uint64_t)ref_len))/(uint64_t)ref_len), pos_in_query, pos_in_ref, items_read_x, items_read_y);

            // Print hits for debug
            //for(i=0; i<n_hits_found; i++){
            //    fprintf(stdout, "%"PRIu32"\n", dict_x_values[i]);
            //}
            //for(i=0; i<n_hits_found; i++){
                //printf("%"PRIu64"\n", diagonals[i]);
                //if(hits[i].p1 > 368000 && hits[i].p2 < 390000)
                    //fprintf(out, "Frag,d:%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,32,32,32,1.0,1.0,0,0\n", (int64_t) hits[i].p1 - (int64_t) hits[i].p2, hits[i].p1, hits[i].p2, hits[i].p1+32, hits[i].p2+32);
            //}

            ////////////////////////////////////////////////////////////////////////////////
            // Sort hits for the current split
            ////////////////////////////////////////////////////////////////////////////////

            //ret = cudaMalloc(&device_diagonals, max_hits * sizeof(uint64_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (1). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&device_diagonals_buf, max_hits * sizeof(uint64_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (2). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&device_hits, max_hits * sizeof(uint32_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (3). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&device_hits_buf, max_hits * sizeof(uint32_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (4). Error: %d\n", ret); exit(-1); }



            // ## POINTER SECTION 5

            start = clock();

            address_checker = 0;
            base_ptr = &data_mem[0];
            address_checker = realign_address(address_checker, 8);
            uint64_t * ptr_device_diagonals = (uint64_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint64_t), 8);

            uint64_t * ptr_device_diagonals_buf = (uint64_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint64_t), 4);

            uint32_t * ptr_device_hits = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            uint32_t * ptr_device_hits_buf = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);
            

            // We will actually sort the diagonals with associated values 0,1,2... to n and use these to index the hits array
            ret = cudaMemcpy(ptr_device_hits, ascending_numbers, n_hits_found*sizeof(uint32_t), cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Uploading forward device hits. Error: %d\n", ret); exit(-1); }

            ret = cudaMemcpy(ptr_device_diagonals, diagonals, n_hits_found*sizeof(uint64_t), cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Uploading device diagonals. Error: %d\n", ret); exit(-1); }

            cub::DoubleBuffer<uint64_t> d_diagonals(ptr_device_diagonals, ptr_device_diagonals_buf);
            cub::DoubleBuffer<uint32_t> d_hits(ptr_device_hits, ptr_device_hits_buf);

            d_temp_storage = NULL;
            temp_storage_bytes = 0;
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_diagonals, d_hits, n_hits_found);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Bad pre-sorting (3). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            // Allocate temporary storage
            ret = cudaMalloc(&d_temp_storage, temp_storage_bytes);
            if(ret != cudaSuccess){ fprintf(stderr, "Bad allocating of temp storage for hits sorting (1). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_diagonals, d_hits, n_hits_found);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "CUB sorting failed on hits. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            
            // Download hits (actually just number indices)
            ret = cudaMemcpy(indexing_numbers, ptr_device_hits_buf, n_hits_found*sizeof(uint32_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device hits. Error: %d\n", ret); exit(-1); }

            ret = cudaMemcpy(diagonals, ptr_device_diagonals_buf, n_hits_found*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device diagonals. Error: %d\n", ret); exit(-1); }

            //
            //for(i=0; i<n_hits_found; i++){
            //    fprintf(stdout, "%"PRIu32" %"PRIu32"\n", hits[indexing_numbers[i]].p1, hits[indexing_numbers[i]].p2);
            //}

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] sort hits Q-R t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            memset(filtered_hits_x, 0x0000, n_hits_found * sizeof(uint32_t));
            memset(filtered_hits_y, 0x0000, n_hits_found * sizeof(uint32_t));

            start = clock();

            uint32_t n_hits_kept = filter_hits_forward(diagonals, indexing_numbers, hits, filtered_hits_x, filtered_hits_y, n_hits_found);

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] filter hits Q-R t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            fprintf(stdout, "[INFO] Remaining hits %"PRIu32"\n", n_hits_kept);

            //for(i=0; i<n_hits_kept; i++){
                //printf("%"PRIu64"\n", diagonals[i]);
                //fprintf(out, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,32,32,32,1.0,1.0,0,0\n", filtered_hits_x[i], filtered_hits_y[i], filtered_hits_x[i]+32, filtered_hits_y[i]+32);
            //}
            
            ret = cudaFree(d_temp_storage);
            //ret = cudaFree(device_hits);
            //ret = cudaFree(device_diagonals);
            //ret = cudaFree(device_diagonals_buf);
            //ret = cudaFree(device_hits_buf);

            

            ////////////////////////////////////////////////////////////////////////////////
            // Generate FORWARD frags
            ////////////////////////////////////////////////////////////////////////////////

            // Allocate both sequences
            //ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for query sequence in device. Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&seq_dev_mem_aux, words_at_once * sizeof(char)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for ref sequence in device. Error: %d\n", ret); exit(-1); }

            //ret = cudaMalloc(&device_filt_hits_x, max_hits * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device filtered hits query. Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&device_filt_hits_y, max_hits * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device filtered hits ref. Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&left_offset, max_hits * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device offset left frags query. Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&right_offset, max_hits * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device offset right frags query. Error: %d\n", ret); exit(-1); }


            // ## POINTER SECTION 6

            start = clock();

            address_checker = 0;
            base_ptr = &data_mem[0];
            ptr_seq_dev_mem = (char *) (base_ptr);
            address_checker += words_at_once;

            ptr_seq_dev_mem_aux = (char *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once, 4);

            uint32_t * ptr_device_filt_hits_x = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            uint32_t * ptr_device_filt_hits_y = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            uint32_t * ptr_left_offset = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            uint32_t * ptr_right_offset = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);



            ret = cudaMemcpy(ptr_seq_dev_mem, &query_seq_host[pos_in_query-words_at_once], MIN(query_len - (pos_in_query - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device for frags. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_seq_dev_mem_aux, &ref_seq_host[pos_in_ref-words_at_once], MIN(ref_len - (pos_in_ref - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device for frags. Error: %d\n", ret); exit(-1); }
            
            ret = cudaMemcpy(ptr_device_filt_hits_x, filtered_hits_x, n_hits_kept * sizeof(uint32_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits x in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_device_filt_hits_y, filtered_hits_y, n_hits_kept * sizeof(uint32_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits y in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(ptr_left_offset, 0x0, n_hits_kept * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy left offset in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(ptr_right_offset, 0x0, n_hits_kept * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy right offset in device. Error: %d\n", ret); exit(-1); }

            //
            //for(i=n_hits_kept-1; i>1; i--){
            //    printf(" Frag %"PRIu64" \t x: %.32s %"PRIu64"\n", i, &query_seq_host[filtered_hits_x[i]], filtered_hits_x[i]);
            //    printf(" \t\t y: %.32s %"PRIu64"\n", &ref_seq_host[filtered_hits_y[i]], filtered_hits_y[i]);
            //}

            number_of_blocks = n_hits_kept; 
            //number_of_blocks = 20; // REMOVE !!

            if(number_of_blocks != 0)
            {
                kernel_frags_forward_register<<<number_of_blocks, threads_number>>>(ptr_device_filt_hits_x, ptr_device_filt_hits_y, ptr_left_offset, ptr_right_offset, ptr_seq_dev_mem, ptr_seq_dev_mem_aux, query_len, ref_len, pos_in_query-words_at_once, pos_in_ref-words_at_once, MIN(pos_in_query, query_len), MIN(pos_in_ref, ref_len));
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Failed on generating forward frags. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            }

            ret = cudaMemcpy(host_left_offset, ptr_left_offset, n_hits_kept * sizeof(uint32_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back left offset. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(host_right_offset, ptr_right_offset, n_hits_kept * sizeof(uint32_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back right offset. Error: %d\n", ret); exit(-1); }

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] frags Q-R t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            /*
            char name[100] = "\0";
            sprintf(name, "onlyfrags-forward_%d", split);

            FILE * anything3 = fopen(name, "wt");
            print_header(anything3, query_len, ref_len);
            for(i=0; i<n_hits_kept; i++){
                uint64_t best_xStart = filtered_hits_x[i] - host_left_offset[i];
                uint64_t best_xEnd = filtered_hits_x[i] + host_right_offset[i];
                uint64_t best_yStart = filtered_hits_y[i] - host_left_offset[i];
                uint64_t best_yEnd = filtered_hits_y[i] + host_right_offset[i];

                int64_t d = (filtered_hits_x[i] - filtered_hits_y[i]);
                //fprintf(anything3, "hitx: %"PRIu64" hity: %"PRIu64" (d: %"PRId64") -> Frag,(%"PRId64"),%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64"\n", filtered_hits_x[i], filtered_hits_y[i], d, (int64_t)best_xStart-(int64_t)best_yStart, best_xStart, best_yStart, best_xEnd, best_yEnd, best_xEnd-best_xStart);
                //fprintf(anything3, "Frag,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,%"PRIu64",32,32,1.0,1.0,0,0\n", d, filtered_hits_x[i], filtered_hits_y[i], best_xStart, best_yStart, best_xEnd, best_yEnd, best_xEnd-best_xStart);
                fprintf(anything3, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,%"PRIu64",32,32,1.0,1.0,0,0\n", best_xStart, best_yStart, best_xEnd, best_yEnd, best_xEnd-best_xStart);
            }
            fclose(anything3);
            */
            

            //cudaFree(seq_dev_mem);
            //cudaFree(seq_dev_mem_aux);
            //cudaFree(device_filt_hits_x);
            //cudaFree(device_filt_hits_y);
            //cudaFree(left_offset);
            //cudaFree(right_offset);
            
            filter_and_write_frags(filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, n_hits_kept, out, 'f', ref_len, min_length);

        }

        // Restart the reference for every block in query
        pos_in_ref = 0;

        ////////////////////////////////////////////////////////////////////////////////
        // Run the reference blocks BUT REVERSED !
        ////////////////////////////////////////////////////////////////////////////////


        while(pos_in_ref < ref_len){

            ////////////////////////////////////////////////////////////////////////////////
            // FORWARD strand in the reference BUT REVERSED !
            ////////////////////////////////////////////////////////////////////////////////

            uint32_t items_read_y = MIN(ref_len - pos_in_ref, words_at_once);

            //ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char));

            // Allocate words table
            //ret = cudaMalloc(&keys, words_at_once * sizeof(uint64_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device reversed (1). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&values, words_at_once * sizeof(uint32_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device reversed (2). Error: %d\n", ret); exit(-1); }
            
            // ## POINTER SECTION 7

            start = clock();

            address_checker = 0;
            base_ptr = &data_mem[0];
            ptr_seq_dev_mem = (char *) (base_ptr);
            address_checker = realign_address(address_checker + words_at_once, 8);

            ptr_keys = (uint64_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 4);

            ptr_values = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 4);


            // Load sequence chunk into ram
            ret = cudaMemcpy(ptr_seq_dev_mem, &ref_rev_seq_host[pos_in_ref], items_read_y, cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device reversed. Error: %d\n", ret); exit(-1); }

            // Run kmers
            ret = cudaMemset(ptr_keys, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(ptr_values, 0xFFFFFFFF, words_at_once * sizeof(uint32_t));
            //number_of_blocks = (((items_read_y - KMER_SIZE + 1)) / (threads_number*4)); 
            //kernel_register_fast_hash_rotational<<<number_of_blocks, threads_number>>>(keys, values, seq_dev_mem, pos_in_ref);
            
            number_of_blocks = ((items_read_y - KMER_SIZE + 1))/64;

            if(number_of_blocks != 0)
            {
                kernel_index_global32<<<number_of_blocks, 64>>>(ptr_keys, ptr_values, ptr_seq_dev_mem, pos_in_ref);
                
                //number_of_blocks = ((items_read_y - KMER_SIZE + 1))/threads_number;
                //kernel_index_global32<<<number_of_blocks, threads_number>>>(ptr_keys, ptr_values, ptr_seq_dev_mem, pos_in_ref);
    
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on ref reversed. Error: %d\n", ret); exit(-1); }
            }

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] words RC t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            //cudaFree(seq_dev_mem);

            ////////////////////////////////////////////////////////////////////////////////
            // Sort reference FORWARD kmers BUT REVERSED !
            ////////////////////////////////////////////////////////////////////////////////

            //ret = cudaMalloc(&keys_buf, words_at_once * sizeof(uint64_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device reversed (3). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&values_buf, words_at_once * sizeof(uint32_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device reversed (4). Error: %d\n", ret); exit(-1); }

            // ## POINTER SECTION 8

            start = clock();

            address_checker = realign_address(address_checker, 8);
            ptr_keys_buf = (uint64_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 4);

            ptr_values_buf = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 4);

            ret = cudaMemset(ptr_keys_buf, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(ptr_values_buf, 0xFFFFFFFF, words_at_once * sizeof(uint32_t));
            cub::DoubleBuffer<uint64_t> d_keys_ref(ptr_keys, ptr_keys_buf);
            cub::DoubleBuffer<uint32_t> d_values_ref(ptr_values, ptr_values_buf);
            d_temp_storage = NULL;
            temp_storage_bytes = 0;
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys_ref, d_values_ref, items_read_y);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Bad pre-sorting (2). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            // Allocate temporary storage
            ret = cudaMalloc(&d_temp_storage, temp_storage_bytes);
            if(ret != cudaSuccess){ fprintf(stderr, "Bad allocating of temp storage for words sorting (2). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys_ref, d_values_ref, items_read_y);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "CUB sorting failed on ref. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            


            // Download kmers
            ret = cudaMemcpy(dict_y_keys, ptr_keys_buf, items_read_y*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (3). Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(dict_y_values, ptr_values_buf, items_read_y*sizeof(uint32_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (4). Error: %d\n", ret); exit(-1); }

            ret = cudaFree(d_temp_storage);
            if(ret != cudaSuccess){ fprintf(stderr, "Bad free of temp storage (2): %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            pos_in_ref += words_at_once;

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] sort words RC t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            ////////////////////////////////////////////////////////////////////////////////
            // Generate hits for the current split BUT REVERSED !
            ////////////////////////////////////////////////////////////////////////////////

            
            //cudaFree(keys);
            //cudaFree(values);
            //cudaFree(keys_buf);
            //cudaFree(values_buf);

            start = clock();

            uint32_t n_hits_found;
            if(fast == 1)
                n_hits_found = generate_hits_fast(max_hits, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len);
            else
                n_hits_found = generate_hits_sensitive(max_hits, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len, max_frequency);
            
            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] hits Q-RC t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            fprintf(stdout, "[INFO] Generated %"PRIu32" hits on reversed split %d -> (%d%%)[%u,%u]{%u,%u}\n", n_hits_found, split, (int)((100*MIN((uint64_t)pos_in_ref, (uint64_t)ref_len))/(uint64_t)ref_len), pos_in_query, pos_in_ref, items_read_x, items_read_y);


            //printf("VALHALA 1\n");
            //continue;

            ////////////////////////////////////////////////////////////////////////////////
            // Sort hits for the current split BUT REVERSED !
            ////////////////////////////////////////////////////////////////////////////////

            //ret = cudaMalloc(&device_diagonals, max_hits * sizeof(uint64_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (1). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&device_diagonals_buf, max_hits * sizeof(uint64_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (2). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&device_hits, max_hits * sizeof(uint32_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (3). Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&device_hits_buf, max_hits * sizeof(uint32_t));
            //if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (4). Error: %d\n", ret); exit(-1); }

            // ## POINTER SECTION 9

            start = clock();

            base_ptr = &data_mem[0];
            address_checker = 0;
            address_checker = realign_address(address_checker, 8);
            uint64_t * ptr_device_diagonals = (uint64_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint64_t), 8);

            uint64_t * ptr_device_diagonals_buf = (uint64_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint64_t), 4);

            uint32_t * ptr_device_hits = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            uint32_t * ptr_device_hits_buf = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            // We will actually sort the diagonals with associated values 0,1,2... to n and use these to index the hits array
            ret = cudaMemcpy(ptr_device_hits, ascending_numbers, n_hits_found*sizeof(uint32_t), cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Uploading device reverse hits. Error: %d\n", ret); exit(-1); }

            ret = cudaMemcpy(ptr_device_diagonals, diagonals, n_hits_found*sizeof(uint64_t), cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Uploading device diagonals. Error: %d\n", ret); exit(-1); }

            cub::DoubleBuffer<uint64_t> d_diagonals(ptr_device_diagonals, ptr_device_diagonals_buf);
            cub::DoubleBuffer<uint32_t> d_hits(ptr_device_hits, ptr_device_hits_buf);
            d_temp_storage = NULL;
            temp_storage_bytes = 0;
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_diagonals, d_hits, n_hits_found);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Bad pre-sorting (3). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            // Allocate temporary storage
            ret = cudaMalloc(&d_temp_storage, temp_storage_bytes);
            if(ret != cudaSuccess){ fprintf(stderr, "Bad allocating of temp storage for hits sorting (1). Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            
            cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_diagonals, d_hits, n_hits_found);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "CUB sorting failed on hits. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            
            // Download hits (actually just number indices)
            ret = cudaMemcpy(indexing_numbers, ptr_device_hits_buf, n_hits_found*sizeof(uint32_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device hits. Error: %d\n", ret); exit(-1); }

            ret = cudaMemcpy(diagonals, ptr_device_diagonals_buf, n_hits_found*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device diagonals. Error: %d\n", ret); exit(-1); }

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] sorting hits Q-RC t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif


            memset(filtered_hits_x, 0x0000, n_hits_found * sizeof(uint32_t));
            memset(filtered_hits_y, 0x0000, n_hits_found * sizeof(uint32_t));

            start = clock();
            uint32_t n_hits_kept = filter_hits_reverse(diagonals, indexing_numbers, hits, filtered_hits_x, filtered_hits_y, n_hits_found);

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] filter hits Q-RC t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            fprintf(stdout, "[INFO] Remaining hits %"PRIu32"\n", n_hits_kept);

            //printf("VALHALA 2\n");
            //continue;
            

            // Filtered hits are in order to diagonal (x-y)*l + x (PROVED)
            // Notice::::: These are sorted equally as with forward hits because the input sequence
            // Y has already been fully reversed and complemented :) 
            /*
            FILE * anything2 = fopen("onlyhits-reverse.csv", "wt");
            print_header(anything2, query_len, ref_len);
            for(i=0; i<n_hits_kept; i++){
                int64_t d = (filtered_hits_x[i] - filtered_hits_y[i]);
                uint64_t best_yStart = ref_len - filtered_hits_y[i] - 1;
                uint64_t best_yEnd = ref_len - (filtered_hits_y[i]+32) - 1;
                //int64_t d = filtered_hits_x[i] + best_yStart;
                fprintf(anything2, "Frag,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,32,32,32,1.0,1.0,0,0\n", d, filtered_hits_x[i], best_yStart, filtered_hits_x[i]+32, best_yEnd);
            }
            fclose(anything2);
            */

            
            
            
            ret = cudaFree(d_temp_storage);
            //ret = cudaFree(device_hits);
            //ret = cudaFree(device_diagonals);
            //ret = cudaFree(device_diagonals_buf);
            //ret = cudaFree(device_hits_buf);

            ////////////////////////////////////////////////////////////////////////////////
            // Generate REVERSE frags
            ////////////////////////////////////////////////////////////////////////////////

            // Allocate both sequences
            //ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for query sequence in device. Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&seq_dev_mem_aux, words_at_once * sizeof(char)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for ref sequence in device. Error: %d\n", ret); exit(-1); }

            //ret = cudaMalloc(&device_filt_hits_x, max_hits * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device filtered hits query. Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&device_filt_hits_y, max_hits * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device filtered hits ref. Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&left_offset, max_hits * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device offset left frags query. Error: %d\n", ret); exit(-1); }
            //ret = cudaMalloc(&right_offset, max_hits * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device offset right frags query. Error: %d\n", ret); exit(-1); }


            // ## POINTER SECTION 10

            start = clock();

            address_checker = 0;
            base_ptr = &data_mem[0];
            ptr_seq_dev_mem = (char *) (base_ptr);
            address_checker += words_at_once;

            ptr_seq_dev_mem_aux = (char *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once, 4);

            uint32_t * ptr_device_filt_hits_x = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            uint32_t * ptr_device_filt_hits_y = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            uint32_t * ptr_left_offset = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);

            uint32_t * ptr_right_offset = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 4);


            ret = cudaMemcpy(ptr_seq_dev_mem, &query_seq_host[pos_in_query-words_at_once], MIN(query_len - (pos_in_query - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device for frags. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_seq_dev_mem_aux, &ref_rev_seq_host[pos_in_ref-words_at_once], MIN(ref_len - (pos_in_ref - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device for frags. Error: %d\n", ret); exit(-1); }
            
            ret = cudaMemcpy(ptr_device_filt_hits_x, filtered_hits_x, n_hits_kept * sizeof(uint32_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits x in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_device_filt_hits_y, filtered_hits_y, n_hits_kept * sizeof(uint32_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits y in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(ptr_left_offset, 0x0, n_hits_kept * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy left offset in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(ptr_right_offset, 0x0, n_hits_kept * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy right offset in device. Error: %d\n", ret); exit(-1); }

            //
            //for(i=n_hits_kept-1; i>1; i--){
            //    printf(" Frag %"PRIu64" \t x: %.32s %"PRIu64"\n", i, &query_seq_host[filtered_hits_x[i]], filtered_hits_x[i]);
            //    printf(" \t\t y: %.32s %"PRIu64"\n", &ref_seq_host[filtered_hits_y[i]], filtered_hits_y[i]);
            //}

            //printf("VALHALA 3\n");
            //continue;


            number_of_blocks = n_hits_kept; 
            //number_of_blocks = 100; 
            //printf("sending blocks: %u\n", number_of_blocks);
            //printf("We are sending: posinquery-wo=%u posinref-wo=%u MIN1=%u MIN2=%u\n", pos_in_query-words_at_once, pos_in_ref-words_at_once, MIN(pos_in_query, query_len), MIN(pos_in_ref, ref_len));

            //number_of_blocks = 20; // REMOVE !!

            if(number_of_blocks != 0)
            {
                kernel_frags_reverse_register<<<number_of_blocks, threads_number>>>(ptr_device_filt_hits_x, ptr_device_filt_hits_y, ptr_left_offset, ptr_right_offset, ptr_seq_dev_mem, ptr_seq_dev_mem_aux, query_len, ref_len, pos_in_query-words_at_once, pos_in_ref-words_at_once, MIN(pos_in_query, query_len), MIN(pos_in_ref, ref_len));

                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Failed on generating forward frags. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            }
                


            ret = cudaMemcpy(host_left_offset, ptr_left_offset, n_hits_kept * sizeof(uint32_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back left offset. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(host_right_offset, ptr_right_offset, n_hits_kept * sizeof(uint32_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back right offset. Error: %d\n", ret); exit(-1); }

            end = clock();
#ifdef SHOWTIME
            fprintf(stdout, "[INFO] frags Q-RC t=%f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

            /*           
            FILE * anything = fopen("onlyfrags-reverse.csv", "wt");
            print_header(anything, query_len, ref_len);
            for(i=0; i<n_hits_kept; i++){
                uint64_t best_xStart = filtered_hits_x[i] - host_left_offset[i];
                uint64_t best_xEnd = filtered_hits_x[i] + host_right_offset[i];
                uint64_t best_yStart = filtered_hits_y[i] - host_left_offset[i];
                uint64_t best_yEnd = filtered_hits_y[i] + host_right_offset[i];
                int64_t d = (filtered_hits_x[i] + filtered_hits_y[i]);


                fprintf(anything, "hitx: %"PRIu64" hity: %"PRIu64" (d: %"PRId64") -> Frag,(%"PRId64"),%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64"\n", filtered_hits_x[i], filtered_hits_y[i], d, (int64_t)best_xStart+(int64_t)best_yStart, best_xStart, best_yStart, best_xEnd, best_yEnd, best_xEnd-best_xStart);
                //fprintf(anything, "Frag,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,%"PRIu64",32,32,1.0,1.0,0,0\n", d, filtered_hits_x[i], filtered_hits_y[i], best_xStart, best_yStart, best_xEnd, best_yEnd, best_xEnd-best_xStart);
            }
            fclose(anything);
            */
            

            //cudaFree(seq_dev_mem);
            //cudaFree(seq_dev_mem_aux);
            //cudaFree(device_filt_hits_x);
            //cudaFree(device_filt_hits_y);
            //cudaFree(left_offset);
            //cudaFree(right_offset);

            
            //printf("VALHALA 4\n");
            //continue;

            filter_and_write_frags(filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, n_hits_kept, out, 'r', ref_len, min_length);


        }

        pos_in_ref = 0;


        ++split;
    }

    fclose(out);
    
    
    
    
    

    fprintf(stdout, "[INFO] Completed\n");

    fclose(query);
    fclose(ref);

    /*
    free(query_seq_host);
    free(ref_seq_host);
    free(ref_rev_seq_host);
    free(dict_x_keys); 
    free(dict_x_values); 
    free(dict_y_keys); 
    free(dict_y_values); 
    free(hits); 
    free(filtered_hits_x); 
    free(filtered_hits_y); 
    free(host_left_offset); 
    free(host_right_offset); 
    free(diagonals); 
    free(ascending_numbers); 
    free(indexing_numbers); 
    */

    cudaFree(data_mem);
    cudaFreeHost(host_pinned_mem);
    /*
    cudaFreeHost(query_seq_host);
    cudaFreeHost(ref_seq_host);
    cudaFreeHost(ref_rev_seq_host);
    cudaFreeHost(dict_x_keys); 
    cudaFreeHost(dict_x_values); 
    cudaFreeHost(dict_y_keys); 
    cudaFreeHost(dict_y_values); 
    cudaFreeHost(hits); 
    cudaFreeHost(filtered_hits_x); 
    cudaFreeHost(filtered_hits_y); 
    cudaFreeHost(host_left_offset); 
    cudaFreeHost(host_right_offset); 
    cudaFreeHost(diagonals); 
    cudaFreeHost(ascending_numbers); 
    cudaFreeHost(indexing_numbers); 
    */

    return 0;
}

void print_header(FILE * out, uint32_t query_len, uint32_t ref_len){

    fprintf(out, "All by-Identity Ungapped Fragments (Hits based approach)\n");
    fprintf(out, "[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>\n");
    fprintf(out, "SeqX filename : undef\n");
    fprintf(out, "SeqY filename : undef\n");
    fprintf(out, "SeqX name : undef\n");
    fprintf(out, "SeqY name : undef\n");
    fprintf(out, "SeqX length : %"PRIu32"\n", query_len);
    fprintf(out, "SeqY length : %"PRIu32"\n", ref_len);
    fprintf(out, "Min.fragment.length : undef\n");
    fprintf(out, "Min.Identity : undef\n");
    fprintf(out, "Tot Hits (seeds) : undef\n");
    fprintf(out, "Tot Hits (seeds) used: undef\n");
    fprintf(out, "Total fragments : undef\n");
    fprintf(out, "========================================================\n");
    fprintf(out, "Total CSB: 0\n");
    fprintf(out, "========================================================\n");
    fprintf(out, "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");
}


uint64_t memory_allocation_chooser(uint64_t total_memory)
{
   

    if(total_memory <= 4340179200) return 100*1024*1024;
    else if(total_memory <= 6442450944) return 150*1024*1024;
    else if(total_memory <= 8689934592) return 200*1024*1024;
    return 300*1024*1024;
 
}





