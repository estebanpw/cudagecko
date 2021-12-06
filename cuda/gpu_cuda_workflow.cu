
// Standard utilities and common systems includes
#include "kernels.cuh"
#include "cpu_functions.c"
#include "cub/cub.cuh"
#include <moderngpu/kernel_mergesort.hxx>
#include <cuda_profiler_api.h>

#define BILLION 1000 * 1000 * 1000;


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////


void print_header(FILE * out, uint32_t query_len, uint32_t ref_len);
float factor_chooser(uint64_t total_memory);
uint64_t memory_allocation_chooser(uint64_t total_memory);
char * dump_memory_region(char * ptr_pointer, uint64_t size);

int main(int argc, char ** argv)
{
    clock_t start = clock(), end = clock();
#ifdef SHOWTIME
    struct timespec HD_start, HD_end;
    uint64_t time_seconds = 0, time_nanoseconds = 0;
#endif
    uint32_t i, min_length = 64, max_frequency = 0, n_frags_per_block = 32;
    float factor = -1;
    int fast = 0; // sensitive is default
    // Parameters for the generation of hits in the gpu
    // Do not alter these 
    uint64_t _u64_SPLITHITS = 1;
    uint64_t global_device_RAM = 0;
    float _f_SECTIONS = 0.2;
    unsigned selected_device = 0;
    FILE * query = NULL, * ref = NULL, * out = NULL;
    init_args(argc, argv, &query, &selected_device, &ref, &out, &min_length, &fast, &max_frequency, &factor, &n_frags_per_block, &_u64_SPLITHITS, &_f_SECTIONS, &global_device_RAM);

    if(fast == 3)
        fprintf(stdout, "[INFO] Using AVX512 intrinsics to compute vector hits.\n");

    ////////////////////////////////////////////////////////////////////////////////
    // Get info of devices
    ////////////////////////////////////////////////////////////////////////////////

    int ret_num_devices;
    int ret;
    
    //int sharedMemPerBlock, multiProcessorCount;
    
    // Query how many devices there are
    if(cudaSuccess != (ret = cudaGetDeviceCount(&ret_num_devices))){ fprintf(stderr, "Failed to query number of devices\n"); exit(-1); }

    cudaDeviceProp device;

    for(i=0; i<ret_num_devices; i++){
        if( cudaSuccess != (ret = cudaGetDeviceProperties(&device, i))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }
        fprintf(stdout, "\tDevice [%" PRIu32"]: %s\n", i, device.name);
        uint64_t ram = device.totalGlobalMem;
        fprintf(stdout, "\t\tGlobal mem   : %" PRIu64" (%" PRIu64" MB)\n", ram, ram / (1024*1024));
    }

    if( cudaSuccess != (ret = cudaSetDevice(selected_device))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }
    fprintf(stdout, "[INFO] Using device %d\n", selected_device);

    if( cudaSuccess != (ret = cudaGetDeviceProperties(&device, selected_device))){ fprintf(stderr, "Failed to get cuda device property: %d\n", ret); exit(-1); }
    
    // If no amount of max memory was specified by the user:
    if(global_device_RAM == 0)
        global_device_RAM = device.totalGlobalMem;

    // Select blocks
    uint32_t insider_kernel_blocks = (uint32_t) device.multiProcessorCount * (uint32_t) (device.sharedMemPerBlock / 768); // 768 b is how much shared memory is used by hits kernel
    
    end = clock();
#ifdef SHOWTIME
    fprintf(stdout, "[INFO] INIT 1 t= %f\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

    start = clock();

    // Calculate how much ram we can use for every chunk
    uint64_t effective_global_ram =  (global_device_RAM - memory_allocation_chooser(global_device_RAM)); //Minus 100 to 300 MBs for other stuff

    // Retune factor
    if(factor < 0) // Only if left by default
        factor = factor_chooser(global_device_RAM);

    // We will do the one-time alloc here
    // i.e. allocate a pool once and used it manually
    char * data_mem;

    // One workload depends on number of words (words + sortwords + generate hits)
    // The other one depends on number of hits (sort hits + filterhits + frags)
    if(fast == 2) factor = 0.45;
    else if(fast == 1) factor = 0.45;
    uint64_t bytes_for_words = (factor * effective_global_ram); // 512 MB for words
    uint64_t words_at_once = bytes_for_words / (8+8+4+4); 
    // We have to subtract the bytes for words as well as the region for storing the DNA sequence
    uint64_t max_hits = (effective_global_ram - bytes_for_words - words_at_once) / (2*8); // Max hits must fit twice because of the sorting
    uint64_t bytes_to_subtract = max_hits * 8;

    // Allocate the pool
    ret = cudaMalloc(&data_mem, (effective_global_ram) * sizeof(char)); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pool memory in device. Error: %d\n", ret); exit(-1); }
    fprintf(stdout, "[INFO] Memory pool at %p of size %lu bytes\n", data_mem, (effective_global_ram) * sizeof(char));

    // Reserve a slice for the pre allocated pool for moderngpu (and for other things as well while it is not used)
    char * pre_alloc = &data_mem[effective_global_ram - bytes_to_subtract]; // points to the last section of the big pool

    fprintf(stdout, "[INFO] You can have %" PRIu64" MB for words (i.e. %" PRIu64" words), and %" PRIu64" MB for hits (i.e. %" PRIu64" hits)\n", 
        bytes_for_words / (1024*1024), words_at_once, (effective_global_ram - bytes_for_words - words_at_once) / (1024*1024), max_hits);

    fprintf(stdout, "[INFO] Filtering at a minimum length of %" PRIu32" bps\n", min_length);
    if(fast == 1) 
        fprintf(stdout, "[INFO] Running on fast mode (some repetitive seeds will be skipped)\n");
    else if(fast == 2)
        fprintf(stdout, "[INFO] Running on hyper fast mode (some repetitive seeds will be skipped)\n");
    else if(fast == 3)
        fprintf(stdout, "[INFO] Running on sensitive mode (all hits computed in CPU, rest in GPU)\n");
    else
        fprintf(stdout, "[INFO] Running on sensitive mode (full execution on GPU[mf:%" PRIu32"])\n", max_frequency);

    // Assign memory pool to moderngpu
    Mem_pool mptr;
    mptr.mem_ptr = pre_alloc;
    mptr.address = 0;
    mptr.limit = bytes_to_subtract;
    mgpu::standard_context_t context(false, 0, &mptr);
    
    // Set working sizes (these will change throughout execution)
    size_t threads_number = 32;
    size_t number_of_blocks;

#ifdef SHOWTIME
    clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

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

#ifdef SHOWTIME
    clock_gettime(CLOCK_MONOTONIC, &HD_end);
    time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
    time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
#endif
    
    fprintf(stdout, "[INFO] Qlen: %" PRIu32"; Rlen: %" PRIu32"\n", query_len, ref_len);


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

    // How about one big alloc (save ~3 seconds on mallocs and improves transfer times)
    char * host_pinned_mem, * base_ptr_pinned;

    // We need for pinned memory:
    // seqx             => query_len
    // seqy             => ref_len
    // seqyrev          => ref_len
    // dict_x_keys      => words_at_once * sizeof(uint64_t)
    // dict_x_values    => words_at_once * sizeof(uint32_t)
    // dict_y_keys      => words_at_once * sizeof(uint64_t) [ONLY FOR CPU PROCESSING]
    // dict_y_values    => words_at_once * sizeof(uint32_t) [ONLY FOR CPU PROCESSING]
    // filtered_hits_x  => max_hits * sizeof(uint32_t)
    // filtered_hits_y  => max_hits * sizeof(uint32_t)
    // host_left        => max_hits * sizeof(uint32_t)
    // host_right       => max_hits * sizeof(uint32_t)
    // diagonals        => max_hits * sizeof(uint64_t)

    uint64_t pinned_bytes_on_host = query_len + 2*ref_len + words_at_once*4 + words_at_once*8 + max_hits*8 + max_hits*16;
    pinned_bytes_on_host += 1024*1024; // Adding 1 MB for the extra padding in the realignments
    if(fast != 0) pinned_bytes_on_host += words_at_once*8 + words_at_once*4;
    uint64_t pinned_address_checker = 0;

    fprintf(stdout, "[INFO] Allocating on host %" PRIu64" bytes (i.e. %" PRIu64" MBs)\n", pinned_bytes_on_host, pinned_bytes_on_host / (1024*1024));
    ret = cudaHostAlloc(&host_pinned_mem, pinned_bytes_on_host, cudaHostAllocMapped); 
    if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate pinned memory for pool. Error: %d\n", ret); exit(-1); }
    
#ifdef SHOWTIME
    clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif
    char * query_seq_host, * ref_seq_host, * ref_rev_seq_host;
    base_ptr_pinned = (char *) &host_pinned_mem[0];
    query_seq_host = (char *) &host_pinned_mem[0];
    pinned_address_checker = realign_address(pinned_address_checker + query_len, 4);

    ref_seq_host = (char *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + ref_len, 4);

    ref_rev_seq_host = (char *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + ref_len, 256);


    if(query_seq_host == NULL || ref_seq_host == NULL || ref_rev_seq_host == NULL) terror("Could not allocate memory for sequences in host");

    ////////////////////////////////////////////////////////////////////////////////
    // Read sequences and reverse the reference
    ////////////////////////////////////////////////////////////////////////////////

    // Create streams to allow concurrent copy and execute
    uint32_t n_streams = 4;
    cudaStream_t streams[n_streams];

    for(i=0; i<n_streams; i++) cudaStreamCreate(&streams[i]);

    fprintf(stdout, "[INFO] Loading query\n");
    memcpy(query_seq_host, pro_q_buffer, query_len);
    fprintf(stdout, "[INFO] Loading reference\n");

    memcpy(ref_seq_host, pro_r_buffer, ref_len);
    fprintf(stdout, "[INFO] Reversing reference\n");

    free(pro_q_buffer);
    free(pro_r_buffer);

    // ## POINTER SECTION 0 [These depict sections of code where the base pointer to the memory pool is changed/used]
    uint64_t address_checker = 0;
    
    char * ptr_seq_dev_mem_aux = &data_mem[0];
    address_checker = realign_address(address_checker + ref_len, 128);
    char * ptr_seq_dev_mem_reverse_aux = &data_mem[address_checker];
    address_checker = realign_address(address_checker + ref_len, 4);

    char * ptr_reverse_write[n_streams];

    uint32_t chunk_size = ref_len/n_streams + 4; // Always force one divisor more cause of decimal loss
    if(chunk_size % 128 != 0){ chunk_size += 128 - (chunk_size % 128); } // Since sequence starts at 0, making the chunks multiple of 128 guarantees 100% GL efficiency


    ptr_reverse_write[0] = ptr_seq_dev_mem_reverse_aux;
    for(i=1; i<n_streams-1; i++){
        ptr_reverse_write[i] = ptr_reverse_write[i-1] + chunk_size;
    }
    ptr_reverse_write[n_streams-1] = ptr_reverse_write[n_streams-2] + chunk_size;
    uintptr_t where_ends = (uintptr_t) (ptr_reverse_write[n_streams-1] + MIN(chunk_size, ref_len - chunk_size*(n_streams-1)));
    if(where_ends % 128 != 0) ptr_reverse_write[n_streams-1] += 128 - (where_ends % 128);

    threads_number = 128;
    if(ref_len > 1024){

        number_of_blocks = chunk_size/threads_number + threads_number; 
        if(number_of_blocks % threads_number != 0){ number_of_blocks += 1; }
        uint32_t offset;

        
        for(i=0; i<n_streams; i++)
        {
            offset = chunk_size * i;
            ret = cudaMemcpyAsync(ptr_seq_dev_mem_aux + offset, ref_seq_host + offset, MIN(chunk_size, ref_len - offset), cudaMemcpyHostToDevice, streams[i]);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy reference sequence to device for reversing. Error: %d\n", ret); exit(-1); }
            kernel_reverse_complement<<<number_of_blocks, threads_number, 0, streams[i]>>>(ptr_seq_dev_mem_aux + offset, ptr_reverse_write[i], MIN(chunk_size, ref_len - offset));

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
        kernel_reverse_complement<<<number_of_blocks, threads_number>>>(ptr_seq_dev_mem_aux, ptr_seq_dev_mem_reverse_aux, ref_len);
        
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "Could not compute reverse on reference. Error: %d\n", ret); exit(-1); }

        ret = cudaMemcpy(ref_rev_seq_host, ptr_seq_dev_mem_reverse_aux, ref_len, cudaMemcpyDeviceToHost);
        if(ret != cudaSuccess){ fprintf(stderr, "Could not copy reference sequence to device for reversing. Error: %d\n", ret); exit(-1); }
    }
    ret = cudaDeviceSynchronize();


#ifdef SHOWTIME
    end = clock();
    clock_gettime(CLOCK_MONOTONIC, &HD_end);
    time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
    time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
    time_seconds *= BILLION;
    fprintf(stdout, "[INFO] rev comp t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
    time_seconds = 0;
    time_nanoseconds = 0;
#endif 

#ifdef SHOWTIME
    clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

    fprintf(stdout, "[INFO] Showing start of reference sequence:\n");
    
    fprintf(stdout, "\t(Begin ref)%.64s\n", ref_seq_host);
    fprintf(stdout, "\t(Begin rev)%.64s\n", ref_rev_seq_host);
    fprintf(stdout, "\t(End   ref)%.64s\n", &ref_seq_host[ref_len-64]);
    fprintf(stdout, "\t(End   rev)%.64s\n", &ref_rev_seq_host[ref_len-64]);

    // Write header to CSV
    print_header(out, query_len, ref_len);

    ////////////////////////////////////////////////////////////////////////////////
    // Allocation of pointers
    ////////////////////////////////////////////////////////////////////////////////


    // Allocate memory in host to download kmers and store hits
    // These depend on the number of words
    uint64_t * dict_x_keys, * dict_y_keys; // Keys are hashes (64b), values are positions (32b)
    uint32_t * dict_x_values, * dict_y_values;

    pinned_address_checker = realign_address(pinned_address_checker, 8);
    dict_x_keys = (uint64_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + words_at_once * sizeof(uint64_t), 4);

    dict_x_values = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + words_at_once * sizeof(uint32_t), 8);

    if(fast != 0){
        dict_y_keys = (uint64_t *) (base_ptr_pinned + pinned_address_checker);
        pinned_address_checker = realign_address(pinned_address_checker + words_at_once * sizeof(uint64_t), 4);

        dict_y_values = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
        pinned_address_checker = realign_address(pinned_address_checker + words_at_once * sizeof(uint32_t), 4);
    }

    // These are now depending on the number of hits
    uint32_t * filtered_hits_x, * filtered_hits_y;

    /*
    // TODO remove
    Hit * hits;
    pinned_address_checker = realign_address(pinned_address_checker, 8);
    hits = (Hit *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(Hit), 4);
    */

    filtered_hits_x = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 4);

    filtered_hits_y = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 4);

    uint32_t * host_left_offset,  * host_right_offset;
    uint64_t * diagonals;

    pinned_address_checker = realign_address(pinned_address_checker, 4);
    host_left_offset = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 4);

    host_right_offset = (uint32_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint32_t), 8);

    diagonals = (uint64_t *) (base_ptr_pinned + pinned_address_checker);
    pinned_address_checker = realign_address(pinned_address_checker + max_hits * sizeof(uint64_t), 4);

    ////////////////////////////////////////////////////////////////////////////////
    // Read the query and reference in blocks
    ////////////////////////////////////////////////////////////////////////////////

#ifdef SHOWTIME
    clock_gettime(CLOCK_MONOTONIC, &HD_end);
    time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
    time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
    time_seconds *= BILLION;
    fprintf(stdout, "[INFO] INIT 3 t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
    time_seconds = 0;
    time_nanoseconds = 0;
#endif

    int split = 0;
    uint32_t pos_in_query = 0, pos_in_ref = 0;

    while(pos_in_query < query_len){

#ifdef SHOWTIME
        clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif
        // ## POINTER SECTION 1
        address_checker = 0;
        char * ptr_seq_dev_mem = &data_mem[0];
        char * base_ptr = ptr_seq_dev_mem;
        address_checker = realign_address(address_checker + words_at_once, 256);
        
        uint64_t * ptr_keys = (uint64_t *) (base_ptr + address_checker); // We have to realign because of the arbitrary length of the sequence chars
        address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 128);

        uint32_t * ptr_values = (uint32_t *) (base_ptr + address_checker); 
        address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 256);
        uint64_t address_CHECKPOINT = address_checker;
        

        fprintf(stdout, "[EXECUTING] Running split %d -> (%d%%)[%u,%u]\n", split, (int)((100*(uint64_t)pos_in_query)/(uint64_t)query_len), pos_in_query, pos_in_ref);

        uint32_t items_read_x = MIN(query_len - pos_in_query, words_at_once);


        ////////////////////////////////////////////////////////////////////////////////
        // Run kmers for query
        ////////////////////////////////////////////////////////////////////////////////
        
        // Load sequence chunk into ram 
        ret = cudaMemcpy(ptr_seq_dev_mem, &query_seq_host[pos_in_query], items_read_x, cudaMemcpyHostToDevice);
        if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device. Error: %d\n", ret); exit(-1); }

        // Initialize space
        ret = cudaMemset(ptr_keys, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
        ret = cudaMemset(ptr_values, 0xFFFFFFFF, words_at_once * sizeof(uint32_t));
        ret = cudaDeviceSynchronize();
        
        number_of_blocks = (items_read_x - KMER_SIZE + 1)/(64) + 1;
        if(number_of_blocks != 0)
        {
            //cudaProfilerStart();
            kernel_index_global32<<<number_of_blocks, 64>>>(ptr_keys, ptr_values, ptr_seq_dev_mem, pos_in_query, items_read_x);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on query. Error: %d\n", ret); exit(-1); }
        }
        else
        {
            fprintf(stdout, "[WARNING] Zero blocks for query words\n");
        }


#ifdef SHOWTIME
        clock_gettime(CLOCK_MONOTONIC, &HD_end);
        time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
        time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
        time_seconds *= BILLION;
        fprintf(stdout, "[INFO] words Q t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
        time_seconds = 0;
        time_nanoseconds = 0;
#endif

        ////////////////////////////////////////////////////////////////////////////////
        // Sort the query kmers
        ////////////////////////////////////////////////////////////////////////////////

        // Notice: Now that we are usiung pooled memory
        // I have not "freed" the part corresponding to the sequence (ptr_seq_dev_mem)
        // And thus next points build upon that
        // But thats no problem because it is a small fraction of memory

        // ## POINTER SECTION 2
#ifdef SHOWTIME
        clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

        mergesort(ptr_keys, ptr_values, items_read_x, mgpu::less_t<uint64_t>(), context);
        ret = cudaDeviceSynchronize();
        if(ret != cudaSuccess){ fprintf(stderr, "MERGESORT sorting failed on query. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
        
        // Download sorted kmers [ They will be reuploaded afterwards ]
        ret = cudaMemcpy(dict_x_keys, ptr_keys, items_read_x*sizeof(uint64_t), cudaMemcpyDeviceToHost);
        if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (1). Error: %d\n", ret); exit(-1); }
        ret = cudaMemcpy(dict_x_values, ptr_values, items_read_x*sizeof(uint32_t), cudaMemcpyDeviceToHost);
        if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers (2). Error: %d\n", ret); exit(-1); }

        
#ifdef SHOWTIME
        clock_gettime(CLOCK_MONOTONIC, &HD_end);
        time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
        time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );

        time_seconds *= BILLION;
        fprintf(stdout, "[INFO] sortwords Q t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
        time_seconds = 0;
        time_nanoseconds = 0;
#endif 

        pos_in_query += words_at_once;

        ////////////////////////////////////////////////////////////////////////////////
        // Run the reference blocks
        ////////////////////////////////////////////////////////////////////////////////

        // These definitions are for the processing of hits - reused in reference and query
        // TODO put these in corresponding pinned parts
        uint64_t * ptr_device_diagonals;
        int32_t * ptr_device_error;
        uint32_t * ptr_hits_log, * ptr_hits_log_extra;
        uint64_t * ptr_keys_2;
        uint32_t * ptr_values_2;

        while(pos_in_ref < ref_len){

            ////////////////////////////////////////////////////////////////////////////////
            // FORWARD strand in the reference
            ////////////////////////////////////////////////////////////////////////////////
#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif
            uint32_t items_read_y = MIN(ref_len - pos_in_ref, words_at_once);

            // ## POINTER SECTION 3
            ptr_seq_dev_mem = &data_mem[0];
            address_checker = address_CHECKPOINT; // The checkpoint is here since pointers are changed for the hits and frags execution

            // Recopy x words because they were overwritten in the hits generation
            ret = cudaMemcpy(ptr_keys, dict_x_keys, items_read_x*sizeof(uint64_t), cudaMemcpyHostToDevice);
            ret = cudaMemcpy(ptr_values, dict_x_values, items_read_x*sizeof(uint32_t), cudaMemcpyHostToDevice);
            // Memset everything that comes after the uploaded keys in case it was overwritten
            
            if(items_read_x < words_at_once){
                ret = cudaMemset(&ptr_keys[items_read_x], 0xFFFFFFFF, (words_at_once-items_read_x) * sizeof(uint64_t));
                ret = cudaMemset(&ptr_values[items_read_x], 0xFFFFFFFF, (words_at_once-items_read_x) * sizeof(uint32_t));
            }
            
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy words back in posterior iterations. Error: %d\n", ret); exit(-1); }

            ptr_keys_2 = (uint64_t *) (base_ptr + address_checker); // We have to realign because of the arbitrary length of the sequence chars
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 256);

            ptr_values_2 = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 256);

            // Load sequence chunk into ram
            ret = cudaMemcpy(ptr_seq_dev_mem, &ref_seq_host[pos_in_ref], items_read_y, cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device. Error: %d\n", ret); exit(-1); }

            // Run kmers
            ret = cudaMemset(ptr_keys_2, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(ptr_values_2, 0xFFFFFFFF, words_at_once * sizeof(uint32_t));
            ret = cudaDeviceSynchronize();
            
            number_of_blocks = ((items_read_y - KMER_SIZE + 1))/(64) + 1;
            if(number_of_blocks != 0)
            {

                kernel_index_global32<<<number_of_blocks, 64>>>(ptr_keys_2, ptr_values_2, ptr_seq_dev_mem, pos_in_ref, items_read_y);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on ref. Error: %d\n", ret); exit(-1); }

            }
            else
            {
                fprintf(stdout, "[WARNING] Zero blocks for ref words\n");
            }
#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] words R t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 
            

            ////////////////////////////////////////////////////////////////////////////////
            // Sort reference FORWARD kmers
            ////////////////////////////////////////////////////////////////////////////////


            // ## POINTER SECTION 4
#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

            mergesort(ptr_keys_2, ptr_values_2, items_read_y, mgpu::less_t<uint64_t>(), context);

            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "MODERNGPU sorting failed on ref words. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            pos_in_ref += words_at_once;

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] sortwords R t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 

            ////////////////////////////////////////////////////////////////////////////////
            // Generate FORWARD hits for the current split
            ////////////////////////////////////////////////////////////////////////////////

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif


            uint32_t n_hits_found = 0;

            if(fast == 3){
                ret = cudaMemcpy(dict_y_keys, ptr_keys_2, items_read_y*sizeof(uint64_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for vectorized hit generation on forward (1). Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(dict_y_values, ptr_values_2, items_read_y*sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for vectorized hit generation on forward (2). Error: %d\n", ret); exit(-1); }
                n_hits_found = generate_hits_sensitive_avx512(max_hits, diagonals, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len);
            }
            else if(fast == 2){
                ret = cudaMemcpy(dict_y_keys, ptr_keys_2, items_read_y*sizeof(uint64_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for fast hit generation on forward (1). Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(dict_y_values, ptr_values_2, items_read_y*sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for fast hit generation on forward (2). Error: %d\n", ret); exit(-1); }
                n_hits_found = generate_hits_fast(max_hits, diagonals, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len);     
            }
            else if(fast == 1){
                ret = cudaMemcpy(dict_y_keys, ptr_keys_2, items_read_y*sizeof(uint64_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for fast hit generation on forward (1). Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(dict_y_values, ptr_values_2, items_read_y*sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for fast hit generation on forward (2). Error: %d\n", ret); exit(-1); }
                n_hits_found = generate_hits_sensitive(max_hits, diagonals, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len, max_frequency, fast);
            }
//#ifdef AVX512CUSTOM
//                n_hits_found = generate_hits_sensitive_avx512(max_hits, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len);
//#else
//                n_hits_found = generate_hits_sensitive(max_hits, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len, max_frequency, fast);
//#endif
            else{


                ////////////////////////////////////////////////////////////////////////////////
                // GPU hits pipeline
                ////////////////////////////////////////////////////////////////////////////////


                uint32_t n_blocks_hits = insider_kernel_blocks; //insider_kernel_blocks; //

                ptr_device_error = (int32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(int32_t), 4);

                ptr_hits_log = (uint32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(uint32_t) * n_blocks_hits, 4);
                ret = cudaMemset(ptr_hits_log, 0x00000000, sizeof(uint32_t) * n_blocks_hits);

                // Set error status to one
                ret = cudaMemset(ptr_device_error, 0x00000000, sizeof(int32_t));
                ret = cudaDeviceSynchronize();
                

                uint64_t hits_in_first_mem_block = max_hits / _u64_SPLITHITS; 
                uint64_t hits_in_second_mem_block = (2 * max_hits - hits_in_first_mem_block) - 1000*1000; // Some has to be removed due to all allocated variables on pool (besides words and seq) TODO: allocate them at the beginning
                uint64_t mem_block = (hits_in_first_mem_block)/n_blocks_hits;
                if(mem_block % 8 != 0) mem_block -= mem_block % 8; // Each section should be a multiple of 8
                uint64_t max_extra_sections = n_blocks_hits * _f_SECTIONS;
                uint64_t extra_large_mem_block = (hits_in_second_mem_block)/max_extra_sections;
                if(extra_large_mem_block % 8 != 0) extra_large_mem_block -= extra_large_mem_block % 8; // Each section should be a multiple of 8

                ptr_hits_log_extra = (uint32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(uint32_t) * max_extra_sections, 4);
                ret = cudaMemset(ptr_hits_log_extra, 0x00000000, sizeof(uint32_t) * max_extra_sections);

                uint32_t * ptr_leftmost_key_x = (uint32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(uint32_t), 4);
                uint32_t * ptr_leftmost_key_y = (uint32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(uint32_t), 4);
                int32_t * ptr_atomic_distributer = (int32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(int32_t), 4);
                ret = cudaMemset(ptr_atomic_distributer, 0x00000000, sizeof(int32_t));
                int32_t * ptr_queue = (int32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(int32_t), 4);
                ret = cudaMemset(ptr_queue, 0x00000000, sizeof(int32_t)); // 0x..2000 is 8192
                ret = cudaDeviceSynchronize();

                kernel_find_leftmost_items<<<1, 1>>>(ptr_keys, ptr_leftmost_key_x, ptr_keys_2, ptr_leftmost_key_y, items_read_x, items_read_y);
                ret = cudaDeviceSynchronize();
                uint32_t leftmost_key_x, leftmost_key_y;
                if(ret != cudaSuccess){ fprintf(stderr, "Error searching true leftmost elements on device on forward strand. Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(&leftmost_key_x, ptr_leftmost_key_x, sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading leftmost element X on forward strand. Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(&leftmost_key_y, ptr_leftmost_key_y, sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading leftmost element Y on forward strand. Error: %d\n", ret); exit(-1); }

                if(leftmost_key_x >= items_read_x || leftmost_key_y >= items_read_y){ fprintf(stderr, "Bad binary search of leftmost items on forward strand. %u, %u out of [%u, %u]\n",  leftmost_key_x, leftmost_key_y, items_read_x, items_read_y); exit(-1); }


                address_checker = realign_address(address_checker, 256);
                ptr_device_diagonals = (uint64_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + hits_in_first_mem_block * sizeof(uint64_t), 256);

                uint64_t * ptr_auxiliary_hit_memory = (uint64_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + hits_in_second_mem_block * sizeof(uint64_t), 8);

                ret = cudaMemset(ptr_device_diagonals, 0xFFFFFFFF, sizeof(uint64_t)*hits_in_first_mem_block);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Setting to 0xFF..FF the device diagonals on forward strand. Error: %d\n", ret); exit(-1); }
                
                ret = cudaMemset(ptr_auxiliary_hit_memory, 0xFFFFFFFF, sizeof(uint64_t)*hits_in_second_mem_block);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Setting to 0xFF..FF the extra device diagonals on forward strand. Error: %d\n", ret); exit(-1); }

                ////////////////////////////////////////////////////////////////////////////////
                // Hits-generation with load balancing (increasingly smaller sets of words)
                ////////////////////////////////////////////////////////////////////////////////

                //cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024000*10);
                
                kernel_hits_load_balancing<<<n_blocks_hits, 32>>>(ptr_keys, ptr_keys_2, ptr_values, ptr_values_2, ptr_device_diagonals, (int32_t) mem_block, 
                    leftmost_key_x, leftmost_key_y, ptr_device_error, ref_len, ptr_hits_log, ptr_atomic_distributer, ptr_auxiliary_hit_memory,
                    (uint32_t) extra_large_mem_block, (uint32_t) max_extra_sections, ptr_hits_log_extra, ptr_queue);//, ptr_messages_log);

                ret = cudaDeviceSynchronize();

                if(ret != cudaSuccess){ fprintf(stdout, "Fatal error generating hits on device on forward strand (Tip: try reducing the factor parameter). Error: %d\n", ret); fflush(stdout); fprintf(stderr, "Fatal error generating hits on device on forward strand (Tip: try reducing the factor parameter). Error: %d\n", ret); exit(-1); }
                int32_t device_error;
                ret = cudaMemcpy(&device_error, ptr_device_error, sizeof(int32_t), cudaMemcpyDeviceToHost);
                int32_t reached_sections = -1;
                ret = cudaMemcpy(&reached_sections, ptr_atomic_distributer, sizeof(int32_t), cudaMemcpyDeviceToHost);

                if(ret != cudaSuccess){ fprintf(stderr, "Downloading error status on hits generation on forward strand. Error: %d\n", ret); exit(-1); }
                if(device_error < 0) { fprintf(stderr, "Error generating hits on device on forward strand (Tip: try reducing the factor parameter). Error: %d\n", device_error); exit(-1); }

                ////////////////////////////////////////////////////////////////////////////////
                // Hits compacting
                ////////////////////////////////////////////////////////////////////////////////

                uint32_t * hits_log = (uint32_t *) malloc(n_blocks_hits*sizeof(uint32_t)); 
                uint32_t * extra_log = (uint32_t *) malloc(max_extra_sections*sizeof(uint32_t)); 
                uint32_t * accum_log = (uint32_t *) malloc(n_blocks_hits*sizeof(uint32_t)); 
                
                ret = cudaMemcpy(hits_log, ptr_hits_log, sizeof(uint32_t)*n_blocks_hits, cudaMemcpyDeviceToHost);
                ret = cudaMemcpy(extra_log, ptr_hits_log_extra, sizeof(uint32_t)*max_extra_sections, cudaMemcpyDeviceToHost);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading hits log on forward strand. Error: %d\n", ret); exit(-1); }
                for(i=0; i<n_blocks_hits; i++) {
                    accum_log[i] = n_hits_found;
                    n_hits_found += hits_log[i];
                }

                // First step: measure how many hits can be stored at once (in the worst case) in the words section
                // This is (consecutive region): ptr_keys,ptr_values,ptr_keys_2,ptr_values_2
                // And amounts for: words_at_once * (8+4+8+4) bytes
                // which equals max number of 8-byte diagonals: 3*words_at_once

                uint64_t * ptr_copy_place_diagonals = (uint64_t *) &ptr_keys[0];
                uint32_t max_copy_diagonals = (2 * words_at_once * sizeof(uint64_t) + 2 * words_at_once * sizeof(uint32_t)) / sizeof(uint64_t) ;
                uint32_t runs = (uint32_t) hits_in_first_mem_block / max_copy_diagonals + 1;

                // Upload accumulated (overwrite sequence data since its no longer needed)
                uint32_t * ptr_accum_log = (uint32_t *) (&data_mem[0]);
                ret = cudaMemcpy(ptr_accum_log, accum_log, sizeof(uint32_t)*n_blocks_hits, cudaMemcpyHostToDevice);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Uploading accum hits log on forward strand. Error: %d\n", ret); exit(-1); }

                if(n_hits_found > 0){
                    i = 1;
                    uint32_t real_i = 0;
                    uint32_t offset_i = 0;
                    while(i <= n_blocks_hits){
                        
                        
                        
                        if(i == n_blocks_hits || accum_log[i] - offset_i > max_copy_diagonals){

                            kernel_compact_hits<<<i-real_i, 512>>>(&ptr_device_diagonals[real_i], &ptr_hits_log[real_i],
                                &ptr_accum_log[real_i], mem_block, ptr_copy_place_diagonals, offset_i);

                            ret = cudaDeviceSynchronize();

                            if(ret != cudaSuccess){ fprintf(stderr, "Could not compact hits on first stage (%u). Error: %d\n", i, ret); exit(-1); }
                            ret = cudaMemcpy(ptr_device_diagonals, ptr_copy_place_diagonals, sizeof(uint64_t)*(accum_log[i-1] - offset_i + hits_log[i-1]), cudaMemcpyDeviceToDevice);
                            ret = cudaDeviceSynchronize();

                            if(i == n_blocks_hits) break;
                            real_i = i;
                            offset_i = accum_log[i];
                            
                        }else{
                            ++i;
                        }
                        
                    }
                }
                

                // Second run
                uint32_t second_hits_found = 0;

                for(i=0; i<max_extra_sections; i++) {
                    accum_log[i] = second_hits_found;
                    second_hits_found += extra_log[i];
                }

                ret = cudaMemcpy(ptr_accum_log, accum_log, sizeof(uint32_t)*max_extra_sections, cudaMemcpyHostToDevice);
                ret = cudaDeviceSynchronize();

                if(reached_sections > 0){
                    uint32_t sections_per_run = max_copy_diagonals / extra_large_mem_block;
                    runs = (uint32_t) reached_sections / sections_per_run + 1;
                
                    uint32_t total_offset = 0;
                    uint32_t total_copied = 0;
                    for(i=0; i<runs; i++){
                        uint32_t j;
                        uint32_t copied = 0;
                        for(j=i*sections_per_run; j<(i+1)*sections_per_run && j<reached_sections; j++){
                            ret = cudaMemcpy(&ptr_copy_place_diagonals[accum_log[j] - total_offset], &ptr_auxiliary_hit_memory[j * extra_large_mem_block], sizeof(uint64_t)*extra_log[j], cudaMemcpyDeviceToDevice);
                            copied += extra_log[j];
                        }
                        ret = cudaDeviceSynchronize();
                        if(ret != cudaSuccess){ fprintf(stderr, "Error copying to safe place on %u. Error: %d\n", i, ret); exit(-1); }
                        total_offset = accum_log[j]; //accum_log[j * sections_per_run];
                        ret = cudaMemcpy(&ptr_device_diagonals[total_copied + n_hits_found], &ptr_copy_place_diagonals[0], sizeof(uint64_t)*copied, cudaMemcpyDeviceToDevice);
                        ret = cudaDeviceSynchronize();
                        total_copied += copied;
                        if(ret != cudaSuccess){ fprintf(stderr, "Error batching safe place on %u. Error: %d\n", i, ret); exit(-1); }
                    }
                    
                }

                // Add together total number of hits
                n_hits_found += second_hits_found;

                // TODO make these pinned to avoid frees
                free(hits_log);
                free(extra_log);
                free(accum_log);

#ifdef SHOWTIME
                fprintf(stdout, "[INFO] Max reached sections on forward = %d out of %" PRIu64"\n", reached_sections, max_extra_sections);
#endif

            }

            if(n_hits_found > max_hits){ fprintf(stderr, "Found too many hits on forward strand (%u). Use a lower factor.\n", n_hits_found); exit(-1); }

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] hits Q-R t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
            fprintf(stdout, "[INFO] Generated %" PRIu32" hits on split %d -> (%d%%)[%u,%u]{%u,%u}\n", n_hits_found, split, (int)((100*MIN((uint64_t)pos_in_ref, (uint64_t)ref_len))/(uint64_t)ref_len), pos_in_query, pos_in_ref, items_read_x, items_read_y);
#endif 

            ////////////////////////////////////////////////////////////////////////////////
            // Sort hits for the current split
            ////////////////////////////////////////////////////////////////////////////////

            
#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif
            // ## POINTER SECTION 5
            // No longer required

            if(fast != 0){
                ptr_device_diagonals = (uint64_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + max_hits * sizeof(uint64_t), 256);
                ret = cudaMemcpy(ptr_device_diagonals, diagonals, n_hits_found*sizeof(uint64_t), cudaMemcpyHostToDevice);
                if(ret != cudaSuccess){ fprintf(stderr, "Could not copy fast hits to device on forward strand. Error: %d\n", ret); exit(-1); }
            }

            mergesort(ptr_device_diagonals, n_hits_found, mgpu::less_t<uint64_t>(), context);
            

            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "MODERNGPU sorting failed on query-ref hits. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] sorthits Q-R t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 

            memset(filtered_hits_x, 0x0000, n_hits_found * sizeof(uint32_t));
            memset(filtered_hits_y, 0x0000, n_hits_found * sizeof(uint32_t));
#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

            //cudaProfilerStart();
            //kernel_filter_hits_parallel_shared<<<n_hits_found/(64)+1, 64>>>(ptr_device_diagonals, ref_len, n_hits_found);
            kernel_filter_hits_parallel<<<n_hits_found/(64)+1, 64>>>(ptr_device_diagonals, ref_len, n_hits_found);
            ret = cudaDeviceSynchronize();
            //cudaProfilerStop();
            //return 0;


            if(ret != cudaSuccess){ fprintf(stderr, "FILTER HITS failed on query-ref hits. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            ret = cudaMemcpy(diagonals, ptr_device_diagonals, n_hits_found*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device diagonals. Error: %d\n", ret); exit(-1); }

            

            // We need to download the filtered hits at some point anyway because they contain the positions of the seeds
            uint32_t n_hits_kept = filter_hits_cpu(diagonals, filtered_hits_x, filtered_hits_y, n_hits_found);
            

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] filterhits Q-R t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
            fprintf(stdout, "[INFO] Remaining hits %" PRIu32"\n", n_hits_kept);
            
#endif

            ////////////////////////////////////////////////////////////////////////////////
            // Generate FORWARD frags
            ////////////////////////////////////////////////////////////////////////////////

            
#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif
            // ## POINTER SECTION 6
            address_checker = 0;
            base_ptr = &data_mem[0];
            ptr_seq_dev_mem = (char *) (base_ptr);
            address_checker = realign_address(address_checker + words_at_once, 4);

            ptr_seq_dev_mem_aux = (char *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once, 128);

            uint32_t * ptr_device_filt_hits_x = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 128);

            uint32_t * ptr_device_filt_hits_y = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 128);
            

            // For half of the hits and frags we use the memory pool and for the other half we use the prealloc pool for sorting
            // Otherwise if there are too many hits they wont fit in just one pool
            uint64_t address_checker_pre_alloc = 0;
            char * base_pre_alloc_ptr = &pre_alloc[0];

            uint32_t * ptr_left_offset = (uint32_t *) (base_pre_alloc_ptr + address_checker_pre_alloc);
            address_checker_pre_alloc = realign_address(address_checker_pre_alloc + max_hits * sizeof(uint32_t), 128);

            uint32_t * ptr_right_offset = (uint32_t *) (base_pre_alloc_ptr + address_checker_pre_alloc);
            address_checker_pre_alloc = realign_address(address_checker_pre_alloc + max_hits * sizeof(uint32_t), 128);

            // And copy all required stuff
            ret = cudaMemcpy(ptr_seq_dev_mem, &query_seq_host[pos_in_query-words_at_once], MIN(query_len - (pos_in_query - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device for frags. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_seq_dev_mem_aux, &ref_seq_host[pos_in_ref-words_at_once], MIN(ref_len - (pos_in_ref - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device for frags. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_device_filt_hits_x, filtered_hits_x, n_hits_kept * sizeof(uint32_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits x in device on forward strand. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_device_filt_hits_y, filtered_hits_y, n_hits_kept * sizeof(uint32_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits y in device on forward strand. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(ptr_left_offset, 0x0, n_hits_kept * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy left offset in device on forward strand. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(ptr_right_offset, 0x0, n_hits_kept * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy right offset in device on forward strand. Error: %d\n", ret); exit(-1); }
            ret = cudaDeviceSynchronize();

            number_of_blocks = (n_hits_kept / n_frags_per_block) + 1;

            if(number_of_blocks != 0)
            {
                kernel_frags_forward_register<<<number_of_blocks, 32>>>(ptr_device_filt_hits_x, ptr_device_filt_hits_y, ptr_left_offset, 
                    ptr_right_offset, ptr_seq_dev_mem, ptr_seq_dev_mem_aux, query_len, ref_len, pos_in_query-words_at_once, pos_in_ref-words_at_once, 
                    MIN(pos_in_query, query_len), MIN(pos_in_ref, ref_len), n_hits_kept, n_frags_per_block);
                
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Failed on generating forward frags. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            }

            ret = cudaMemcpy(host_left_offset, ptr_left_offset, n_hits_kept * sizeof(uint32_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back left offset. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(host_right_offset, ptr_right_offset, n_hits_kept * sizeof(uint32_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back right offset. Error: %d\n", ret); exit(-1); }
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Error copying offsets back on forward frags. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] frags Q-R t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

            filter_and_write_frags(filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, n_hits_kept, out, 'f', ref_len, min_length);

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] filterFrags Q-R t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 

        } 

        ////////////////////////////////////////////////////////////////////////////////
        // This concludes the execution for the forward strand
        // Now we perform the same operations but in reverse
        ////////////////////////////////////////////////////////////////////////////////

        

        // Restart the reference for every block in query
        pos_in_ref = 0;

        ////////////////////////////////////////////////////////////////////////////////
        // Run the reference blocks BUT REVERSED !
        ////////////////////////////////////////////////////////////////////////////////

        while(pos_in_ref < ref_len){

            ////////////////////////////////////////////////////////////////////////////////
            // FORWARD strand in the reference BUT REVERSED !
            ////////////////////////////////////////////////////////////////////////////////

            // In case some was overwritten during hits generation
            ret = cudaMemset(base_ptr, 0xFFFFFFFF, (words_at_once-items_read_x) * sizeof(char));

            uint32_t items_read_y = MIN(ref_len - pos_in_ref, words_at_once);

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif
            // ## POINTER SECTION 7
            // We will give it the spot of the reference kmers/keys
            base_ptr = &data_mem[0];
            ptr_seq_dev_mem = (char *) (base_ptr);
            address_checker = address_CHECKPOINT;

            ptr_keys_2 = (uint64_t *) (base_ptr + address_checker); // We have to realign because of the arbitrary length of the sequence chars
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint64_t), 128);

            ptr_values_2 = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once * sizeof(uint32_t), 128);

            // Load sequence chunk into ram
            ret = cudaMemcpy(ptr_seq_dev_mem, &ref_rev_seq_host[pos_in_ref], items_read_y, cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device reversed. Error: %d\n", ret); exit(-1); }

            // Recopy because it gets overwritten in the generation of hits
            ret = cudaMemcpy(ptr_keys, dict_x_keys, items_read_x*sizeof(uint64_t), cudaMemcpyHostToDevice);
            ret = cudaMemcpy(ptr_values, dict_x_values, items_read_x*sizeof(uint32_t), cudaMemcpyHostToDevice);
            
            // Memset everything that comes after the uploaded keys in case it was overwritten
            if(items_read_x < words_at_once){
                ret = cudaMemset(&ptr_keys[items_read_x], 0xFFFFFFFF, (words_at_once-items_read_x) * sizeof(uint64_t));
                ret = cudaMemset(&ptr_values[items_read_x], 0xFFFFFFFF, (words_at_once-items_read_x) * sizeof(uint32_t));
            }
            
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy words back in posterior iterations for the reverse. Error: %d\n", ret); exit(-1); }

            ret = cudaMemset(ptr_keys_2, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(ptr_values_2, 0xFFFFFFFF, words_at_once * sizeof(uint32_t));
            ret = cudaDeviceSynchronize();
            
            number_of_blocks = ((items_read_y - KMER_SIZE + 1))/64 + 1;

            if(number_of_blocks != 0)
            {
                kernel_index_global32<<<number_of_blocks, 64>>>(ptr_keys_2, ptr_values_2, ptr_seq_dev_mem, pos_in_ref, items_read_y);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Could not compute kmers on ref reversed. Error: %d\n", ret); exit(-1); }
            }
#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );

            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] words RC t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 

            ////////////////////////////////////////////////////////////////////////////////
            // Sort reference FORWARD kmers BUT REVERSED !
            ////////////////////////////////////////////////////////////////////////////////

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

            // ## POINTER SECTION 8
            // Not required anymore

            mergesort(ptr_keys_2, ptr_values_2, items_read_y, mgpu::less_t<uint64_t>(), context);
            ret = cudaDeviceSynchronize();

            if(ret != cudaSuccess){ fprintf(stderr, "MODERNGPU sorting failed on words reverse. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            // Increment position 
            pos_in_ref += words_at_once;

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] sortwords RC t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 

            ////////////////////////////////////////////////////////////////////////////////
            // Generate hits for the current split BUT REVERSED !
            ////////////////////////////////////////////////////////////////////////////////

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

            uint32_t n_hits_found = 0;

            if(fast == 3){
                ret = cudaMemcpy(dict_y_keys, ptr_keys_2, items_read_y*sizeof(uint64_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for vectorized hit generation on reverse (1). Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(dict_y_values, ptr_values_2, items_read_y*sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for vectorized hit generation on reverse (2). Error: %d\n", ret); exit(-1); }
                n_hits_found = generate_hits_sensitive_avx512(max_hits, diagonals, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len);
            }
            else if(fast == 2){
                ret = cudaMemcpy(dict_y_keys, ptr_keys_2, items_read_y*sizeof(uint64_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for fast hit generation on reverse (1). Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(dict_y_values, ptr_values_2, items_read_y*sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for fast hit generation on reverse (2). Error: %d\n", ret); exit(-1); }
                n_hits_found = generate_hits_fast(max_hits, diagonals, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len);
            }else if(fast == 1){
                ret = cudaMemcpy(dict_y_keys, ptr_keys_2, items_read_y*sizeof(uint64_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for fast hit generation on reverse (1). Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(dict_y_values, ptr_values_2, items_read_y*sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading device kmers for fast hit generation on reverse (2). Error: %d\n", ret); exit(-1); }
                n_hits_found = generate_hits_sensitive(max_hits, diagonals, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len, max_frequency, fast);
            }
//            else
//#ifdef AVX512CUSTOM
//                n_hits_found = generate_hits_sensitive_avx512(max_hits, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len);
//#else
//                n_hits_found = generate_hits_sensitive(max_hits, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, items_read_x, items_read_y, query_len, ref_len, max_frequency, fast);
//#endif
            else{
                uint32_t n_blocks_hits = insider_kernel_blocks;

                ptr_device_error = (int32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(int32_t), 4);

                ptr_hits_log = (uint32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(uint32_t) * n_blocks_hits, 4);
                ret = cudaMemset(ptr_hits_log, 0x00000000, sizeof(uint32_t) * n_blocks_hits);

                // Set error status to zero (negative returns mean error)
                ret = cudaMemset(ptr_device_error, 0x00000000, sizeof(int32_t));
                ret = cudaDeviceSynchronize();

                uint64_t hits_in_first_mem_block = max_hits / _u64_SPLITHITS;
                uint64_t hits_in_second_mem_block = (2 * max_hits - hits_in_first_mem_block) - 1000*1000; // Some has to be removed due to all allocated variables on pool (besides words and seq) TODO: allocate them at the beginning
                uint64_t mem_block = (hits_in_first_mem_block)/n_blocks_hits;
                uint64_t max_extra_sections = n_blocks_hits * _f_SECTIONS;
                uint64_t extra_large_mem_block = (hits_in_second_mem_block)/max_extra_sections;

                ptr_hits_log_extra = (uint32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(uint32_t) * max_extra_sections, 4);
                ret = cudaMemset(ptr_hits_log_extra, 0x00000000, sizeof(uint32_t) * max_extra_sections);
                

                uint32_t * ptr_leftmost_key_x = (uint32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(uint32_t), 4);
                uint32_t * ptr_leftmost_key_y = (uint32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(uint32_t), 4);
                int32_t * ptr_atomic_distributer = (int32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(int32_t), 4);
                ret = cudaMemset(ptr_atomic_distributer, 0x00000000, sizeof(int32_t));
                int32_t * ptr_queue = (int32_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + sizeof(int32_t), 4);
                ret = cudaMemset(ptr_queue, 0x00000000, sizeof(int32_t)); // 0x..2000 is 8192
                ret = cudaDeviceSynchronize();


                kernel_find_leftmost_items<<<1, 1>>>(ptr_keys, ptr_leftmost_key_x, ptr_keys_2, ptr_leftmost_key_y, items_read_x, items_read_y);
                ret = cudaDeviceSynchronize();
                uint32_t leftmost_key_x, leftmost_key_y;
                if(ret != cudaSuccess){ fprintf(stderr, "Error searching true leftmost elements on device on the reverse. Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(&leftmost_key_x, ptr_leftmost_key_x, sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading leftmost element X on the reverse. Error: %d\n", ret); exit(-1); }
                ret = cudaMemcpy(&leftmost_key_y, ptr_leftmost_key_y, sizeof(uint32_t), cudaMemcpyDeviceToHost);
                if(ret != cudaSuccess){ fprintf(stderr, "Downloading leftmost element Y on the reverse. Error: %d\n", ret); exit(-1); }

                if(leftmost_key_x >= items_read_x || leftmost_key_y >= items_read_y){ fprintf(stderr, "Bad binary search of leftmost items on reverse strand. %u, %u out of [%u, %u]\n",  leftmost_key_x, leftmost_key_y, items_read_x, items_read_y); exit(-1); }

                address_checker = realign_address(address_checker, 256);
                ptr_device_diagonals = (uint64_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + hits_in_first_mem_block * sizeof(uint64_t), 256);

                uint64_t * ptr_auxiliary_hit_memory = (uint64_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + hits_in_second_mem_block * sizeof(uint64_t), 128);

                ret = cudaMemset(ptr_device_diagonals, 0xFFFFFFFF, sizeof(uint64_t)*hits_in_first_mem_block);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Setting to 0xFF..FF the device diagonals on the reverse. Error: %d\n", ret); exit(-1); }
                
                ret = cudaMemset(ptr_auxiliary_hit_memory, 0xFFFFFFFF, sizeof(uint64_t)*hits_in_second_mem_block);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Setting to 0xFF..FF the extra device diagonals on the reverse. Error: %d\n", ret); exit(-1); }
                
                ////////////////////////////////////////////////////////////////////////////////
                // Hits-generation with load balancing (increasingly smaller sets of words)
                ////////////////////////////////////////////////////////////////////////////////


                kernel_hits_load_balancing<<<n_blocks_hits, 32>>>(ptr_keys, ptr_keys_2, ptr_values, ptr_values_2, ptr_device_diagonals, (int32_t) mem_block, 
                    leftmost_key_x, leftmost_key_y, ptr_device_error, ref_len, ptr_hits_log, ptr_atomic_distributer, ptr_auxiliary_hit_memory,
                    (uint32_t) extra_large_mem_block, (uint32_t) max_extra_sections, ptr_hits_log_extra, ptr_queue);//, ptr_messages_log);

                ret = cudaDeviceSynchronize();

                if(ret != cudaSuccess){ fprintf(stdout, "Fatal error generating hits on device on reverse strand (Tip: try reducing the factor parameter). Error: %d\n", ret); fflush(stdout); fprintf(stderr, "Fatal error generating hits on device on reverse (Tip: try reducing the factor parameter). Error: %d\n", ret); exit(-1); }
                int32_t device_error;
                ret = cudaMemcpy(&device_error, ptr_device_error, sizeof(int32_t), cudaMemcpyDeviceToHost);
                int32_t reached_sections = -1;
                ret = cudaMemcpy(&reached_sections, ptr_atomic_distributer, sizeof(int32_t), cudaMemcpyDeviceToHost);

                if(ret != cudaSuccess){ fprintf(stderr, "Downloading error status on hits generation on reverse. Error: %d\n", ret); exit(-1); }
                if(device_error < 0) { fprintf(stderr, "Error generating hits on device on reverse (Tip: try reducing the factor parameter). Error: %d\n", device_error); exit(-1); }

                ////////////////////////////////////////////////////////////////////////////////
                // Hits compacting but REVERSED!
                ////////////////////////////////////////////////////////////////////////////////

                uint32_t * hits_log = (uint32_t *) malloc(n_blocks_hits*sizeof(uint32_t)); 
                uint32_t * extra_log = (uint32_t *) malloc(max_extra_sections*sizeof(uint32_t)); 
                uint32_t * accum_log = (uint32_t *) malloc(n_blocks_hits*sizeof(uint32_t)); 
                
                ret = cudaMemcpy(hits_log, ptr_hits_log, sizeof(uint32_t)*n_blocks_hits, cudaMemcpyDeviceToHost);
                ret = cudaMemcpy(extra_log, ptr_hits_log_extra, sizeof(uint32_t)*max_extra_sections, cudaMemcpyDeviceToHost);
                for(i=0; i<n_blocks_hits; i++) {
                    accum_log[i] = n_hits_found;
                    n_hits_found += hits_log[i];
                }

                // First step: measure how many hits can be stored at once (in the worst case) in the words section
                // This is (consecutive region): ptr_keys,ptr_values,ptr_keys_2,ptr_values_2
                // And amounts for: words_at_once * (8+4+8+4) bytes
                // which equals max number of 8-byte diagonals: 3*words_at_once

                uint64_t * ptr_copy_place_diagonals = (uint64_t *) &ptr_keys[0];
                uint32_t max_copy_diagonals = (2 * words_at_once * sizeof(uint64_t) + 2 * words_at_once * sizeof(uint32_t)) / sizeof(uint64_t) ;
                uint32_t runs = (uint32_t) hits_in_first_mem_block / max_copy_diagonals + 1;


                // Upload accumulated (overwrite sequence data since its no longer needed)
                uint32_t * ptr_accum_log = (uint32_t *) (&data_mem[0]);
                ret = cudaMemcpy(ptr_accum_log, accum_log, sizeof(uint32_t)*n_blocks_hits, cudaMemcpyHostToDevice);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Uploading accum hits log on reverse strand. Error: %d\n", ret); exit(-1); }

                if(n_hits_found > 0){
                    i = 1;
                    uint32_t real_i = 0;
                    uint32_t offset_i = 0;
                    while(i <= n_blocks_hits){
                        
                        
                        
                        if(i == n_blocks_hits || accum_log[i] - offset_i > max_copy_diagonals){

                            kernel_compact_hits<<<i-real_i, 512>>>(&ptr_device_diagonals[real_i], &ptr_hits_log[real_i],
                                &ptr_accum_log[real_i], mem_block, ptr_copy_place_diagonals, offset_i);

                            ret = cudaDeviceSynchronize();

                            if(ret != cudaSuccess){ fprintf(stderr, "Could not compact hits on first stage (%u) on reverse. Error: %d\n", i, ret); exit(-1); }
                            ret = cudaMemcpy(ptr_device_diagonals, ptr_copy_place_diagonals, sizeof(uint64_t)*(accum_log[i-1] - offset_i + hits_log[i-1]), cudaMemcpyDeviceToDevice);
                            ret = cudaDeviceSynchronize();

                            if(i == n_blocks_hits) break;
                            real_i = i;
                            offset_i = accum_log[i];
                            
                        }else{
                            ++i;
                        }
                        
                    }
                }

                // Second run
                uint32_t second_hits_found = 0;

                for(i=0; i<max_extra_sections; i++) {
                    accum_log[i] = second_hits_found;
                    second_hits_found += extra_log[i];
                }

                ret = cudaMemcpy(ptr_accum_log, accum_log, sizeof(uint32_t)*max_extra_sections, cudaMemcpyHostToDevice);
                ret = cudaDeviceSynchronize();

                if(reached_sections > 0){
                    uint32_t sections_per_run = max_copy_diagonals / extra_large_mem_block;
                    runs = (uint32_t) reached_sections / sections_per_run + 1;

                    uint32_t total_offset = 0;
                    uint32_t total_copied = 0;
                    for(i=0; i<runs; i++){
                        uint32_t j;
                        uint32_t copied = 0;
                        for(j=i*sections_per_run; j<(i+1)*sections_per_run && j<reached_sections; j++){
                            ret = cudaMemcpy(&ptr_copy_place_diagonals[accum_log[j] - total_offset], &ptr_auxiliary_hit_memory[j * extra_large_mem_block], sizeof(uint64_t)*extra_log[j], cudaMemcpyDeviceToDevice);
                            copied += extra_log[j];
                        }
                        ret = cudaDeviceSynchronize();
                        if(ret != cudaSuccess){ fprintf(stderr, "Error copying to safe place on %u on reverse. Error: %d\n", i, ret); exit(-1); }
                        total_offset = accum_log[j]; //accum_log[j * sections_per_run];
                        ret = cudaMemcpy(&ptr_device_diagonals[total_copied + n_hits_found], &ptr_copy_place_diagonals[0], sizeof(uint64_t)*copied, cudaMemcpyDeviceToDevice);
                        ret = cudaDeviceSynchronize();
                        total_copied += copied;
                        if(ret != cudaSuccess){ fprintf(stderr, "Error batching safe place on %u on reverse. Error: %d\n", i, ret); exit(-1); }
                    }
                    
                }

                // Add together total number of hits
                n_hits_found += second_hits_found;


                // TODO add this to pinned
                free(hits_log);
                free(extra_log);
                free(accum_log);
                
#ifdef SHOWTIME
                fprintf(stdout, "[INFO] Max reached sections on reverse = %d out of %" PRIu64"\n", reached_sections, max_extra_sections);
#endif

            }
            
            if(n_hits_found > max_hits){ fprintf(stderr, "Found too many hits on reverse strand (%u). Use a lower factor.\n", n_hits_found); exit(-1); }


#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] hits Q-RC t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
            fprintf(stdout, "[INFO] Generated %" PRIu32" hits on reversed split %d -> (%d%%)[%u,%u]{%u,%u}\n", n_hits_found, split, (int)((100*MIN((uint64_t)pos_in_ref, (uint64_t)ref_len))/(uint64_t)ref_len), pos_in_query, pos_in_ref, items_read_x, items_read_y);
#endif 
            
            ////////////////////////////////////////////////////////////////////////////////
            // Sort hits for the current split BUT REVERSED !
            ////////////////////////////////////////////////////////////////////////////////

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

            // ## POINTER SECTION 9
            //No longer required

            if(fast != 0){
                ptr_device_diagonals = (uint64_t *) (base_ptr + address_checker);
                address_checker = realign_address(address_checker + max_hits * sizeof(uint64_t), 128);
                ret = cudaMemcpy(ptr_device_diagonals, diagonals, n_hits_found*sizeof(uint64_t), cudaMemcpyHostToDevice);
                if(ret != cudaSuccess){ fprintf(stderr, "Could not copy fast hits to device on reverse strand. Error: %d\n", ret); exit(-1); }
            }

            mergesort(ptr_device_diagonals, n_hits_found, mgpu::less_t<uint64_t>(), context);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "MODERNGPU sorting failed on hits rev. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );

            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] sortinghits Q-RC t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 

            memset(filtered_hits_x, 0x0000, n_hits_found * sizeof(uint32_t));
            memset(filtered_hits_y, 0x0000, n_hits_found * sizeof(uint32_t));
#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif


            kernel_filter_hits_parallel<<<n_hits_found/(64)+1, 64>>>(ptr_device_diagonals, ref_len, n_hits_found);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "FILTER HITS failed on query-ref-comp hits on reverse. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            ret = cudaMemcpy(diagonals, ptr_device_diagonals, n_hits_found*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device diagonals on reverse. Error: %d\n", ret); exit(-1); }

            uint32_t n_hits_kept = filter_hits_cpu(diagonals, filtered_hits_x, filtered_hits_y, n_hits_found);


#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );

            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] filterhits Q-RC t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
            fprintf(stdout, "[INFO] Remaining hits %" PRIu32"\n", n_hits_kept);
#endif 


            ////////////////////////////////////////////////////////////////////////////////
            // Generate REVERSE frags
            ////////////////////////////////////////////////////////////////////////////////

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif
            // ## POINTER SECTION 10
            address_checker = 0;
            base_ptr = &data_mem[0];
            ptr_seq_dev_mem = (char *) (base_ptr);
            address_checker += words_at_once;

            ptr_seq_dev_mem_aux = (char *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + words_at_once, 128);

            uint32_t * ptr_device_filt_hits_x = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 128);

            uint32_t * ptr_device_filt_hits_y = (uint32_t *) (base_ptr + address_checker);
            address_checker = realign_address(address_checker + max_hits * sizeof(uint32_t), 128);

            // For half of the hits and frags we use the memory pool and for the other half we use the prealloc pool for sorting
            // Otherwise if there are too many hits they wont fit in just one pool
            uint64_t address_checker_pre_alloc = 0;
            char * base_pre_alloc_ptr = &pre_alloc[0];

            uint32_t * ptr_left_offset = (uint32_t *) (base_pre_alloc_ptr + address_checker_pre_alloc);
            address_checker_pre_alloc = realign_address(address_checker_pre_alloc + max_hits * sizeof(uint32_t), 128);

            uint32_t * ptr_right_offset = (uint32_t *) (base_pre_alloc_ptr + address_checker_pre_alloc);
            address_checker_pre_alloc = realign_address(address_checker_pre_alloc + max_hits * sizeof(uint32_t), 128);

            ret = cudaMemcpy(ptr_seq_dev_mem, &query_seq_host[pos_in_query-words_at_once], MIN(query_len - (pos_in_query - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device for frags. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_seq_dev_mem_aux, &ref_rev_seq_host[pos_in_ref-words_at_once], MIN(ref_len - (pos_in_ref - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device for frags. Error: %d\n", ret); exit(-1); }
            
            ret = cudaMemcpy(ptr_device_filt_hits_x, filtered_hits_x, n_hits_kept * sizeof(uint32_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits x in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(ptr_device_filt_hits_y, filtered_hits_y, n_hits_kept * sizeof(uint32_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits y in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(ptr_left_offset, 0x0, n_hits_kept * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy left offset in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(ptr_right_offset, 0x0, n_hits_kept * sizeof(uint32_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy right offset in device. Error: %d\n", ret); exit(-1); }
            ret = cudaDeviceSynchronize();

            number_of_blocks = (n_hits_kept / n_frags_per_block) + 1;

            if(number_of_blocks != 0)
            {
                // Plot twist: its the same kernel for forward and reverse since sequence is completely reversed
                kernel_frags_forward_register<<<number_of_blocks, 32>>>(ptr_device_filt_hits_x, ptr_device_filt_hits_y, ptr_left_offset, ptr_right_offset, ptr_seq_dev_mem, ptr_seq_dev_mem_aux, query_len, ref_len, pos_in_query-words_at_once, pos_in_ref-words_at_once, MIN(pos_in_query, query_len), MIN(pos_in_ref, ref_len), n_hits_kept, n_frags_per_block);
                ret = cudaDeviceSynchronize();
                if(ret != cudaSuccess){ fprintf(stderr, "Failed on generating forward frags. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }
            }

            ret = cudaMemcpy(host_left_offset, ptr_left_offset, n_hits_kept * sizeof(uint32_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back left offset. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(host_right_offset, ptr_right_offset, n_hits_kept * sizeof(uint32_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back right offset. Error: %d\n", ret); exit(-1); }

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] frags Q-RC t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 


#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_start);
#endif

            filter_and_write_frags(filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, n_hits_kept, out, 'r', ref_len, min_length);

#ifdef SHOWTIME
            clock_gettime(CLOCK_MONOTONIC, &HD_end);
            time_seconds += ( (uint64_t) HD_end.tv_sec - (uint64_t) HD_start.tv_sec ) ;
            time_nanoseconds += ( (uint64_t) HD_end.tv_nsec - (uint64_t) HD_start.tv_nsec );
            time_seconds *= BILLION;
            fprintf(stdout, "[INFO] filterFrags Q-RC t= %" PRIu64 " ns\n", time_seconds + time_nanoseconds);
            time_seconds = 0;
            time_nanoseconds = 0;
#endif 

        }

        // Restart reference for next query section
        pos_in_ref = 0;
        ++split;

    }

    // Close file where frags are written
    fclose(out);
    fprintf(stdout, "[INFO] Completed all executions\n");

    fclose(query);
    fclose(ref);

    cudaFree(data_mem);
    cudaFreeHost(host_pinned_mem);
    return 0;
}

void print_header(FILE * out, uint32_t query_len, uint32_t ref_len){

    fprintf(out, "All by-Identity Ungapped Fragments (Hits based approach)\n");
    fprintf(out, "[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>\n");
    fprintf(out, "SeqX filename : undef\n");
    fprintf(out, "SeqY filename : undef\n");
    fprintf(out, "SeqX name : undef\n");
    fprintf(out, "SeqY name : undef\n");
    fprintf(out, "SeqX length : %" PRIu32"\n", query_len);
    fprintf(out, "SeqY length : %" PRIu32"\n", ref_len);
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

float factor_chooser(uint64_t total_memory)
{
    // Hits are generated quadratically -> linearly more memory -> quadratically more hits
    // TODO test these
    if(total_memory <=  5000000000) return 0.20;
    if(total_memory <=  7000000000) return 0.20;
    if(total_memory <= 10000000000) return 0.20;
    if(total_memory <= 13000000000) return 0.20;
    if(total_memory <= 20000000000) return 0.17;
    return 0.15;
}

uint64_t memory_allocation_chooser(uint64_t total_memory)
{
    if(total_memory <= 4340179200) return 100*1024*1024;
    else if(total_memory <= 6442450944) return 150*1024*1024;
    else if(total_memory <= 8689934592) return 200*1024*1024;
    return 300*1024*1024;
}

char * dump_memory_region(char * ptr_pointer, uint64_t size)
{
    char * anything = (char *) malloc(size*sizeof(char));
    int ret = cudaMemcpy(anything, ptr_pointer, size * sizeof(char), cudaMemcpyDeviceToHost);
    if(ret != cudaSuccess){ fprintf(stderr, "Dumping data region. Error: %d\n", ret); exit(-1); }
    return anything;
}





