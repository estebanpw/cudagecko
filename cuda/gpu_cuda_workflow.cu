
// Standard utilities and common systems includes
#include "kernels.cuh"
#include "cpu_functions.c"
#include "cub/cub.cuh"

//#define DIMENSION 1000


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

//uint64_t filter_hits(hit_advanced * hits_in, ulong kmer_size, ulong n_hits_found);

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


    uint64_t * keys, * values, * keys_buf, * values_buf;
    fprintf(stdout, "[INFO] I will use %"PRIu64" bytes for hash (you can have %"PRIu64" entries) and their buffers\n", words_at_once * 2 * sizeof(Word), words_at_once);

    

    
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

    char * query_seq_host = (char *) malloc(query_len * sizeof(char));
    char * ref_seq_host = (char *) malloc(ref_len * sizeof(char));

    if(query_seq_host == NULL || ref_seq_host == NULL) terror("Could not allocate memory for sequences in host");

    fprintf(stdout, "[INFO] Loading query\n");
    load_seq(query, query_seq_host);
    fprintf(stdout, "[INFO] Loading reference\n");
    load_seq(ref, ref_seq_host);

    
     
    
    clock_t begin;

    // Pointer to device memory allocating the query sequence
    char * seq_dev_mem = NULL, * seq_dev_mem_aux = NULL;


    // Allocate memory in host to download kmers and store hits
    uint64_t * dict_x_keys, * dict_x_values, * dict_y_keys, * dict_y_values;
    dict_x_keys = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    dict_x_values = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    if(dict_x_keys == NULL || dict_x_values == NULL) { fprintf(stderr, "Allocating for kmer download in query. Error: %d\n", ret); exit(-1); }
    dict_y_keys = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    dict_y_values = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    if(dict_y_keys == NULL || dict_y_values == NULL) { fprintf(stderr, "Allocating for kmer download in ref. Error: %d\n", ret); exit(-1); }
    Hit * hits = (Hit *) malloc(words_at_once*sizeof(Hit));
    uint64_t * filtered_hits_x = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    uint64_t * filtered_hits_y = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    uint64_t * device_filt_hits_x, * device_filt_hits_y, * left_offset, * right_offset;
    uint64_t * host_left_offset = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    uint64_t * host_right_offset = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    if(host_left_offset == NULL || host_right_offset == NULL) terror("Could not allocate host offsets");
    if(hits == NULL || filtered_hits_x == NULL || filtered_hits_y == NULL) terror("Could not allocate hits");
    uint64_t * diagonals = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    uint64_t * device_diagonals, * device_diagonals_buf;
    uint64_t * device_hits, * device_hits_buf; // These will actually be just indices to redirect the hits sorting
    uint64_t * ascending_numbers = (uint64_t *) malloc(words_at_once*sizeof(uint64_t)); for(i=0; i<words_at_once; i++) ascending_numbers[i] = i;
    uint64_t * indexing_numbers = (uint64_t *) malloc(words_at_once*sizeof(uint64_t));
    if(hits == NULL) { fprintf(stderr, "Allocating for hits download. Error: %d\n", ret); exit(-1); }


    ////////////////////////////////////////////////////////////////////////////////
    // Read the query and reference in blocks
    ////////////////////////////////////////////////////////////////////////////////

    
    fprintf(out, "All by-Identity Ungapped Fragments (Hits based approach)\n");
    fprintf(out, "[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>\n");
    fprintf(out, "SeqX filename : undef\n");
    fprintf(out, "SeqY filename : undef\n");
    fprintf(out, "SeqX name : undef\n");
    fprintf(out, "SeqY name : undef\n");
    fprintf(out, "SeqX length : %"PRIu64"\n", query_len);
    fprintf(out, "SeqY length : %"PRIu64"\n", ref_len);
    fprintf(out, "Min.fragment.length : undef\n");
    fprintf(out, "Min.Identity : undef\n");
    fprintf(out, "Tot Hits (seeds) : undef\n");
    fprintf(out, "Tot Hits (seeds) used: undef\n");
    fprintf(out, "Total fragments : undef\n");
    fprintf(out, "========================================================\n");
    fprintf(out, "Total CSB: 0\n");
    fprintf(out, "========================================================\n");
    fprintf(out, "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");

    int split = 0;
    uint64_t pos_in_query = 0, pos_in_ref = 0;
    while(pos_in_query < query_len){

        // Allocate memory in device for sequence chunk
        // We have to this here since later on we will have to free all memory to load the hits
        ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char));
        if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for query sequence in device (Attempted %"PRIu64" bytes). Error: %d\n", words_at_once * sizeof(char), ret); exit(-1); }



        // Allocate words table
        ret = cudaMalloc(&keys, words_at_once * sizeof(uint64_t));
        if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (1). Error: %d\n", ret); exit(-1); }
        ret = cudaMalloc(&values, words_at_once * sizeof(uint64_t));
        if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (2). Error: %d\n", ret); exit(-1); }
        ret = cudaMalloc(&keys_buf, words_at_once * sizeof(uint64_t));
        if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (3). Error: %d\n", ret); exit(-1); }
        ret = cudaMalloc(&values_buf, words_at_once * sizeof(uint64_t));
        if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (4). Error: %d\n", ret); exit(-1); }

        fprintf(stdout, "[EXECUTING] Running split %d -> (%d%%)\n", split, (int)((100*pos_in_query)/query_len));
        uint64_t items_read = MIN(query_len - pos_in_query, words_at_once);


        ////////////////////////////////////////////////////////////////////////////////
        // Run kmers for query
        ////////////////////////////////////////////////////////////////////////////////
        
        // Load sequence chunk into ram
        ret = cudaMemcpy(seq_dev_mem, &query_seq_host[pos_in_query], items_read, cudaMemcpyHostToDevice);
        if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device. Error: %d\n", ret); exit(-1); }

        // Run kmers
        ret = cudaMemset(keys, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
        ret = cudaMemset(values, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
        //number_of_blocks = (((items_read - KMER_SIZE + 1)) / (threads_number*4));
        //kernel_register_fast_hash_rotational<<<number_of_blocks, threads_number>>>(keys, values, seq_dev_mem, pos_in_query);
        number_of_blocks = ((items_read - KMER_SIZE + 1));
        kernel_index_global32<<<number_of_blocks, threads_number>>>(keys, values, seq_dev_mem, pos_in_query);
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

        // Print hits for debug
        //for(i=0; i<items_read; i++){
        //    fprintf(out, "%"PRIu64"\n", dict_x_values[i]);
        //}

        ret = cudaFree(d_temp_storage);
        if(ret != cudaSuccess){ fprintf(stderr, "Bad free of temp storage (1): %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

        pos_in_query += words_at_once;

        cudaFree(seq_dev_mem);
        cudaFree(keys);
        cudaFree(values);
        cudaFree(keys_buf);
        cudaFree(values_buf);        

        ////////////////////////////////////////////////////////////////////////////////
        // Run the reference blocks
        ////////////////////////////////////////////////////////////////////////////////

        while(pos_in_ref < ref_len){

            items_read = MIN(ref_len - pos_in_ref, words_at_once);

            ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char));

            // Allocate words table
            ret = cudaMalloc(&keys, words_at_once * sizeof(uint64_t));
            if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (1). Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&values, words_at_once * sizeof(uint64_t));
            if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (2). Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&keys_buf, words_at_once * sizeof(uint64_t));
            if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (3). Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&values_buf, words_at_once * sizeof(uint64_t));
            if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for table in device (4). Error: %d\n", ret); exit(-1); }

            // Load sequence chunk into ram
            ret = cudaMemcpy(seq_dev_mem, &ref_seq_host[pos_in_ref], items_read, cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device. Error: %d\n", ret); exit(-1); }

            // Run kmers
            ret = cudaMemset(keys, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            ret = cudaMemset(values, 0xFFFFFFFF, words_at_once * sizeof(uint64_t));
            //number_of_blocks = (((items_read - KMER_SIZE + 1)) / (threads_number*4)); 
            //kernel_register_fast_hash_rotational<<<number_of_blocks, threads_number>>>(keys, values, seq_dev_mem, pos_in_ref);
            number_of_blocks = ((items_read - KMER_SIZE + 1));
            kernel_index_global32<<<number_of_blocks, threads_number>>>(keys, values, seq_dev_mem, pos_in_ref);
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

            
            //read_kmers(query_len, query_seq_host, dict_x_keys, dict_x_values);
            //Qsort(dict_x_keys, dict_x_values, 0, (int64_t) query_len);
            //for(i=0; i<words_at_once; i++) printf("%" PRIu64" %"PRIu64"\n", dict_x_keys[i], dict_x_values[i]);
            //read_kmers(ref_len, ref_seq_host, dict_y_keys, dict_y_values);
            //Qsort(dict_y_keys, dict_y_values, 0, (int64_t) ref_len);
            //for(i=0; i<words_at_once; i++) printf("%" PRIu64" %"PRIu64"\n", dict_x_keys[i], dict_x_values[i]);

            cudaFree(seq_dev_mem);
            cudaFree(keys);
            cudaFree(values);
            cudaFree(keys_buf);
            cudaFree(values_buf);


            uint64_t n_hits_found = generate_hits(words_at_once, diagonals, hits, dict_x_keys, dict_y_keys, dict_x_values, dict_y_values, query_len, ref_len);
            
            fprintf(stdout, "[INFO] Generated %"PRIu64" hits on split %d -> (%d%%)\n", n_hits_found, split, (int)((100*MIN(pos_in_ref, ref_len))/ref_len));

            // Print hits for debug
            //for(i=0; i<items_read; i++){
                //fprintf(out, "%"PRIu64"\n", dict_x_values[i]);
            //}
            //for(i=0; i<n_hits_found; i++){
                //printf("%"PRIu64"\n", diagonals[i]);
                //if(hits[i].p1 > 368000 && hits[i].p2 < 390000)
                    //fprintf(out, "Frag,d:%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,32,32,32,1.0,1.0,0,0\n", (int64_t) hits[i].p1 - (int64_t) hits[i].p2, hits[i].p1, hits[i].p2, hits[i].p1+32, hits[i].p2+32);
            //}

            ////////////////////////////////////////////////////////////////////////////////
            // Sort hits for the current split
            ////////////////////////////////////////////////////////////////////////////////

            ret = cudaMalloc(&device_diagonals, words_at_once * sizeof(uint64_t));
            if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (1). Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&device_diagonals_buf, words_at_once * sizeof(uint64_t));
            if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (2). Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&device_hits, words_at_once * sizeof(uint64_t));
            if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (3). Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&device_hits_buf, words_at_once * sizeof(uint64_t));
            if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate memory for hits in device (4). Error: %d\n", ret); exit(-1); }

            // We will actually sort the diagonals with associated values 0,1,2... to n and use these to index the hits array
            ret = cudaMemcpy(device_hits, ascending_numbers, n_hits_found*sizeof(uint64_t), cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Uploading device hits. Error: %d\n", ret); exit(-1); }

            ret = cudaMemcpy(device_diagonals, diagonals, n_hits_found*sizeof(uint64_t), cudaMemcpyHostToDevice);
            if(ret != cudaSuccess){ fprintf(stderr, "Uploading device diagonals. Error: %d\n", ret); exit(-1); }

            cub::DoubleBuffer<uint64_t> d_diagonals(device_diagonals, device_diagonals_buf);
            cub::DoubleBuffer<uint64_t> d_hits(device_hits, device_hits_buf);
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
            ret = cudaMemcpy(indexing_numbers, device_hits, n_hits_found*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device hits. Error: %d\n", ret); exit(-1); }

            ret = cudaMemcpy(diagonals, device_diagonals, n_hits_found*sizeof(uint64_t), cudaMemcpyDeviceToHost);
            if(ret != cudaSuccess){ fprintf(stderr, "Downloading device diagonals. Error: %d\n", ret); exit(-1); }


            uint64_t n_hits_kept = filter_hits(diagonals, indexing_numbers, hits, filtered_hits_x, filtered_hits_y, n_hits_found);

            fprintf(stdout, "[INFO] Remaining hits %"PRIu64"\n", n_hits_kept);

            //for(i=0; i<n_hits_kept; i++){
                //printf("%"PRIu64"\n", diagonals[i]);
                //fprintf(out, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,32,32,32,1.0,1.0,0,0\n", filtered_hits_x[i], filtered_hits_y[i], filtered_hits_x[i]+32, filtered_hits_y[i]+32);
            //}
            
            ret = cudaFree(d_temp_storage);
            ret = cudaFree(device_hits);
            ret = cudaFree(device_diagonals);
            ret = cudaFree(device_diagonals_buf);
            ret = cudaFree(device_hits_buf);

            ////////////////////////////////////////////////////////////////////////////////
            // Generate frags
            ////////////////////////////////////////////////////////////////////////////////

            // Allocate both sequences
            ret = cudaMalloc(&seq_dev_mem, words_at_once * sizeof(char)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for query sequence in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&seq_dev_mem_aux, words_at_once * sizeof(char)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for ref sequence in device. Error: %d\n", ret); exit(-1); }

            ret = cudaMalloc(&device_filt_hits_x, words_at_once * sizeof(uint64_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device filtered hits query. Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&device_filt_hits_y, words_at_once * sizeof(uint64_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device filtered hits ref. Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&left_offset, words_at_once * sizeof(uint64_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device offset left frags query. Error: %d\n", ret); exit(-1); }
            ret = cudaMalloc(&right_offset, words_at_once * sizeof(uint64_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not allocate for device offset right frags query. Error: %d\n", ret); exit(-1); }

            ret = cudaMemcpy(seq_dev_mem, &query_seq_host[pos_in_query-words_at_once], MIN(query_len - (pos_in_query - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy query sequence to device for frags. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(seq_dev_mem_aux, &ref_seq_host[pos_in_ref-words_at_once], MIN(ref_len - (pos_in_ref - words_at_once), words_at_once), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy ref sequence to device for frags. Error: %d\n", ret); exit(-1); }
            
            ret = cudaMemcpy(device_filt_hits_x, filtered_hits_x, n_hits_kept * sizeof(uint64_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits x in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(device_filt_hits_y, filtered_hits_y, n_hits_kept * sizeof(uint64_t), cudaMemcpyHostToDevice); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy filtered hits y in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(left_offset, 0x0, n_hits_kept * sizeof(uint64_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy left offset in device. Error: %d\n", ret); exit(-1); }
            ret = cudaMemset(right_offset, 0x0, n_hits_kept * sizeof(uint64_t)); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy right offset in device. Error: %d\n", ret); exit(-1); }


            number_of_blocks = n_hits_kept; 
            //number_of_blocks = 20; // REMOVE !!
            kernel_frags_forward_register<<<number_of_blocks, threads_number>>>(device_filt_hits_x, device_filt_hits_y, left_offset, right_offset, seq_dev_mem, seq_dev_mem_aux, query_len, ref_len);
            ret = cudaDeviceSynchronize();
            if(ret != cudaSuccess){ fprintf(stderr, "Failed on generating forward frags. Error: %d -> %s\n", ret, cudaGetErrorString(cudaGetLastError())); exit(-1); }

            ret = cudaMemcpy(host_left_offset, left_offset, n_hits_kept * sizeof(uint64_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back left offset. Error: %d\n", ret); exit(-1); }
            ret = cudaMemcpy(host_right_offset, right_offset, n_hits_kept * sizeof(uint64_t), cudaMemcpyDeviceToHost); if(ret != cudaSuccess){ fprintf(stderr, "Could not copy back right offset. Error: %d\n", ret); exit(-1); }

            cudaFree(seq_dev_mem);
            cudaFree(seq_dev_mem_aux);
            cudaFree(device_filt_hits_x);
            cudaFree(device_filt_hits_y);
            cudaFree(left_offset);
            cudaFree(right_offset);

            /*
            for(i=0; i<n_hits_kept; i++){
                //printf("Hit (%"PRIu64", %"PRIu64") yields offsets (%"PRIu64", %"PRIu64") diff: \n", filtered_hits_x[i], filtered_hits_y[i], host_left_offset[i], host_right_offset[i]);
                uint64_t xStart = filtered_hits_x[i] - host_left_offset[i];
                uint64_t xEnd = filtered_hits_x[i] + host_right_offset[i];
                uint64_t yStart = filtered_hits_y[i] - host_left_offset[i];
                uint64_t yEnd = filtered_hits_y[i] + host_right_offset[i];
                
                if((xEnd - xStart) > 99){
                    fprintf(out, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,%"PRIu64",75,75,0.75,0.75,0,0\n", xStart, yStart, xEnd, yEnd, xEnd - xStart);
                }
            
                
            }
            */
            filter_and_write_frags(filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, n_hits_kept, out);

        }

        // Restart the reference for every block in query
        pos_in_ref = 0;




        ++split;
    }

    fclose(out);
    
    
    
    
    

    fprintf(stdout, "[INFO] Completed\n");

    fclose(query);
    fclose(ref);
    free(query_seq_host);
    free(ref_seq_host);
    free(diagonals);
    free(ascending_numbers);
    

    return 0;
}




