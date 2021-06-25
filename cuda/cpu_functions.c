#include "cpu_functions.h"

void terror(const char * s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

void build_frag(uint32_t * xStart, uint32_t * xEnd, uint32_t * yStart, uint32_t * yEnd, uint32_t * curr_l, char strand,
    uint32_t * filtered_hits_x, uint32_t * filtered_hits_y, uint32_t * host_left_offset, uint32_t * host_right_offset, uint32_t id){
    if(strand == 'f'){
        *xStart = filtered_hits_x[id] - host_left_offset[id];
        *xEnd = filtered_hits_x[id] + host_right_offset[id];
        *yStart = filtered_hits_y[id] - host_left_offset[id];
        *yEnd = filtered_hits_y[id] + host_right_offset[id];
    }else{
        *xStart = filtered_hits_x[id] - host_left_offset[id];
        *xEnd = filtered_hits_x[id] + host_right_offset[id];
        *yStart = filtered_hits_y[id] - host_left_offset[id];
        *yEnd = filtered_hits_y[id] + host_right_offset[id];
    }
    //printf("Fraggo %u %u %u %u\n", *xStart, *xEnd, *yStart, *yEnd);
    *curr_l = *xEnd - *xStart;
}

void filter_and_write_frags(uint32_t * filtered_hits_x, uint32_t * filtered_hits_y, uint32_t * host_left_offset, 
    uint32_t * host_right_offset, uint32_t n_frags, FILE * out, char strand, uint32_t ref_len, uint32_t min_length){

    
    uint32_t current = 0;

    uint32_t xStart, xEnd, yStart, yEnd, curr_l;

    build_frag(&xStart, &xEnd, &yStart, &yEnd, &curr_l, strand, filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, current);


    uint32_t next_xStart, next_yStart, next_xEnd, next_yEnd, next_l;

    uint32_t max_id = 0;
    uint32_t written_frags = 0;
    uint32_t last_unwritten = 0;

    if(n_frags == 1 && xEnd-xStart >= min_length){
        if(strand == 'f'){

            uint32_t score = (uint32_t)(curr_l * MIN_P_IDENT);
            fprintf(out, "Frag,%" PRIu32",%" PRIu32",%" PRIu32",%" PRIu32",%c,0,%" PRIu32",%" PRIu32",%" PRIu32",%.2f,%.2f,0,0\n", xStart, yStart, xEnd, yEnd, strand, curr_l, score, score, OUTPUT_P_IDENT, OUTPUT_P_IDENT);
        }else{
            uint32_t best_yStart = ref_len - yStart - 1;
            uint32_t best_yEnd = ref_len - yEnd - 1;
            uint32_t score = (uint32_t)(curr_l * MIN_P_IDENT);
            fprintf(out, "Frag,%" PRIu32",%" PRIu32",%" PRIu32",%" PRIu32",%c,0,%" PRIu32",%" PRIu32",%" PRIu32",%.2f,%.2f,0,0\n", xStart, best_yStart, xEnd, best_yEnd, strand, curr_l, score, score, OUTPUT_P_IDENT, OUTPUT_P_IDENT);

        }
        ++written_frags;
    }


    while(current + 1 < n_frags){

        build_frag(&next_xStart, &next_xEnd, &next_yStart, &next_yEnd, &next_l, strand, filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, current+1);
        //fprintf(stdout, "Looking for [%" PRIu64"]: Frag,%" PRIu64",%" PRIu64",%" PRIu64",%" PRIu64",%c,0,%" PRIu64", his hit [%" PRIu64", %" PRIu64"]\n", current+1, next_xStart, next_yStart, next_xEnd, next_yEnd, strand, next_l, filtered_hits_x[current+1], filtered_hits_y[current+1]);

        /*
        if(strand == 'r'){
            next_yStart = ref_len - next_yStart - 1;
            next_yEnd = ref_len - next_yEnd - 1;
            fprintf(out, "Frag,%" PRIu64",%" PRIu64",%" PRIu64",%" PRIu64",%c,0,%" PRIu64",75,75,0.75,0.75,0,0\n", next_xStart, next_yStart, next_xEnd, next_yEnd, strand, next_l);
        }
        */

        // If they are overlapping (on both x and y)
        if(strand == 'f'){
            int64_t dprev = (int64_t) xStart - (int64_t) yStart;
            int64_t dpost = (int64_t) next_xStart - (int64_t) next_yStart;

            //fprintf(stdout, "Their diags [%"PRId64":%"PRId64"]\n", dprev, dpost);

            if(dprev == dpost && xStart <= next_xEnd && next_xStart <= xEnd && yStart <= next_yEnd && next_yStart <= yEnd){

                // If the new one is bigger
                if(next_l >= curr_l){

                    //fprintf(stdout, "Yea, next one is bigger\n");
                    curr_l = next_l;
                    max_id = current+1;
                    last_unwritten = 1;
                }


            }else{

                //fprintf(stdout, "Not overlapping or different diag now so lets go save [%" PRIu64"]\n", max_id);

                uint32_t best_xStart, best_xEnd, best_yStart, best_yEnd, best_l;

                build_frag(&best_xStart, &best_xEnd, &best_yStart, &best_yEnd, &best_l, strand, filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, max_id);

                if((best_xEnd - best_xStart) >= min_length){
                    
                    // Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY
                    uint32_t score = (uint32_t)(best_l * MIN_P_IDENT);
                    fprintf(out, "Frag,%" PRIu32",%" PRIu32",%" PRIu32",%" PRIu32",%c,0,%" PRIu32",%" PRIu32",%" PRIu32",%.2f,%.2f,0,0\n", best_xStart, best_yStart, best_xEnd, best_yEnd, strand, best_l, score, score, OUTPUT_P_IDENT, OUTPUT_P_IDENT);
                    ++written_frags;
                }
                max_id = current+1;
                last_unwritten = 0;


                
                build_frag(&xStart, &xEnd, &yStart, &yEnd, &curr_l, strand, filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, max_id);
                //fprintf(stdout, "I am [%" PRIu64"]: Frag,%" PRIu64",%" PRIu64",%" PRIu64",%" PRIu64",%c,0,%" PRIu64", his hit [%" PRIu64", %" PRIu64"]\n", max_id, xStart, yStart, xEnd, yEnd, strand, curr_l, filtered_hits_x[max_id], filtered_hits_y[max_id]);

            }
        }

        if(strand == 'r'){

            int64_t dprev = (int64_t) xStart - (int64_t) yStart;
            int64_t dpost = (int64_t) next_xStart - (int64_t) next_yStart;

            //fprintf(stdout, "Their diags [%"PRId64":%"PRId64"]\n", dprev, dpost);

            if(dprev == dpost && xStart <= next_xEnd && next_xStart <= xEnd && yStart <= next_yEnd && next_yStart <= yEnd){

                // If the new one is bigger
                if(next_l >= curr_l){
                    //fprintf(stdout, "Yea, next one is bigger\n");
                    curr_l = next_l;
                    max_id = current+1;
                    last_unwritten = 1;
                }

            }else{


                //fprintf(stdout, "Not overlapping or different diag now so lets go save [%" PRIu64"]\n", max_id);

                uint32_t best_xStart, best_xEnd, best_yStart, best_yEnd, best_l;

                build_frag(&best_xStart, &best_xEnd, &best_yStart, &best_yEnd, &best_l, strand, filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, max_id);

                if((best_xEnd - best_xStart) >= min_length){
                    //int64_t d = best_xStart + best_yStart;
                    best_yStart = ref_len - best_yStart - 1;
                    best_yEnd = ref_len - best_yEnd - 1;
                    uint32_t score = (uint32_t)(best_l * MIN_P_IDENT);
                    fprintf(out, "Frag,%" PRIu32",%" PRIu32",%" PRIu32",%" PRIu32",%c,0,%" PRIu32",%" PRIu32",%" PRIu32",%.2f,%.2f,0,0\n", best_xStart, best_yStart, best_xEnd, best_yEnd, strand, best_l, score, score, OUTPUT_P_IDENT, OUTPUT_P_IDENT);
                    ++written_frags;
                    
                }
                max_id = current+1;
                last_unwritten = 0;

                build_frag(&xStart, &xEnd, &yStart, &yEnd, &curr_l, strand, filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, max_id);
                //fprintf(stdout, "I am [%" PRIu64"]: Frag,%" PRIu64",%" PRIu64",%" PRIu64",%" PRIu64",%c,0,%" PRIu64", his hit [%" PRIu64", %" PRIu64"]\n", max_id, xStart, yStart, xEnd, yEnd, strand, curr_l, filtered_hits_x[max_id], filtered_hits_y[max_id]);
            }
        }

        ++current;
    }

    if(last_unwritten != 0)
    {
        if(strand == 'f')
        {
            uint32_t best_xStart, best_xEnd, best_yStart, best_yEnd, best_l;

            build_frag(&best_xStart, &best_xEnd, &best_yStart, &best_yEnd, &best_l, strand, filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, max_id);

            if((best_xEnd - best_xStart) >= min_length){

                uint32_t score = (uint32_t)(best_l * MIN_P_IDENT);
                fprintf(out, "Frag,%" PRIu32",%" PRIu32",%" PRIu32",%" PRIu32",%c,0,%" PRIu32",%" PRIu32",%" PRIu32",%.2f,%.2f,0,0\n", best_xStart, best_yStart, best_xEnd, best_yEnd, strand, best_l, score, score, OUTPUT_P_IDENT, OUTPUT_P_IDENT);
                ++written_frags;
            }

        }
        else
        {

            uint32_t best_xStart, best_xEnd, best_yStart, best_yEnd, best_l;

            build_frag(&best_xStart, &best_xEnd, &best_yStart, &best_yEnd, &best_l, strand, filtered_hits_x, filtered_hits_y, host_left_offset, host_right_offset, max_id);

            if((best_xEnd - best_xStart) >= min_length){
                best_yStart = ref_len - best_yStart - 1;
                best_yEnd = ref_len - best_yEnd - 1;
                uint32_t score = (uint32_t)(best_l * MIN_P_IDENT);
                fprintf(out, "Frag,%" PRIu32",%" PRIu32",%" PRIu32",%" PRIu32",%c,0,%" PRIu32",%" PRIu32",%" PRIu32",%.2f,%.2f,0,0\n", best_xStart, best_yStart, best_xEnd, best_yEnd, strand, best_l, score, score, OUTPUT_P_IDENT, OUTPUT_P_IDENT);
                ++written_frags;

            }

        }
    }

#ifdef SHOWTIME
    fprintf(stdout, "[INFO] Remaining frags %" PRIu32" out of %" PRIu32" on strand %c\n", written_frags, n_frags, strand);
#endif

}

uint32_t generate_hits_fast(uint32_t max_hits, uint64_t * diagonals, Hit * hits, uint64_t * keys_x, 
    uint64_t * keys_y, uint32_t * values_x, uint32_t * values_y, uint32_t items_x, uint32_t items_y, uint32_t query_len, uint32_t ref_len){

    // Nota: para generar TODOS los hits hay que tener en cuenta que si hay hits repetidos en
    // ambos diccionarios se debe volver atrás cuando se encuentra uno distinto
    // Si no, solo saldra ruido horizontal o vertical 

    uint64_t id_x = 0, id_y = 0, n_hits_found = 0;
    uint64_t diff_offset = ref_len;
    uint64_t diag_len = MAX(query_len, ref_len);
    //int64_t last_position_y = -1;

    
    while(id_x < items_x || id_y < items_y) {

        // Compare
        if(id_x == items_x || id_y == items_y){ break; }
        
        
        if(keys_x[id_x] == keys_y[id_y] && values_x[id_x] != 0xFFFFFFFF && values_y[id_y] != 0xFFFFFFFF) {
            
            //if(last_position_y == -1) last_position_y = (int64_t) id_y;
            //if(last_position_x == -1) last_position_x = (int64_t) id_x;

            // This is a hit
            //printf("Made hit: ");
            hits[n_hits_found].p1 = values_x[id_x];
            hits[n_hits_found].p2 = values_y[id_y];
            // Compute diagonal value with associated x value -> (x - y) * Ld + x 
            // Le estoy sumando el diff_offset (i.e. lo mas largo que puede ser la y) para que siempre sean positivos (como x es positivo)
            // pues solo la y resta
            // This is good enough for sequences length up to 2,147,483,648 bp 
            diagonals[n_hits_found] =  ((diff_offset + (uint64_t) values_x[id_x]) - (uint64_t) values_y[id_y]) * diag_len + (diff_offset + (uint64_t) values_x[id_x]);

            //printf("Matching hash %" PRIu64" with %" PRIu64" @ (%" PRIu64", %" PRIu64")\n", keys_x[id_x], keys_y[id_y], values_x[id_x], values_y[id_y]);

            ++n_hits_found;
            if(n_hits_found == max_hits){ fprintf(stderr, "Reached maximum limit of hits\n"); }

            ++id_y;

            //printf("next comp is %" PRIu64" with %" PRIu64"\n", keys_x[id_x], keys_y[id_y]);
        }
        else if(keys_x[id_x] < keys_y[id_y]){
            // Hits are lost if the same hashes are on both dictionaries repeated
            // Because once you skip, you skip forever
            //id_y = (uint64_t) last_position_y;
            ++id_x;
        } 
        else {
            ++id_y; 
            
        }
    }

    //printf("Generated %" PRIu64" hits \n", n_hits_found);
    return n_hits_found;

}

uint32_t generate_hits_quadratic(uint32_t words_at_once, uint64_t * diagonals, Hit * hits, uint64_t * keys_x, 
    uint64_t * keys_y, uint32_t * values_x, uint32_t * values_y, uint32_t items_x, uint32_t items_y, uint32_t query_len, uint32_t ref_len){

    // Nota: para generar TODOS los hits hay que tener en cuenta que si hay hits repetidos en
    // ambos diccionarios se debe volver atrás cuando se encuentra uno distinto
    // Si no, solo saldra ruido horizontal o vertical 

    uint64_t id_x = 0, id_y = 0, n_hits_found = 0;
    uint64_t diff_offset = ref_len;
    uint64_t diag_len = MAX(query_len, ref_len);
    //int64_t last_position_y = -1;

    
    for(id_x=0; id_x < items_x ; id_x++)
    {

        for(id_y=0; id_y < items_y; id_y++) 
        {

            // Compare
            if(keys_x[id_x] == keys_y[id_y] && values_x[id_x] != 0xFFFFFFFF && values_y[id_y] != 0xFFFFFFFF) {
                
                //if(last_position_y == -1) last_position_y = (int64_t) id_y;
                //if(last_position_x == -1) last_position_x = (int64_t) id_x;

                // This is a hit
                //printf("Made hit: ");
                hits[n_hits_found].p1 = values_x[id_x];
                hits[n_hits_found].p2 = values_y[id_y];
                // Compute diagonal value with associated x value -> (x - y) * Ld + x 
                // Le estoy sumando el diff_offset (i.e. lo mas largo que puede ser la y) para que siempre sean positivos (como x es positivo)
                // pues solo la y resta
                // This is good enough for sequences length up to 2,147,483,648 bp 
                diagonals[n_hits_found] =  ((diff_offset + (uint64_t) values_x[id_x]) - (uint64_t) values_y[id_y]) * diag_len + (diff_offset + (uint64_t) values_x[id_x]);

                //printf("Matching hash %" PRIu64" with %" PRIu64" @ (%" PRIu64", %" PRIu64")\n", keys_x[id_x], keys_y[id_y], values_x[id_x], values_y[id_y]);

                ++n_hits_found;
                if(n_hits_found == words_at_once){ fprintf(stderr, "Reached maximum limit of hits\n"); }

                //printf("next comp is %" PRIu64" with %" PRIu64"\n", keys_x[id_x], keys_y[id_y]);
            }
        }
    }

    //printf("Generated %" PRIu64" hits \n", n_hits_found);
    return n_hits_found;

}




uint32_t generate_hits_sensitive(uint32_t max_hits, uint64_t * diagonals, Hit * hits, uint64_t * keys_x, 
    uint64_t * keys_y, uint32_t * values_x, uint32_t * values_y, uint32_t items_x, uint32_t items_y, uint32_t query_len, uint32_t ref_len, uint32_t max_frequency, int fast){

    // Nota: para generar TODOS los hits hay que tener en cuenta que si hay hits repetidos en
    // ambos diccionarios se debe volver atrás cuando se encuentra uno distinto
    // Si no, solo saldra ruido horizontal o vertical 

    uint64_t id_x = 0, id_y = 0, n_hits_found = 0;
    //uint64_t diff_offset = ref_len; //Necessary only for square root split mode
    //uint64_t diag_len = MAX(query_len, ref_len); //Same as before; mode double 32-bits does not need this
    uint64_t current_hits = 0;
    //int64_t last_position_y = -1;



    
    while(id_x < items_x || id_y < items_y) {

        // Compare
        if(id_x == items_x || id_y == items_y){ break; }

        
        if(keys_x[id_x] == keys_y[id_y] && values_x[id_x] != 0xFFFFFFFF && values_y[id_y] != 0xFFFFFFFF) {

		uint64_t step_x = 1, step_y = 1;
		if(fast > 0) find_consecutive_seeds(id_x, id_y, keys_x, keys_y, values_x, values_y, items_x, items_y, &step_x, &step_y);


            uint64_t curr_id_y;

            for(curr_id_y = id_y; curr_id_y < items_y; curr_id_y += step_y){

                if(keys_x[id_x] != keys_y[curr_id_y] || values_x[id_x] == 0xFFFFFFFF || values_y[curr_id_y] == 0xFFFFFFFF) break;

                hits[n_hits_found].p1 = values_x[id_x];
                hits[n_hits_found].p2 = values_y[curr_id_y];

                // Original way
                //diagonals[n_hits_found] =  ((diff_offset + (uint64_t) values_x[id_x]) - (uint64_t) values_y[curr_id_y]) * diag_len + (diff_offset + (uint64_t) values_x[id_x]);
                //printf("1 x %" PRIu32 " y %" PRIu32 " -> %" PRIu64 "\n", values_x[id_x], values_y[curr_id_y], diagonals[n_hits_found]);
                // New way
                diagonals[n_hits_found] = ((( ref_len + (uint64_t) values_x[id_x]) - (uint64_t) values_y[curr_id_y]) << 32) + (uint64_t) values_x[id_x];
                //printf("2 x %" PRIu32 " y %" PRIu32 " -> %" PRIu64 "\n", values_x[id_x], values_y[curr_id_y], diagonals[n_hits_found]);

                ++n_hits_found;
                ++current_hits;
                if(current_hits == max_frequency) break;
                if(n_hits_found == max_hits){ fprintf(stderr, "Reached maximum limit of hits (max %" PRIu64")\n", n_hits_found); }

            }

            id_x += step_x;
            current_hits = 0;
            
        }
        else if(keys_x[id_x] < keys_y[id_y]){
            ++id_x;
        } 
        else {
            ++id_y; 
        }

    }


    //printf("Generated %" PRIu64" hits \n", n_hits_found);
    return n_hits_found;

}

uint32_t filter_hits_cpu(uint64_t * diagonals, uint32_t * filtered_hits_x, uint32_t * filtered_hits_y, uint32_t n_hits_found){
    if(n_hits_found == 0) return 0;
    uint32_t t_kept = 0;

    for(uint32_t i=0; i<n_hits_found; i++){
        if(diagonals[i] != 0xFFFFFFFFFFFFFFFF){
            filtered_hits_x[t_kept] = (uint32_t) ((diagonals[i] & 0xFFFFFFFF00000000) >> 32);
            filtered_hits_y[t_kept] = (uint32_t) (diagonals[i] & 0x00000000FFFFFFFF);
            
            ++t_kept;
        }
    }
    return t_kept;
}

uint32_t filter_hits_forward(uint64_t * diagonals, uint32_t * indexing_numbers, Hit * hits, uint32_t * filtered_hits_x, uint32_t * filtered_hits_y, uint32_t n_hits_found){
   

 
    if(n_hits_found == 0) return 0;
    int64_t diagonal = (int64_t) hits[indexing_numbers[0]].p1 - (int64_t) hits[indexing_numbers[0]].p2;
    uint32_t last_position = hits[indexing_numbers[0]].p1, t = 1, t_kept = 0;

    // First hit is saved no matter what
    filtered_hits_x[t_kept] = hits[indexing_numbers[0]].p1;
    filtered_hits_y[t_kept] = hits[indexing_numbers[0]].p2;
    ++t_kept;


    while (t < (n_hits_found)) {

        if(diagonal != ((int64_t) hits[indexing_numbers[t]].p1 - (int64_t) hits[indexing_numbers[t]].p2) || hits[indexing_numbers[t]].p1 > (last_position+63)){
            
            last_position = hits[indexing_numbers[t]].p1;
            diagonal = (int64_t) hits[indexing_numbers[t]].p1 - (int64_t) hits[indexing_numbers[t]].p2;
            filtered_hits_x[t_kept] = hits[indexing_numbers[t]].p1;
            filtered_hits_y[t_kept] = hits[indexing_numbers[t]].p2;
            ++t_kept;
        }
        ++t;
    }
    return t_kept;
}

uint32_t filter_hits_reverse(uint64_t * diagonals, uint32_t * indexing_numbers, Hit * hits, uint32_t * filtered_hits_x, uint32_t * filtered_hits_y, uint32_t n_hits_found){
    
    if(n_hits_found == 0) return 0;
    int64_t diagonal = (int64_t) hits[indexing_numbers[0]].p1 - (int64_t) hits[indexing_numbers[0]].p2;
    uint32_t last_position = hits[indexing_numbers[0]].p1, t = 1, t_kept = 0;

    filtered_hits_x[t_kept] = hits[indexing_numbers[0]].p1;
    filtered_hits_y[t_kept] = hits[indexing_numbers[0]].p2;
    ++t_kept;

    while (t < (n_hits_found)) {

        //fprintf(stdout, "dprev: % "PRId64", dnew: %" PRId64"\n", diagonal, (int64_t) hits[indexing_numbers[t]].p1 + (int64_t) hits[indexing_numbers[t]].p2);
        //printf("d %" PRId64 " xy: % " PRIu32 " %" PRIu32 "\n", (int64_t) hits[indexing_numbers[t]].p1 - (int64_t) hits[indexing_numbers[t]].p2,  hits[indexing_numbers[t]].p1, hits[indexing_numbers[t]].p2);

        if(diagonal != ((int64_t) hits[indexing_numbers[t]].p1 - (int64_t) hits[indexing_numbers[t]].p2) || hits[indexing_numbers[t]].p1 > (last_position+63)){

            last_position = hits[indexing_numbers[t]].p1;
            diagonal = (int64_t) hits[indexing_numbers[t]].p1 - (int64_t) hits[indexing_numbers[t]].p2;
            filtered_hits_x[t_kept] = hits[indexing_numbers[t]].p1;
            filtered_hits_y[t_kept] = hits[indexing_numbers[t]].p2;
            ++t_kept;
        }
        ++t;
    }
    return t_kept;
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

void init_args(int argc, char ** av, FILE ** query, unsigned * selected_device, FILE ** ref, FILE ** out, uint32_t * min_length, int * fast, uint32_t * max_frequency, float * factor, uint32_t * n_frags_per_block){
    
    int pNum = 0;
    char * p1 = NULL, * p2 = NULL;
    char outname[2048]; outname[0] = '\0';
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           GPUGECKO -query [file] -ref [file] -dev [device]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -len        Minimum length of a frag (default 32)\n");
            fprintf(stdout, "           -max_freq   [only works in --sensitive] Maximum frequency per hit (default: unlimited)\n");
            fprintf(stdout, "                       (fast mode can skip highly repeated seeds)\n");
            fprintf(stdout, "           -factor     Fraction of GPU Ram dedicated to words (default: 0.125)\n");
            fprintf(stdout, "                       (The bigger the fraction, the faster it will run - however highly similar sequences\n");
            fprintf(stdout, "                       such as human and gorilla chromosomes require a smaller fractions because of the\n");
            fprintf(stdout, "                       huge number of hits that are generated)\n");
            fprintf(stdout, "           --fast      Runs in fast mode (sensitive is default)\n");
            fprintf(stdout, "           --hyperfast Runs in hyper fast mode \n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
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

        if(strcmp(av[pNum], "-seeds_pb") == 0){
            *n_frags_per_block = (uint32_t) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 1) { fprintf(stderr, "Seeds per block must be >0\n"); exit(-1); }
        }

		if(strcmp(av[pNum], "-len") == 0){
            *min_length = (uint32_t) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 1) { fprintf(stderr, "Length must be >0\n"); exit(-1); }
        }


        if(strcmp(av[pNum], "-factor") == 0){
            *factor = atof(av[pNum+1]);
            if(atof(av[pNum+1]) <= 0) { fprintf(stderr, "Factor must be >0\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "--fast") == 0){
            *fast = 1;
        }

        if(strcmp(av[pNum], "--hyperfast") == 0){
            *fast = 2;
        }

        if(strcmp(av[pNum], "-max_freq") == 0){
            *max_frequency = (uint32_t) atoi(av[pNum+1]);
            if(atoi(av[pNum+1]) < 1) { fprintf(stderr, "Frequency must be >0\n"); exit(-1); }
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
    if(*fast != 0 && *max_frequency != 0){ fprintf(stderr, "Sensitive mode must be enabled to use max frequency per hits (use --sensitive)\n"); exit(-1);}
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
        fprintf(fout, "#%" PRIu64", %s\n", i, word);
        fprintf(fout, "#%" PRIu64", %" PRIu64" at %"  PRIu64"\n", i, keys[i], values[i]);
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

uint32_t get_seq_len(FILE * f) {
    char c = '\0';
    uint32_t l = 0;

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

uint32_t load_seq(FILE * f, char * seq) {
    char c = '\0';
    uint32_t l = 0;

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

uint32_t from_ram_load(char * ram, char * dst, uint32_t size) {
    char c = ram[0];
    uint32_t l = 0;
    uint32_t real_l = 0;

	while(c != '>')  c = ram[l++];
	while(c != '\n') c = ram[l++];

    while(l<size){
        c = ram[l++];

        if(c == '>'){
			dst[real_l++] = 'N';
            while(c != '\n') c = ram[l++];
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            dst[real_l++] = c;
        }
    }


    return real_l;
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

uint64_t realign_address(uint64_t address, uint64_t align)
{

    if(address % align == 0) return address;

    return address + align - (address % align);

}

void find_consecutive_seeds(uint64_t i, uint64_t j, uint64_t * keys_x, uint64_t * keys_y, uint32_t * values_x, uint32_t * values_y, uint32_t items_x, uint32_t items_y, uint64_t * step_x, uint64_t * step_y)
{
	uint64_t i_cons = i + 1, j_cons = j + 1;

	while(i_cons < items_x)
	{

		if(keys_x[i_cons] != keys_x[i] || keys_x[i_cons] == 0xFFFFFFFFFFFFFFFF) break;
		++i_cons;

	}

	while(j_cons < items_y)
	{

		if(keys_y[j_cons] != keys_y[j] || keys_y[j_cons] == 0xFFFFFFFFFFFFFFFF) break;
		++j_cons;

	}

	// This is for uniform mode with fixed factor

	/*

	*step_x = (i_cons - i) / 10;
	*step_y = (j_cons - j) / 10;
	
	*/

	// This is for uniform mode relative to the number of seeds
	// This is MAX / MIN = RATIO; where MIN is the step in Y and RATIO is the step in X

	uint64_t xn = (i_cons - i);
	uint64_t yn = (j_cons - j);

	*step_x = sqrt(xn);
	*step_y = sqrt(yn);

	// This is the same but reducing some numbers

	*step_x = *step_x * 2;
	*step_y = *step_y * 2;


	if(*step_x == 0) *step_x = 1;
	if(*step_y == 0) *step_y = 1;



}








