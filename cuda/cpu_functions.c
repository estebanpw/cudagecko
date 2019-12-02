#include "cpu_functions.h"

void terror(const char * s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

void filter_and_write_frags(uint64_t * filtered_hits_x, uint64_t * filtered_hits_y, uint64_t * host_left_offset, uint64_t * host_right_offset, uint64_t n_frags, FILE * out, char strand, uint64_t ref_len){

    
    uint64_t current = 0;

    uint64_t xStart = filtered_hits_x[current] - host_left_offset[current];
    uint64_t xEnd = filtered_hits_x[current] + host_right_offset[current];
    uint64_t yStart = filtered_hits_y[current] - host_left_offset[current];
    uint64_t yEnd = filtered_hits_y[current] + host_right_offset[current];
    uint64_t curr_l = xEnd - xStart;

    uint64_t next_xStart, next_yStart, next_xEnd, next_yEnd, next_l;

    uint64_t max_id = 0;
    uint64_t written_frags = 0;

    while(current + 1 < n_frags){

        next_xStart = filtered_hits_x[current+1] - host_left_offset[current+1];
        next_xEnd = filtered_hits_x[current+1] + host_right_offset[current+1];
        next_yStart = filtered_hits_y[current+1] - host_left_offset[current+1];
        next_yEnd = filtered_hits_y[current+1] + host_right_offset[current+1];
        next_l = next_xEnd - next_xStart;

        // If they are overlapping (on both x and y)
        if(xStart <= next_xEnd && next_xStart <= xEnd && yStart <= next_yEnd && next_yStart <= yEnd){

            // If the new one is bigger
            if(next_l > curr_l){
                curr_l = next_l;
                max_id = current+1;
            }

        }else{
            //printf("%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64" l:%"PRIu64" does not overlap with \n", xStart, yStart, xEnd, yEnd, xEnd-xStart);
            //printf("%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64" l:%"PRIu64"\n", next_xStart, next_yStart, next_xEnd, next_yEnd, next_xEnd-next_xStart);
            uint64_t best_xStart = filtered_hits_x[max_id] - host_left_offset[max_id];
            uint64_t best_xEnd = filtered_hits_x[max_id] + host_right_offset[max_id];
            uint64_t best_yStart = filtered_hits_y[max_id] - host_left_offset[max_id];
            uint64_t best_yEnd = filtered_hits_y[max_id] + host_right_offset[max_id];
            uint64_t best_l = best_xEnd - best_xStart;

            //if(frag.strand=='r'){
			//frag.yStart = ytotal - frag.yStart - 1;
			//frag.yEnd = ytotal - frag.yEnd - 1;

            //printf("So I write: Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",f,0,%"PRIu64",75,75,0.75,0.75,0,0\n", best_xStart, best_yStart, best_xEnd, best_yEnd, best_l);
            if(strand == 'f'){
                fprintf(out, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,0,%"PRIu64",75,75,0.75,0.75,0,0\n", best_xStart, best_yStart, best_xEnd, best_yEnd, strand, best_l);
            } 
            else {
                
                best_yStart = ref_len - best_yStart - 1;
                best_yEnd = ref_len - best_yEnd - 1;
                fprintf(out, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,0,%"PRIu64",75,75,0.75,0.75,0,0\n", best_xStart, best_yStart, best_xEnd, best_yEnd, strand, best_l);
            }
            max_id = current+1;

            xStart = filtered_hits_x[max_id] - host_left_offset[max_id];
            xEnd = filtered_hits_x[max_id] + host_right_offset[max_id];
            yStart = filtered_hits_y[max_id] - host_left_offset[max_id];
            yEnd = filtered_hits_y[max_id] + host_right_offset[max_id];
            curr_l = xEnd - xStart;

            ++written_frags;
        }

        ++current;
    }
    fprintf(stdout, "[INFO] Remaining frags %"PRIu64" out of %"PRIu64" on strand %c\n", written_frags, n_frags, strand);
    

}

uint64_t generate_hits(uint64_t words_at_once, uint64_t * diagonals, Hit * hits, uint64_t * keys_x, 
    uint64_t * keys_y, uint64_t * values_x, uint64_t * values_y, uint64_t query_len, uint64_t ref_len){

    // Nota: para generar TODOS los hits hay que tener en cuenta que si hay hits repetidos en
    // ambos diccionarios se debe volver atrás cuando se encuentra uno distinto
    // Si no, solo saldra ruido horizontal o vertical 

    uint64_t id_x = 0, id_y = 0, n_hits_found = 0;
    uint64_t diff_offset = ref_len;
    uint64_t diag_len = MAX(query_len, ref_len);
    
    while(id_x < words_at_once || id_y < words_at_once) {

        // Compare
        if(id_x == words_at_once || id_y == words_at_once){ break; }
        
        
        if(keys_x[id_x] == keys_y[id_y] && values_x[id_x] != 0xFFFFFFFFFFFFFFFF && values_y[id_y] != 0xFFFFFFFFFFFFFFFF) {
            
            // This is a hit
            //printf("Made hit: ");
            hits[n_hits_found].p1 = values_x[id_x];
            hits[n_hits_found].p2 = values_y[id_y];
            // Compute diagonal value with associated x value -> (x - y) * Ld + x 
            // This is good enough for sequences length up to 2,147,483,648 bp 
            diagonals[n_hits_found] =  ((diff_offset + values_x[id_x]) - values_y[id_y]) * diag_len + (diff_offset + values_x[id_x]);

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


uint64_t filter_hits(uint64_t * diagonals, uint64_t * indexing_numbers, Hit * hits, uint64_t * filtered_hits_x, uint64_t * filtered_hits_y, uint64_t n_hits_found){
    
    if(n_hits_found == 0) return 0;
    int64_t diagonal = (int64_t) hits[indexing_numbers[0]].p1 - (int64_t) hits[indexing_numbers[0]].p2;
    uint64_t last_position = hits[indexing_numbers[0]].p1, t = 1, t_kept = 0;

    while (t < (n_hits_found-1)) {

        if(diagonal == ((int64_t) hits[indexing_numbers[t]].p1 - (int64_t) hits[indexing_numbers[t]].p2) && hits[indexing_numbers[t]].p1 < (last_position+63)){
            
        }else{

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

void init_args(int argc, char ** av, FILE ** query, unsigned * selected_device, FILE ** ref, FILE ** out, unsigned * write){
    
    int pNum = 0;
    char * p1 = NULL, * p2 = NULL;
    char outname[2048]; outname[0] = '\0';
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           GECKOCUDA -query [file] -ref [file] -dev [device]\n");
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