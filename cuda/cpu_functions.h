#include "structs.h"

uint64_t generate_hits(uint64_t words_at_once, uint64_t * diagonals, Hit * hits, uint64_t * keys_x, uint64_t * keys_y, uint64_t * values_x, uint64_t * values_y, uint64_t query_len, uint64_t ref_len);
uint64_t filter_hits(uint64_t * diagonals, uint64_t * indexing_numbers, Hit * hits, uint64_t * filtered_hits_x, uint64_t * filtered_hits_y, uint64_t n_hits_found);
void read_kmers(uint64_t query_l, char * seq_x, uint64_t * keys_x, uint64_t * values_x);
void init_args(int argc, char ** av, FILE ** query, unsigned * selected_device, FILE ** ref, FILE ** out, unsigned * write);
void perfect_hash_to_word(char * word, uint64_t hash, uint64_t k);
void print_kmers_to_file(uint64_t * keys, uint64_t * values, uint64_t table_size, FILE * fout);
char * get_dirname(char * path);
char * get_basename(char * path);
uint64_t load_seq(FILE * f, char * seq);
uint64_t get_seq_len(FILE * f);
void terror(const char * s);
void Qsort(uint64_t * keys, uint64_t * values, int64_t x, int64_t y);
