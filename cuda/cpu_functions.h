#include "structs.h"

void build_frag(uint64_t * xStart, uint64_t * xEnd, uint64_t * yStart, uint64_t * yEnd, uint64_t * curr_l, char strand, uint64_t * filtered_hits_x, uint64_t * filtered_hits_y, uint64_t * host_left_offset, uint64_t * host_right_offset, uint64_t id);
void filter_and_write_frags(uint32_t * filtered_hits_x, uint32_t * filtered_hits_y, uint32_t * host_left_offset, uint32_t * host_right_offset, uint32_t n_frags, FILE * out, char strand, uint32_t ref_len, uint32_t min_length);
uint32_t generate_hits_quadratic(uint32_t words_at_once, uint64_t * diagonals, Hit * hits, uint64_t * keys_x, uint64_t * keys_y, uint32_t * values_x, uint32_t * values_y, uint32_t items_x, uint32_t items_y, uint32_t query_len, uint32_t ref_len);
uint32_t generate_hits_fast(uint32_t max_hits, uint64_t * diagonals, Hit * hits, uint64_t * keys_x, uint64_t * keys_y, uint32_t * values_x, uint32_t * values_y, uint32_t items_x, uint32_t items_y, uint32_t query_len, uint32_t ref_len);
uint32_t generate_hits_sensitive(uint32_t max_hits, uint64_t * diagonals, Hit * hits, uint64_t * keys_x, uint64_t * keys_y, uint32_t * values_x, uint32_t * values_y, uint32_t items_x, uint32_t items_y, uint32_t query_len, uint32_t ref_len, uint32_t max_frequency);
uint32_t filter_hits_forward(uint64_t * diagonals, uint32_t * indexing_numbers, Hit * hits, uint32_t * filtered_hits_x, uint32_t * filtered_hits_y, uint32_t n_hits_found);
uint32_t filter_hits_reverse(uint64_t * diagonals, uint32_t * indexing_numbers, Hit * hits, uint32_t * filtered_hits_x, uint32_t * filtered_hits_y, uint32_t n_hits_found);
void read_kmers(uint64_t query_l, char * seq_x, uint64_t * keys_x, uint64_t * values_x);
void init_args(int argc, char ** av, FILE ** query, unsigned * selected_device, FILE ** ref, FILE ** out, uint32_t * min_length, int * fast, uint32_t * max_frequency, float * factor);
void perfect_hash_to_word(char * word, uint64_t hash, uint64_t k);
void print_kmers_to_file(uint64_t * keys, uint64_t * values, uint64_t table_size, FILE * fout);
char * get_dirname(char * path);
char * get_basename(char * path);
uint32_t load_seq(FILE * f, char * seq);
uint32_t get_seq_len(FILE * f);
void terror(const char * s);
void Qsort(uint64_t * keys, uint64_t * values, int64_t x, int64_t y);
