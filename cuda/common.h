//////////////////////////////////////////////////////
//
//                  common.h
//
//      Definitions of structures and methods
//
//
//          Author(s): estebanpw, ortrelles
//
//////////////////////////////////////////////////////



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <vector>
#include <iostream>

#define KSIZE 32                // Used for word size in dictionaries
#define MAX_NAME 1024           // Used for names of automatically generated files
#define MAX_LINE 1024           // Used for maximum length per break line in text files
#define BATCH_SIZE 1024 * 1024  // Used for progressively loading words into memory


uint32_t estimate_seq_size(FILE * f);

char * allocate_memory_for_sequence(FILE * f);

char * load_seq(FILE * f, uint32_t * l, std::vector<uint64_t> * index);

uint32_t get_seq_len(FILE * f);

char * reverse_complement_sequence(char * s, uint64_t l);

void get_alignments(char * s_x, char * s_y, char * r_y, uint64_t l_fastax, uint64_t l_fastay, std::vector<uint64_t> * index_x, std::vector<uint64_t> * index_y, std::vector<uint64_t> * index_r, FILE * csv);

uint64_t search(uint64_t value, std::vector<uint64_t> * a, uint64_t l, uint64_t * pos);
