#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>
#include <CL/cl.h>

//Structs required for the dotplot workflow
#define MAXLID 200
//#define READBUF 2000000000 //2 GB
#define READBUF 50000000 //50MB
#define INITSEQS 300000 //Number of starting sequences (in database)
#define REALLOC_FREQ 40000
#define POINT 4
#define BYTES_IN_WORD 8
#define MAX_DEVICES 32
#define MAX_FRAMES 10000


#define MAXLID 200

//Struct for words program
typedef struct {
    //Each letter is stored using 2 bits
    //We have 4 letters per byte and a
    //maximum of 32 in 'b'
    unsigned char b[8];
} word;


//Struct for words and sort program
typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
    //Using this field to store the sme word in reverse
    //in the same iteration. Either 'f' or 'r'
    char strand;
} wentryR;

//Struct for w2hd program
typedef struct {
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
    //Using this field to store the sme word in reverse
    //in the same iteration. Either 'f' or 'r'
    char strand;
} locationR;

// Parameters for GPU words
typedef struct param_word{
    ulong kmer_size;
	ulong seq_length;
	ulong t_work_items;
	ulong kmers_per_work_item;
    ulong offset;
} parameter_word;

//Struct for GPU words
typedef struct {
    unsigned char b[8];     // The actual word
    uint64_t pos;           // The position in the sequence file
                            // The sequence it belongs to will come later
} wordGPU;

typedef struct {
    uint64_t pos_x;
    uint64_t pos_y;
    uint64_t seq_x;
    uint64_t seq_y;
} hitGPU;

typedef struct {
    ulong pos_x;
    ulong pos_y;
    ulong length;
    ulong identities;
    // Notice that seq_x and seq_y can be retrieved from the hits when CPU is writing them
} reduced_frag;

typedef struct {
    ulong kmer_size;
    ulong size_x;
    ulong size_y;
    ulong t_hits; // Tell how many hits in total there will be as opposed to work items which are divisible by local item size
} parameter_frags;

//Struct for w2hd program
typedef struct {
    //Word compressed in binary format
    word w;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint64_t num;
    //The ocurrences with position and
    //sequence
    locationR *locs;
} hashentryR;

//Struct for words and sort program
typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} wentryF;

//Struct for w2hd program
typedef struct {
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} locationF;

//Struct for w2hd program
typedef struct {
    //Word compressed in binary format
    word w;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint64_t num;
    //The ocurrences with position and
    //sequence
    locationF *locs;
} hashentryF;

//Struct for hits, sortHits and filterHits programs
typedef struct {
    int64_t diag;
    //Ocurrence position in sequence X
    uint64_t posX;
    //Ocurrence position in sequence Y
    uint64_t posY;
    //For multiple sequence files this var
    //reflects in what sequence of X file
    //occurs the word
    uint64_t seqX;
    //For multiple sequence files this var
    //reflects in what sequence of Y file
    //occurs the word
    uint64_t seqY;
} hit;



//Struct for FragHits, af2png and leeFrag programs
struct FragFile {
    //Diagonal where the frag is located
    //This value is calculated as:
    //posX - posY
    int64_t diag;
    //Start position in sequence X
    uint64_t xStart;
    //Start position in Sequence Y
    uint64_t yStart;
    //End position in Sequence X
    uint64_t xEnd;
    //End position in Sequence Y
    uint64_t yEnd;
    //Fragment Length
    //For ungaped aligment is:
    //xEnd-xStart+1
    uint64_t length;
    //Number of identities in the
    //fragment
    uint64_t ident;
    //Score of the fragment. This
    //depends on the score matrix
    //used
    uint64_t score;
    //Percentage of similarity. This
    //is calculated as score/scoreMax
    //Where score max is the maximum
    //score possible
    float similarity;
    //sequence number in the 'X' file
    uint64_t seqX;
    //sequence number in the 'Y' file
    uint64_t seqY;
    //synteny block id
    int64_t block;
    //'f' for the forward strain and 'r' for the reverse
    char strand;
    //E-value of fragment
    long double evalue;
};

//Struct for leeSeqDB function
struct Sequence {
    char ident[MAXLID + 1];
    char *datos;
};

//Struct for holding nucleotide frequencies
struct tFreqs{
    uint64_t A;
    uint64_t C;
    uint64_t G;
    uint64_t T;
    uint64_t total;
};

//Struct for calculating karlin and lambda parameters for different query/sequence and PAM matrix
struct statsHSP{
    struct tFreqs tf;
    double karlin;
    double lambda;
    
};

//Struct for reads index tuple
struct rIndex2 {
    char id[MAXLID];
    uint64_t  rNumber;
    uint64_t  rLen;
    uint64_t  rLmasked; //Masked positions
    uint64_t  nonACGT;  //N's
    uint64_t  pos;      //Start position of sequence
    uint64_t  Lac;      //Accumulated length
};

// For the next generation, GPU-like gecko

//Struct for GPU words
typedef struct {
    uint64_t h;        // The actual word
    uint64_t pos;      // The position in the sequence file
} word_hash_GPU;

typedef struct {
    ulong offset;
    ulong end;
    ulong kmer_size;
    ulong kmers_in_work_item;
    ulong t_work_items;
} param_words_advanced;

typedef struct param_sort{
    ulong N;
    ulong comparators_per_wi;
    ulong stage;
    ulong step;
    ulong kmer_size;
} parameter_sort;

typedef struct {
    uint64_t pos_x;
    uint64_t pos_y;
} hit_advanced;


#endif
