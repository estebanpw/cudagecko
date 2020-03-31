#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>
//Structs required for the dotplot workflow
#define MAXLID 200
#define MAXLS 1000000000
#define READBUF 50000000 //50MB

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
} wentry;

//Struct for w2hd program
typedef struct {
	//Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint64_t num;
} hashentry;

//Struct for w2hd program
typedef struct {
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} location;

//Struct for hits, sortHits and filterHits programs
typedef struct {
	//Diagonal where the hit is located
	//This value is calculated as:
	//posX - posY
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
};

//Struct for leeSeqDB function
struct Sequence{
    char ident[MAXLID+1];
    char *datos;
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

#endif
