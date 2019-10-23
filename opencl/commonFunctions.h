#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H
#include "structs.h"
/**
 * Print the error message 's' and exit(-1)
 */
void terror(char *s);

void read_sequence_length(uint64_t *length, FILE *f);

void write_sequence_length(uint64_t *length, FILE *f);

unsigned long timestart();

unsigned long timestop(unsigned long start);

/**
 * Function to read char by char buffered from a FILE
 */
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

/**
 * Read the sequence from disk and store it in a list of Sequence struct
 * n: sequence length
 * ns: number of nodes in the list
 */
struct Sequence*LeeSeqDB(FILE *f, uint64_t *n, uint64_t *nSeqs, int fAst);

/**
 * Get the value of the sequence in a given position of the list node ns
 */
char getValue(struct Sequence *s, uint64_t pos);

/**
 * Get the length of the sequence 's'
 */
long getSeqLength(struct Sequence *s, uint64_t start, int ns);

/**
 *	Compute the power of n to for with a lookup table
 */
uint64_t quick_pow4(uint32_t n);

/*
 *	Compute the power of 0,1,2,3 to the power of n quickly by lookup table
*/
uint64_t quick_pow4byLetter(uint32_t n, const char c);

/**
 *	Compute the hash of a word of letters of length k (in ascii)
 */
uint64_t hashOfWord(const char * word, uint32_t k);

/*
	Computes the chi squared test yielding the best pvalue possible (out of a table)
*/
long double chiSquaredAlfaTest(struct tFreqs of);

/*
	Compute the frequencies of the residues
*/
void computeFrequencies(FILE * db, struct tFreqs * tf);

/*
	Read each sequence's length and store them in an array. Return total length (sum).
	If seqsLen is null then individual sequences length will not be stored (to only calculate total length)
	nSeqsCounted will hold how many sequences were counted (should be equal to length of seqsLen)
*/
uint64_t readLengthOfSequences(FILE * f, uint64_t * seqsLen, uint64_t * nSeqsCouted);

/*
	Counts the number of '>' in a FASTA file
*/
uint64_t getNumberOfSequences(FILE * f);

/*
	Reverses a string p in d
*/
inline void strrev(char *p, char *d, uint32_t k);

/* 
	Returns the probability of x, given the distribution described by mu and sigma.
*/
long double pdf(long double x, long double mu, long double sigma);

/*
	Computes the (approximated) cumulative distribution function in the interval [-inf,x] of a gaussian distribution
*/
long double cdf(long double x, long double mu, long double sigma);



#endif /* COMMON_FUNCTIONS_H */
