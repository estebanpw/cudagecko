#ifndef COMPARISON_FUNCTIONS_H
#define COMPARISON_FUNCTIONS_H

#define min(x,y)    (((x) < (y)) ? (x) : (y))

/**
 * Read a new hash entry from 'f' with no more occurences
 * than 'freqThr'
 */
int readHashEntry(hashentry *h, FILE *f);

/**
 * Read the ocurrences of the given hash entry from the
 * ocurrences file 'f' and stored them in 'pos'
 */
void loadWordOcurrences(hashentry he, location** pos, FILE** f);

/**
 * Get the maximum score possible between the two sequences
 */
unsigned long scoreMax(char *seq, char *seq2, uint64_t len, int point);

/**
 * Read the sequence from disk and store it in a list of Sequence struct
 * n: sequence length
 * ns: number of nodes in the list
 */
struct Sequence* LeeSeqDB(FILE *f, uint64_t *n, uint64_t *nSeqs, int fAst);

/**
 * Get the value of the sequence in a given position of the list node ns
 */
char getValue(struct Sequence *s, uint64_t pos);

/**
 * Get the length of the sequence 's'
 */
long getSeqLength(struct Sequence *s);

/**
 * Function to read a fragment from the specified file
 */
void readFragment(struct FragFile *frag, FILE *f);

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f);

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f);

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f);

/**
 * Function to return the sizeof a fragment.
 * Due to architeture issues the value of sizeof built-in
 * function could be different
 */
long int sizeofFragment();

#endif /* COMPARISON_FUNCTIONS_H */
