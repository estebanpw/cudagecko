#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"

void terror(const char * s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

/**
 * Function to read the sequence length
 */
void read_sequence_length(uint64_t *length, FILE *f) {
    
    if (fread(length, sizeof(uint64_t), 1, f) != 1) {
        terror("Error reading sequence length");
    }
}

/**
 * Function to write the sequence length
 */
void write_sequence_length(uint64_t *length, FILE *f) {
    
    fwrite(length, sizeof(uint64_t), 1, f);
    
}

unsigned long timestart() {
    struct timeval tv;

    gettimeofday(&tv, NULL);

    return (tv.tv_usec / 1000) + (tv.tv_sec * 1000);
}

unsigned long timestop(unsigned long start) {
    struct timeval tv;

    gettimeofday(&tv, NULL);

    return (tv.tv_usec / 1000) + (tv.tv_sec * 1000) - start;
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

struct Sequence *LeeSeqDB(FILE *f, uint64_t *n, uint64_t *nSeqs,
                          int fAst) {
    char c;
    uint64_t lon = 0, k = 0, seqs = 0;
    uint64_t lonFinal = 0;
    uint64_t SIZE = 0;
    uint64_t i = 0, r = 0;
    struct Sequence *sX;
    char *seq = NULL;

    //Initialize
    *n = 0;

    //Memory
    if ((sX = (struct Sequence *) calloc(1, sizeof(struct Sequence))) == NULL)
        terror("Memory...");

    fseek(f, 0, SEEK_END);
    SIZE = ftell(f);
    fseek(f, 0, SEEK_SET);

    if ((sX->datos = (char *) calloc(SIZE, sizeof(char))) == NULL) {
        terror("Memory for sequence...");
    }

    if ((seq = (char *) calloc(READBUF, sizeof(char))) == NULL) {
        terror("Memory for sequence...");
    }

    i = READBUF + 1;

    if (!fAst)
        while ((c = buffered_fgetc(seq, &i, &r, f)) != '>' && (!feof(f) || (feof(f) &&  i < r ) )); //start seq
    if (feof(f) && i >= r)
        return 0;
    while ((c = buffered_fgetc(seq, &i, &r, f)) == ' ');

    while (k < MAXLID && c != '\n' && c != ' ') {
        if (feof(f) && i >= r)
            return 0;

        sX->ident[k++] = c;
        c = buffered_fgetc(seq, &i, &r, f);
    }

    sX->ident[k] = 0; //end of data.
    while (c != '\n')
        c = buffered_fgetc(seq, &i, &r, f);
    c = buffered_fgetc(seq, &i, &r, f);

    //start list with sX2
    while (!feof(f) || (feof(f) &&  i < r ) ) {
        c = toupper(c);
        if (c == '>') {
            fAst = 1;
            seqs++;
            sX->datos[lon++] = '*';
            while (c != '\n') {
                if ((feof(f) &&  i >= r ))
                    return 0;
                c = buffered_fgetc(seq, &i, &r, f);
            }
            //break;
        }
        if (isupper(c))
            sX->datos[lon++] = c;
        if (c == '*') {
            sX->datos[lon++] = c;
        }
        c = buffered_fgetc(seq, &i, &r, f);
    }

    free(seq);

    sX->datos[lon] = 0x00;

    lonFinal += lon;
    *nSeqs = seqs + 1;
    *n = lonFinal - seqs;
    return sX;
}

char getValue(struct Sequence *s, uint64_t pos) {
    return s->datos[pos];
}

long getSeqLength(struct Sequence *s, uint64_t start, int ns) {
    uint64_t s1 = 0;
    while (s1 > 0 && s->datos[s1] != '*') {
        s1--;
    }
    s1++;
    char *tmp = strchr(s->datos + s1, '*');
    if (tmp == NULL) {
        return strlen(s->datos) - s1 + 1;
    }
    return tmp - (s->datos + s1) + 1;
}

uint64_t quick_pow4(uint32_t n){
    static uint64_t pow4[33]={1L, 4L, 16L, 64L, 256L, 1024L, 4096L, 16384L, 65536L, 
    262144L, 1048576L, 4194304L, 16777216L, 67108864L, 268435456L, 1073741824L, 4294967296L, 
    17179869184L, 68719476736L, 274877906944L, 1099511627776L, 4398046511104L, 17592186044416L, 
    70368744177664L, 281474976710656L, 1125899906842624L, 4503599627370496L, 18014398509481984L, 
    72057594037927936L, 288230376151711744L, 1152921504606846976L, 4611686018427387904L};
    return pow4[n];
}

uint64_t quick_pow4byLetter(uint32_t n, const char c){
    static uint64_t pow4_G[33]={2*1L, 2*4L, 2*16L, 2*64L, 2*256L, 2*1024L, 2*4096L, 2*16384L, 2*65536L, 
    (uint64_t)2*262144L, (uint64_t)2*1048576L,(uint64_t)2*4194304L, (uint64_t)2*16777216L, (uint64_t)2*67108864L, (uint64_t)2*268435456L, (uint64_t)2*1073741824L, (uint64_t)2*4294967296L, 
    (uint64_t)2*17179869184L, (uint64_t)2*68719476736L, (uint64_t)2*274877906944L, (uint64_t)2*1099511627776L, (uint64_t)2*4398046511104L, (uint64_t)2*17592186044416L, 
    (uint64_t)2*70368744177664L, (uint64_t)2*281474976710656L, (uint64_t)2*1125899906842624L, (uint64_t)2*4503599627370496L, (uint64_t)2*18014398509481984L, 
    (uint64_t)2*72057594037927936L, (uint64_t) 2*288230376151711744L, (uint64_t) 2*1152921504606846976L, (uint64_t) 2*4611686018427387904L};
    
    static uint64_t pow4_T[33]={3*1L, 3*4L, 3*16L, 3*64L, 3*256L, 3*1024L, 3*4096L, 3*16384L, 3*65536L, 
    (uint64_t)3*262144L, (uint64_t) 3*1048576L, (uint64_t)3*4194304L, (uint64_t)3*16777216L, (uint64_t)3*67108864L, (uint64_t)3*268435456L, (uint64_t)3*1073741824L, (uint64_t)3*4294967296L, 
    (uint64_t)3*17179869184L, (uint64_t)3*68719476736L, (uint64_t)3*274877906944L, (uint64_t)3*1099511627776L, (uint64_t)3*4398046511104L, (uint64_t)3*17592186044416L, 
    (uint64_t)3*70368744177664L, (uint64_t)3*281474976710656L, (uint64_t)3*1125899906842624L, (uint64_t)3*4503599627370496L, (uint64_t)3*18014398509481984L, 
    (uint64_t)3*72057594037927936L, (uint64_t) 3*288230376151711744L, (uint64_t) 3*1152921504606846976L, (uint64_t) 3*4611686018427387904L};
    
    if(c == 'A') return 0;
    if(c == 'C') return quick_pow4(n);
    if(c == 'G') return pow4_G[n];
    if(c == 'T') return pow4_T[n];
    return 0;
}

uint64_t hashOfWord(const char * word, uint32_t k){
    /*
    uint64_t value = 0, jIdx;
    for(jIdx=0;jIdx<k;jIdx++){
        value += quick_pow4byLetter(k-(jIdx+1), word[jIdx]);
    }
    return value;
    */
    if(word[0]=='A') return 0;
    if(word[1]=='C') return 1;
    if(word[2]=='G') return 2;
    return 3;
}

long double chiSquaredAlfaTest(struct tFreqs of){

    long double chi2value[11]={0.07172177,  0.21579528,  4.64162768,  6.25138863,  7.81472790,
    9.34840360,  9.83740931, 11.34486673, 12.83815647, 14.79551705, 16.26623620};
    long double pvalue[11]={0.95, 0.90, 0.80, 0.70, 0.50, 0.30, 0.20, 0.10, 0.05, 0.01, 0.001};
    long double chi2 = 0;

    //Chi squared test
    
    long double Aef,Cef,Gef,Tef;
    long double Aof,Cof,Gof,Tof;
    Aef = Cef = Gef = Tef = 0.25;
    Aof = (long double)of.A/of.total;
    Cof = (long double)of.C/of.total;
    Gof = (long double)of.G/of.total;
    Tof = (long double)of.T/of.total;

    //Differences between obtained and expected values
    // FAKE
    //Aof = 0.25; Cof = 0.25; Gof = 0.25; Tof = 0.25;


    chi2 += ((Aof-Aef)*(Aof-Aef))/(Aef);
    chi2 += ((Cof-Cef)*(Cof-Cef))/(Cef);
    chi2 += ((Gof-Gef)*(Gof-Gef))/(Gef);
    chi2 += ((Tof-Tef)*(Tof-Tef))/(Tef);
    chi2 *= of.total;

    

    //Find best pvalue
    int i;
    long double confidenceLevelCHI = 0.001;
    for(i=10;i>=0;i--){
        if(chi2 <= chi2value[i]) confidenceLevelCHI = pvalue[i];
    }

    fprintf(stdout, "[INFO] CHI-squared test using Pearson's: %Le :->: Using p-value: %Le\n", chi2, confidenceLevelCHI);
    return confidenceLevelCHI;
}

void computeFrequencies(FILE * db, struct tFreqs * tf){
    char c;
    //Count absolute number of letter frequencies
    while(!feof(db)){
        c = fgetc(db);
        if(c == '>') while(c != '\n') c = fgetc(db);
        switch(c){
            case 'A': tf->A++;
            break;
            case 'C': tf->C++;
            break;
            case 'G': tf->G++;
            break;
            case 'T': tf->T++;
            break;
            default:
            break;
        }
    }
    tf->total = tf->A + tf->C + tf->G + tf->T;
    //Rewind back to start
    fseeko64(db, 0L, SEEK_SET);
}

uint64_t readLengthOfSequences(FILE * f, uint64_t * seqsLen, uint64_t * nSeqsCounted){

    //Number of sequences counted
    if(nSeqsCounted != NULL) *nSeqsCounted = 0;

    uint64_t i = 0;
    uint64_t t_length = 0;
    char c = '0'; //Some character

    while(!feof(f)){
        
        if(c == '>'){
            while(c != '\n') c = fgetc(f); //Skip id

            while(c != '>' && !feof(f)){ //Until next id
                c = fgetc(f);
                if((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')){
                    if(seqsLen != NULL) seqsLen[i]++;
                    t_length++;
                }
            }

            //printf("SEQ LEN %"PRIu64": %"PRIu64"\n", i, seqsLen[i]);
            if(nSeqsCounted != NULL)*nSeqsCounted += 1;
            i++; //Next sequence
        }else{
            c = fgetc(f);    
        }
        
    }

    //Rewind
    fseeko64(f, 0L, SEEK_SET);
    return t_length;
}

uint64_t getNumberOfSequences(FILE * f){
    uint64_t t_seqs = 0;
    char c = 'Z';
    while(!feof(f)){
        if(c=='>') t_seqs++;
        c=fgetc(f);
    }
    //Rewind
    fseeko64(f, 0L, SEEK_SET);
    return t_seqs;
}


inline void strrev(char *p, char *d, uint32_t k){
    int i;
    char c;
    for(i=0;i<(int)k;i++){
	c = p[k-i-1];
	switch(c){
	case 'A': c='T';
	break;
	case 'C': c='G';
	break;
	case 'G': c='C';
	break;
	case 'T': c='A';
	break;
	}
        d[i] = c; 
    }
  
}

// Returns the probability of x, given the distribution described by mu and sigma.
long double pdf(long double x, long double mu, long double sigma){
        return expl( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
}

// Computes the (approximated) cumulative distribution function in the interval [-inf,x] of a gaussian distribution
long double cdf(long double x, long double mu, long double sigma){
        return 0.5 * (1 + erfl((x - mu) / (sigma * sqrtl(2.))));
}

