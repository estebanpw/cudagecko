#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <arpa/inet.h>
#include "structs.h"
#include "commonFunctions.h"

#define MAXBUF 10000000

int readHashEntry(hashentry *h, FILE *f, uint64_t freqThr) {
	do {
		if(fread(h, sizeof(hashentry), 1, f)!=1){
			if(ferror(f))terror("Error reading hash entry from the file");
		}
	} while (!feof(f) && h->num > freqThr);

	if (feof(f))
		return -1;

	return h->num;
}

void loadWordOcurrences(hashentry he, location** pos, FILE** f) {
	// Load word position for he.word
	if (he.num > MAXBUF) {
		*pos = (location*) realloc(pos, he.num * sizeof(location));
		if (pos == NULL)
			terror("memory for position array");
	}
	fseek(*f, he.pos, SEEK_SET);
	if(fread(*pos, sizeof(location), he.num, *f)!=he.num){
		terror("Not possible to read the number of elements described");
	}
}

unsigned long scoreMax(char *seq, char *seq2, uint64_t len, int point) {
	//CHANGE WHEN USING PAM MATRIX
	//Scor=0; for (i=0;ii<len;i++) if (Seq[i]==Seq2[i])scor+=Ptos
	/* TO DELETE WARNINGS*/
	if(seq+1){
	}

	if(seq2+1){
	}
	/**********************/

	return len * point;
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
    if ((sX = calloc(1, sizeof(struct Sequence))) == NULL)
        terror("Memory...");

    fseek(f, 0, SEEK_END);
    SIZE = ftell(f);
    fseek(f, 0, SEEK_SET);

    if ((sX->datos = calloc(SIZE, sizeof(char))) == NULL) {
        terror("Memory for sequence...");
    }

    if ((seq = calloc(READBUF, sizeof(char))) == NULL) {
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

long getSeqLength(struct Sequence *s) {
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

void endianessConversion(char *source, char *target, int numberOfBytes){
	int i,j;
	for(i=numberOfBytes-1;i>=0;i--){
		j=numberOfBytes-1-i;
		target[j]=source[i];
	}
}


/**
 * Function to read a fragment from the specified file
 */
void readFragment(struct FragFile *frag, FILE *f) {
    char tmpArray[sizeof(long double)];

    if (htons(1) == 1) {
        //big endian
        if (fread(&frag->diag, sizeof(int64_t), 1, f) != 1) {
            if (feof(f))return;
            terror("Error reading the HSP diagonal");
        }
        if (fread(&frag->xStart, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xStart");
        }
        if (fread(&frag->yStart, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yStart");
        }
        if (fread(&frag->xEnd, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xEnd");
        }
        if (fread(&frag->yEnd, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yEnd");
        }
        if (fread(&frag->length, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP length");
        }
        if (fread(&frag->ident, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP identities");
        }
        if (fread(&frag->score, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP score");
        }
        if (fread(&frag->similarity, sizeof(float), 1, f) != 1) {
            terror("Error reading the HSP similarity");
        }
        if (fread(&frag->seqX, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqX");
        }
        if (fread(&frag->seqY, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqY");
        }
        if (fread(&frag->block, sizeof(int64_t), 1, f) != 1) {
            terror("Error reading the HSP block");
        }
        frag->strand = fgetc(f);
        if (fread(&frag->evalue, sizeof(long double), 1, f) != 1) {
            terror("Error reading the HSP evalue");
        }
    } else {
        //little endian
        if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
            if (feof(f))return;
            terror("Error reading the HSP diagonal");
        }
        endianessConversion(tmpArray, (char *) (&frag->diag), sizeof(int64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xStart");
        }
        endianessConversion(tmpArray, (char *) (&frag->xStart), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yStart");
        }
        endianessConversion(tmpArray, (char *) (&frag->yStart), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xEnd");
        }
        endianessConversion(tmpArray, (char *) (&frag->xEnd), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yEnd");
        }
        endianessConversion(tmpArray, (char *) (&frag->yEnd), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP length");
        }
        endianessConversion(tmpArray, (char *) (&frag->length), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP identity");
        }
        endianessConversion(tmpArray, (char *) (&frag->ident), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP score");
        }
        endianessConversion(tmpArray, (char *) (&frag->score), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(float), 1, f) != 1) {
            terror("Error reading the HSP float");
        }
        endianessConversion(tmpArray, (char *) (&frag->similarity), sizeof(float));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqX");
        }
        endianessConversion(tmpArray, (char *) (&frag->seqX), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqY");
        }
        endianessConversion(tmpArray, (char *) (&frag->seqY), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
            terror("Error reading the HSP block");
        }
        endianessConversion(tmpArray, (char *) (&frag->block), sizeof(int64_t));
        frag->strand = fgetc(f);
        if (fread(tmpArray, sizeof(long double), 1, f) != 1) {
            terror("Error reading the HSP evalue");
        }
        endianessConversion(tmpArray, (char *) (&frag->evalue), sizeof(long double));
    }
}

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f) {
    char tmpArray[sizeof(long double)];
    if (htons(1) == 1) {
        //Big endian
        fwrite(&frag->diag, sizeof(int64_t), 1, f);
        fwrite(&frag->xStart, sizeof(uint64_t), 1, f);
        fwrite(&frag->yStart, sizeof(uint64_t), 1, f);
        fwrite(&frag->xEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag->yEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag->length, sizeof(uint64_t), 1, f);
        fwrite(&frag->ident, sizeof(uint64_t), 1, f);
        fwrite(&frag->score, sizeof(uint64_t), 1, f);
        fwrite(&frag->similarity, sizeof(float), 1, f);
        fwrite(&frag->seqX, sizeof(uint64_t), 1, f);
        fwrite(&frag->seqY, sizeof(uint64_t), 1, f);
        fwrite(&frag->block, sizeof(int64_t), 1, f);
        fputc(frag->strand, f);
        fwrite(&frag->evalue, sizeof(long double), 1, f);
    } else {
        //Little endian
        endianessConversion((char *) (&frag->diag), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        endianessConversion((char *) (&frag->xStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->yStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->xEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->yEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->length), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->ident), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->score), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->similarity), tmpArray, sizeof(float));
        fwrite(tmpArray, sizeof(float), 1, f);
        endianessConversion((char *) (&frag->seqX), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->seqY), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->block), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        fputc(frag->strand, f);
        endianessConversion((char *) (&frag->evalue), tmpArray, sizeof(long double));
        fwrite(tmpArray, sizeof(long double), 1, f);
    }
}

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//big endian
		if(fread(length, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading sequence length");
		}
	} else {
		//little endian
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading sequence length");
		}
		endianessConversion(tmpArray, (char *)length, sizeof(uint64_t));
	}
}

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//big endian
		fwrite(length, sizeof(uint64_t), 1, f);
	} else {
		//little endian
		endianessConversion((char *)length, tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
	}
}

/**
 * Function to return the sizeof a fragment.
 * Due to architeture issues the value of sizeof built-in
 * function could be different
 */
long int sizeofFragment(){
	return 9*sizeof(uint64_t)+2*sizeof(int64_t)+1*sizeof(float)+1*sizeof(char)+sizeof(long double);
}

/**************** ARJONA *******************/
void cpyFrag2(struct FragFile *f, struct FragFile g){

	f->diag = g.diag;
	f->xStart = g.xStart;
	f->xEnd = g.xEnd;
	f->yStart = g.yStart;
	f->yEnd = g.yEnd;
	f->length = g.length;
	f->score = g.score;
	f->similarity = g.similarity;
	f->seqX = g.seqX;
	f->seqY = g.seqY;
	f->ident = g.ident;
	f->block = g.block;
	f->strand = g.strand;
	f->evalue = g.evalue;
}

/******/
struct FragFile* readFragments(char* s,int* nf,uint64_t *xtotal,uint64_t *ytotal){
//Fragment* readFragments(char* s,int* nf,int *xtotal,int *ytotal){

	FILE* fe;

	struct FragFile* fs,f;
	int n;

	if((fe=fopen(s,"rb"))==NULL){
		printf("***ERROR Opening input file");
		exit(-1);
	}
	n=0;
	readSequenceLength(xtotal, fe);
	readSequenceLength(ytotal, fe);
//	long int longFile;
//	fseek(fe,0,SEEK_END);
//	longFile=ftell(fe);
//	n=(int)(longFile-2*sizeof(uint64_t))/sizeof(struct FragFile);
//	printf("\n ReadFragments Complete\nnum: %d\n",n);

	// Alternativa +++++++++++++++++
	n=0;
	readFragment(&f, fe);
	while(!feof(fe)){
		readFragment(&f, fe);
		n++;
	}
//	printf("\n ReadFragments Complete\nnum: %d\n",n);
	//+++++++++++++
	rewind(fe);
	fs=(struct FragFile*)malloc(sizeof(struct FragFile)*(n));
	if(fs==NULL){printf("****ERROR: Out of memory\n");exit(-1);}

	*nf=n;
	readSequenceLength(xtotal, fe);
	readSequenceLength(ytotal, fe);
	n=0;
	readFragment(&f, fe);
	while(!feof(fe)){



		if(f.length>0){
			cpyFrag2(&fs[n],f);
		}
		n++;
		readFragment(&f, fe);
		//fprintf(stdout,"%d\t%" PRId64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%c" "\t%" PRIu64 "\n",n,fs[n].xStart, fs[n].yStart, fs[n].xEnd, fs[n].yEnd, fs[n].length, fs[n].strand, fs[n].ident);

	}
	*nf=n;


	fclose(fe);
	return fs;
}
/************************/

