#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"

int seq2word(char* buf, int wsize, word* w) {
	int i;
	int b = 6;
	memset(w, 0, sizeof(word));

	for (i = 0; i < wsize; i++) {
		if (buf[i] >= 4)
			return -1;
		w->b[i / 4] |= buf[i] << b;
		b -= 2;
		if (b < 0)
			b = 6;
	}
	return 0;
}

void skipIDLine(FILE *fIn) {
	char c;
	// first line (skip ID fasta Line)
	c = fgetc(fIn);
	while (c != '\n')
		c = fgetc(fIn);
}

int letterToIndex(char c) {
	// coding (a=0,c=1,g=2,t=3,'>'=4 others=9 )
	switch (c) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	case '>':
		return 4;
	default:
		return 9;
	}
}

int wordcmp(unsigned char *w1, unsigned char *w2, int n) {

	int i = 0, limit;

	if(n%4 != 0){
		w1[n/4] = w1[n/4] >> (2*(3-((n-1)%4)));
		w1[n/4] = w1[n/4] << (2*(3-((n-1)%4)));
		w2[n/4] = w2[n/4] >> (2*(3-((n-1)%4)));
		w2[n/4] = w2[n/4] << (2*(3-((n-1)%4)));
		limit=(n/4)+1;
	} else {
		limit = n/4;
	}

	for (i=0;i<limit;i++) {
		if (w1[i]<w2[i]) return -1;
		if (w1[i]>w2[i]) return +1;
	}
	return 0;
}

void showWord(word* w, char *ws) {
	char Alf[] = { 'A', 'C', 'G', 'T' };
	int i;
	int wsize = 8;
	unsigned char c;
	for (i = 0; i < wsize; i++) {
		c = w->b[i];
		c = c >> 6;
		ws[4*i] = Alf[(int) c];
		c = w->b[i];
		c = c << 2;
		c = c >> 6;
		ws[4*i+1] = Alf[(int) c];
		c = w->b[i];
		c = c << 4;
		c = c >> 6;
		ws[4*i+2] = Alf[(int) c];
		c = w->b[i];
		c = c << 6;
		c = c >> 6;
		ws[4*i+3] = Alf[(int) c];
	}
}
