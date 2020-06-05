//////////////////////////////////////////////////////
//
//                  common.c
//
//	
//
//
//          Author(s): estebanpw, ortrelles
//
//////////////////////////////////////////////////////


#include "common.h"

uint32_t estimate_seq_size(FILE * f)
{

    if(f == NULL) {fprintf(stderr, "File is null. Check inputs.\n"); exit(-1);}
    
    fseek(f, 0L, SEEK_END);
    uint32_t l = (uint32_t) ftell(f);
    rewind(f);
    return l;

}

char * allocate_memory_for_sequence(FILE * f)
{

    uint32_t l = estimate_seq_size(f);
    if(l == 0) {fprintf(stderr, "Sequence is 0-sized. Check inputs.\n"); exit(-1);}
    char * s = NULL;
    s = (char *) malloc(l * sizeof(char));
    if(s == NULL) {fprintf(stderr, "Could not allocate memory for sequence.\n"); exit(-1);}

    return s;

}

char * load_seq(FILE * f, uint32_t * l, std::vector<uint64_t> * index)
{

    char c = '\0';
    *l = 0;
    char * s = allocate_memory_for_sequence(f);    

	// Skip first sequence
	while(c != '>')  c = getc(f);
	while(c != '\n') c = getc(f);


    while(!feof(f))
    {
        c = getc(f);
        if(c == '>')
        {
			// Add an 'N' between sequences to break kmers
			s[*l] = 'N';
			index->push_back(*l);
			*l = *l + 1;
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z')
        {
            s[*l] = c;
            *l = *l + 1;
        }
    }
    rewind(f);

    return s;

}

uint32_t get_seq_len(FILE * f) 
{

    char c = '\0';
    uint32_t l = 0;
    while(!feof(f))
    {
        c = getc(f);
        if(c == '>')
        {
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z') ++l;
    }
    rewind(f);

    return l;

}

char * reverse_complement_sequence(char * s, uint64_t l)
{

	char * r = NULL;
	r = (char *) malloc(l * sizeof(char));
	if(r == NULL) {fprintf(stderr, "Could not allocate memory for reverse sequence.\n"); exit(-1);}

	int64_t i, j = 0;
	for(i=(int64_t) l-1; i>=0; i--)
	{

		char c = s[i];
		if(c == 'A') r[j++] = 'T';
		else if(c == 'C') r[j++] = 'G';
		else if(c == 'G') r[j++] = 'C';
		else if(c == 'T') r[j++] = 'A';
		else r[j++] = 'N';

	}

	return r;

}

void get_alignments(char * s_x, char * s_y, char * r_y, uint64_t l_fastax, uint64_t l_fastay, std::vector<uint64_t> * index_x, std::vector<uint64_t> * index_y, std::vector<uint64_t> * index_r, FILE * csv)
{

	uint64_t i = 0;
	char buff[MAX_LINE];
	char * none, * up, * bottom;

	while(i < 17) { none = fgets(buff, MAX_LINE, csv); ++i; }

	uint64_t xstart, ystart, xend, yend, len;
	char strand;

	while(!feof(csv))
	{
		// Frag,443223,2424422,443383,2424582,f,0,160,128,128,80.00,80.00,0,0

		uint64_t idents = 0;
		if(fgets(buff, MAX_LINE, csv) == NULL) break;
		uint64_t t = strlen(buff);
		for(i=0; i<t; i++) if(buff[i] == ',') buff[i] = ' ';

		sscanf(buff, "%*s %lu %lu %lu %lu %c %*u %lu %*u %*u %*f %*f %*u %*u", &xstart, &ystart, &xend, &yend, &strand, &len);

		if(strand == 'f')
		{
			up = &s_x[xstart];
			bottom = &s_y[ystart];
		}
		else
		{
			ystart = l_fastay - ystart;
			--ystart;
			up = &s_x[xstart];
			bottom = &r_y[ystart];
		}

		fprintf(stdout, "%.*s\n", (int) len, up);
		for(i=0; i<len; i++)
		{
			if(up[i] == bottom[i])
			{
				fprintf(stdout, "|");
				++idents;
			}
			else 
				fprintf(stdout, " ");
		}

		fprintf(stdout, "\n%.*s\n", (int) len, bottom);
		uint64_t posX, posY, closest;
		if(index_x->size() > 0)	closest = search(xstart, index_x, index_x->size(), &posX); else posX = 0;
		if(index_y->size() > 0)
		{
			if(strand == 'f') closest = search(ystart, index_y, index_y->size(), &posY); else closest = search(ystart, index_r, index_r->size(), &posY);
		}else{
			posY = 0;
		}
		fprintf(stdout, "Query %lu Reference %lu HSP strand %c @ %lu %lu Idents %lu / %lu ( %.2f %%)\n", posX, posY, strand, xstart, ystart, idents, len, 100*(float)idents/(float)len);

		fprintf(stdout, "------------\n");

	}
	

}


uint64_t search(uint64_t value, std::vector<uint64_t> * a, uint64_t l, uint64_t * pos) {


	if(value < a->at(0)) {
		*pos = 0;
		return a->at(0);
	}
	if(value > a->at(l-1)) {
		*pos = l-1;
		return a->at(l-1);
	}

	int64_t lo = 0;
	int64_t hi = l - 1;


	while (lo <= hi) {
		int64_t mid = (hi + lo) / 2;

		if (value < a->at(mid)) {
			hi = mid - 1;
		} else if (value > a->at(mid)) {
			lo = mid + 1;
		} else {
			return a->at(mid);
		}
	}
	((int64_t) a->at(lo) - (int64_t) value) < ((int64_t)value - (int64_t) a->at(hi)) ? (*pos = (uint64_t) lo) : (*pos = (uint64_t) hi);
	return (uint64_t) ((int64_t) a->at(lo) - (int64_t) value) < ((int64_t)value - (int64_t) a->at(hi)) ? a->at(lo) : a->at(hi);
}






