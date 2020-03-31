/*
 * fragmentv3.h
 *
 *  Created on: 23/02/2014
 *      Author: jarjonamedina
 */

#ifndef FRAGMENTV3_H_
#define FRAGMENTV3_H_
/*
 typedef struct {
    unsigned long xIni;
    unsigned long yIni;
    unsigned long xFin;
    unsigned long yFin;
    unsigned long diag;
    unsigned long length;
    unsigned long ident;
    unsigned long score;
    //float similarity;
    unsigned long seqX; //sequence number in the 'X' file
    unsigned long seqY; //sequence number in the 'Y' file
    unsigned long block;            //synteny block id
    char strand;        //'f' for the forward strain and 'r' for the reverse
}Fragment;

*/


/* The new one*/

typedef struct {
    unsigned long diag;
    unsigned long xIni;
    unsigned long yIni;
    unsigned long xFin;
    unsigned long yFin;
    unsigned long length;
    unsigned long ident;
    unsigned long score;
    float similarity;
    unsigned long seqX; //sequence number in the 'X' file
    unsigned long seqY; //sequence number in the 'Y' file
    int block;          //synteny block id
    char strand;        //'f' for the forward strain and 'r' for the reverse
		// New fields
		long double pvalue;
		float COVLocusX;
		float COVLocusY;
		// Modelling E.E.
		int id;
		int gap;
		int translocations;
		int reversions;
		int duplications;
		int events;
}Fragmentv3;

Fragmentv3* readFragmentsv3(char* s,int* nf,int *xtotal,int *ytotal);
int writeFragmentsv3(Fragmentv3* f,char* s,int nf, int xtotal, int ytotal);


#endif /* FRAGMENTV3_H_ */
