/*
 * fragmentv2.h
 *
 *  Created on: 17/07/2013
 *      Author: jarjonamedina
 */

#ifndef FRAGMENTV2_H_
#define FRAGMENTV2_H_
#include <inttypes.h>
#include "structs.h"


struct FragFile* readFragmentsv2(char* s,int* nf,uint64_t *xtotal,uint64_t *ytotal);

int writeFragments(struct FragFile* f,char* s,int nf,uint64_t  xtotal, uint64_t ytotal);


#endif /* FRAGMENTV2_H_ */
