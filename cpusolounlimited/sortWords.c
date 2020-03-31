#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "commonFunctions.h"
#include "quicksort.h"

int GT(BaseType a1, BaseType a2){
	if(a1.w.b[0] > a2.w.b[0]) return 1;
	else if(a1.w.b[0] < a2.w.b[0]) return 0;

	if(a1.w.b[1] > a2.w.b[1]) return 1;
	else if(a1.w.b[1] < a2.w.b[1]) return 0;

	if(a1.w.b[2] > a2.w.b[2]) return 1;
	else if(a1.w.b[2] < a2.w.b[2]) return 0;

	if(a1.w.b[3] > a2.w.b[3]) return 1;
	else if(a1.w.b[3] < a2.w.b[3]) return 0;

	if(a1.w.b[4] > a2.w.b[4]) return 1;
	else if(a1.w.b[4] < a2.w.b[4]) return 0;

	if(a1.w.b[5] > a2.w.b[5]) return 1;
	else if(a1.w.b[5] < a2.w.b[5]) return 0;

	if(a1.w.b[6] > a2.w.b[6]) return 1;
	else if(a1.w.b[6] < a2.w.b[6]) return 0;

	if(a1.w.b[7] > a2.w.b[7]) return 1;
	else if(a1.w.b[7] < a2.w.b[7]) return 0;

	if(a1.seq > a2.seq) return 1;
	else if(a1.seq < a2.seq) return 0;

	if(a1.pos > a2.pos) return 1;
	return 0;
}

int main(int ac, char** av){
	if(ac<4) {
		terror("USE: sortWords <max_size> <num_proc> <input_file> <output_file>\n");
	}

	return psort(atoi(av[1]),atoi(av[2]),av[3],av[4]);
}

