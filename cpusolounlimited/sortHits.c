#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "commonFunctions.h"
#include "quicksort.h"

int GT(BaseType a1, BaseType a2) {
	if(a1.diag > a2.diag)
		return 1;
	else if (a1.diag < a2.diag)
		return 0;
	if (a1.posX > a2.posX)
		return 1;

	return 0;
}

int main(int ac, char** av) {
	if (ac < 4) {
		printf("USE: sortHits <max_size> <num_proc> <input_file> <output_file>\n");
		exit(1);
	}

	return psort(atoi(av[1]), atoi(av[2]), av[3], av[4]);
}

