#ifndef QUICKSORT_H
#define QUICKSORT_H

/**
 * Function to order the ifile using a buffer size of
 * 'maxsize' and a thread number of 'nproc'. the result
 * will be written to ofile
 */
int psort(int maxsize, int nproc, char* ifile, char* ofile);

/**
 * Function to determine if object 1 is strictly greater than 2.
 * Returns 1 if true, 0 otherwise
 * Has to be defined in the main procedure
 */
int GT(BaseType a1, BaseType a2);

#endif /* QUICKSORT_H */
