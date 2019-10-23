#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

/**
 * Print the error message 's' and exit(-1)
 */
void terror(char *s);

/**
 * Function to buffer file reading
 */
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

#endif /* COMMON_FUNCTIONS_H */
