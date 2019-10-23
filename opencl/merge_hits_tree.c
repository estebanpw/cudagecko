
// Standard utilities and common systems includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "structs.h"


#define BUFFER_SIZE 2048
#define HITS_BUFFER 400000  // times 40 and for 300 frames aprox it is 4.8 GB
#define MAX_KERNEL_SIZE BUFFER_SIZE*100
#define CORES_PER_COMPUTE_UNIT 32
#define MIN(a,b) (((a)<(b))?(a):(b))

char * get_dirname(char * path);
char * get_basename(char * path);
void init_args(int argc, char ** av, FILE ** hits_file, FILE ** hits_file_headers, FILE ** hits_filtered);


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
    

    FILE * hits_file = NULL, * hits_filtered = NULL, * hits_file_headers = NULL;
    FILE * hits_temp = NULL;

    init_args(argc, argv, &hits_file, &hits_file_headers, &hits_filtered);

    
    uint64_t i;

    // Determine number of frames
    uint64_t n_frames = 0;
    char c = 'z';
    while(!feof(hits_file_headers)){
        c = fgetc(hits_file_headers);
        if(c == '\n') n_frames++;
    }
    rewind(hits_file_headers);

    fprintf(stdout, "[INFO] Detected %"PRIu64" frames\n", n_frames);

    // Determine frame distance and load frames
    uint64_t * frame_positions = (uint64_t *) malloc(n_frames * sizeof(uint64_t));
    uint64_t max_distance = 0;
    if(1 != fscanf(hits_file_headers, "%"PRIu64, &frame_positions[0])){ fprintf(stderr, "Could not read header\n"); exit(-1); }
    for(i=1; i<n_frames; i++){
        if(1 != fscanf(hits_file_headers, "%"PRIu64, &frame_positions[i])){ fprintf(stderr, "Could not read header\n"); exit(-1); }
        fprintf(stdout, "Frames: %"PRIu64"\n", frame_positions[i]);
        if(frame_positions[i] - frame_positions[i-1] > max_distance) max_distance = frame_positions[i] - frame_positions[i-1];
    }

    // Allocate buffers
    
    
    hit * buffer_A = (hit *) malloc(max_distance * sizeof(hit));
    if(buffer_A == NULL) { fprintf(stderr, "Could not allocate second level loop of buffer frames on A\n"); exit(-1); }
    hit * buffer_B = (hit *) malloc(max_distance * sizeof(hit));
    if(buffer_B == NULL) { fprintf(stderr, "Could not allocate second level loop of buffer frames on B\n"); exit(-1); }

    /*
    fseek(hits_file, 0, SEEK_END);
    uint64_t t_hits = ftell(hits_file)/sizeof(hit);
    rewind(hits_file);
    */

    // Get first divisible by two
    uint64_t first_by_two = n_frames - (n_frames % 2);

    uint64_t i, t_to_write, write_A, write_B;
    char temp_name[2048];
    for(i=0; i<first_by_two; i++){

        // Read hits
        t_to_write = 0;
        write_A = fread(buffer_A, sizeof(hit), frame_positions[i+1]-frame_positions[i], hits_file);
        t_to_write += write_A;
        write_B = fread(buffer_B, sizeof(hit), frame_positions[i+2]-frame_positions[i+1], hits_file);
        t_to_write += write_B;

        // Open output temp file
        temp_name[0] = '\0';
        sprintf(temp_name, "hits-%"PRIu64, i);
        hits_temp = fopen(temp_name, "wb");
        if(hits_temp == NULL) { fprintf(stderr, "Could not open temp hits %s\n", temp_name); exit(-1); }

        // Merge
        uint64_t j1 = 0, j2 = 0;
        while((j1+j2) < t_to_write){

            if(j1 == write_A){
                fwrite(&buffer_A[j2], sizeof(hit), 1, hits_temp); ++j2;
            }else if(j2 == write_B){
                fwrite(&buffer_A[j1], sizeof(hit), 1, hits_temp); ++j1;
            }else{
                if(buffer_A[j1].diag < buffer_B[j2].diag){
                    fwrite(&buffer_A[j1], sizeof(hit), 1, hits_temp); ++j1;
                }else{
                    fwrite(&buffer_A[j2], sizeof(hit), 1, hits_temp); ++j2;
                }
            }
        }

    }



















    // Read items
    uint64_t read;
    for(i=0; i<n_frames-1; i++){
        buffers[i].curr_pos = frame_positions[i];
        buffers[i].i = 0;
        buffers[i].limit = frame_positions[i+1];
        fseek(hits_file, buffers[i].curr_pos * sizeof(hit), SEEK_SET);
        read = fread(buffers[i].h, sizeof(hit), buffer_size, hits_file);
        if(read == 0) fprintf(stdout, "Read zero items\n");
    }
    buffers[n_frames-1].curr_pos = frame_positions[n_frames-1];
    fseek(hits_file, buffers[i].curr_pos * sizeof(hit), SEEK_SET);
    read = fread(buffers[i].h, sizeof(hit), buffer_size, hits_file);
    buffers[n_frames-1].limit = t_hits;
    buffers[n_frames-1].i = 0;
    
    // Create queue
    fixed_priority_queue * fpq = new fixed_priority_queue(queue_size, hits_filtered);

    // Loop
    uint64_t buffers_finished = 0;
    while( buffers_finished < n_frames ){

        int64_t minimum_value = 9223372036854775807;
        int64_t remember_queue = -1;
        // Check if a buffer is done 
        uint64_t t;
        for(t=0; t<n_frames; t++){
            // If we reached end of buffer and there is still more to read
            if((buffers[t].i + buffers[t].curr_pos) == buffers[t].limit) buffers_finished++;
            if(buffers[t].i == buffer_size && (buffers[t].i + buffers[t].curr_pos) < buffers[t].limit){
                // Read another piece
                buffers[t].curr_pos += buffers[t].i;
                fseek(hits_file, buffers[t].curr_pos * sizeof(hit), SEEK_SET);
                read = fread(buffers[i].h, sizeof(hit), buffer_size, hits_file);
                buffers[t].i = 0;
            }

            

            if(minimum_value > buffers[t].h[buffers[t].i].diag){ minimum_value = buffers[t].h[buffers[t].i].diag; remember_queue = t; }
            fpq->insert(&buffers[t].h[buffers[t].i]);
            //printf("prev buf[t].i = %"PRIu64"\n", buffers[t].i);
            buffers[t].i++;
            //printf("buf[t].i = %"PRIu64"\n", buffers[t].i);
            
        }
        //printf("minimum found: \t%"PRId64"\n", minimum_value);
        // Re apply to queue with the minimum
        if(remember_queue != -1){
                fpq->insert(&buffers[remember_queue].h[buffers[remember_queue].i]); buffers[remember_queue].i++;
        }
        
        // Write to disk
        fpq->liberate_hits(minimum_value);

    }


    
    free(frame_positions);

    for(i=0; i<n_frames; i++){
        free(&buffers[i]);
    }
    free(buffers);
    fclose(hits_filtered);
    
    fprintf(stdout, "[INFO] Completed processing query splits\n");
    
    return 0;
}


void init_args(int argc, char ** av, FILE ** hits_file, FILE ** hits_file_headers, FILE ** hits_filtered){
    
    int pNum = 0;
    
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           mergehits -f [file] -o [file]\n");
            fprintf(stdout, "OPTIONAL:\n");
            
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }

        if(strcmp(av[pNum], "-f") == 0){
            *hits_file = fopen(av[pNum+1], "rb");
            if(*hits_file==NULL){ fprintf(stderr, "Could not open input hits file\n"); exit(-1); }
            char temp[2048]; temp[0] = '\0';
            strcpy(temp, av[pNum+1]);
            strcat(temp, ".headers");
            *hits_file_headers = fopen(temp, "rt");
            if(*hits_file_headers==NULL){ fprintf(stderr, "Could not open input hits headers file\n"); exit(-1); }
        }

        if(strcmp(av[pNum], "-o") == 0){
            *hits_filtered = fopen(av[pNum+1], "wb");
            if(*hits_filtered==NULL){ fprintf(stderr, "Could not open output filtered hits file\n"); exit(-1); }
        }

        pNum++;

    }   
    
    if(*hits_file==NULL || *hits_filtered==NULL){ fprintf(stderr, "You have to include an input and output hits file!\n"); exit(-1); }
}

char * get_dirname(char * path){
    int pos_last = 0, i = 0;
    while(path[i] != '\0'){
        if(path[i] == '/'){
            pos_last = i;
        }
        ++i;
    }
    char * dirname = (char *) malloc(BUFFER_SIZE * sizeof(char));
    if(dirname == NULL){ fprintf(stderr, "Could not allocate dirname char\n"); exit(-1); }

    memcpy(&dirname[0], &path[0], pos_last);
    dirname[pos_last] = '\0';

    return dirname;
}

char * get_basename(char * path){
    char * s = strrchr(path, '/');
    if (!s) return strdup(path); else return strdup(s + 1);
}