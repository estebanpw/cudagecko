
CC=gcc
CXX=g++
OPENCL=/usr/local/cuda-8.0/include/
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: gpu_words gpu_sort_words gpu_sort_hits gpu_frags gpu_print_words gpu_print_hits gpu_convert_frags_csv gpu_merge_hits gpu_sort_hits_plus_frags gpu_workflow

gpu_words: gpu_words.c
	$(CXX) $(CFLAGS) -I$(OPENCL) gpu_words.c -l OpenCL -o $(BIN)/gpu_words

gpu_sort_words: gpu_sort.c
	$(CXX) $(CFLAGS) -I$(OPENCL) gpu_sort.c -l OpenCL -o $(BIN)/gpu_sort_words

gpu_sort_hits: gpu_sort_hits.c
	$(CXX) $(CFLAGS) -I$(OPENCL) gpu_sort_hits.c -l OpenCL -o $(BIN)/gpu_sort_hits

gpu_sort_hits_plus_frags: gpu_sort_hits_plus_frags.c
	$(CXX) $(CFLAGS) -I$(OPENCL) gpu_sort_hits_plus_frags.c -l OpenCL commonFunctions.c -o $(BIN)/gpu_sort_hits_plus_frags

gpu_merge_hits: merge_hits.c
	$(CXX) $(CFLAGS) -I$(OPENCL) merge_hits.c -l OpenCL -o $(BIN)/merge_hits

gpu_frags: gpu_frags.c
	$(CXX) $(CFLAGS) -I$(OPENCL) gpu_frags.c -l OpenCL commonFunctions.c -o $(BIN)/gpu_frags

gpu_print_words: gpu_print_words_file.c
	$(CXX) $(CFLAGS) -I$(OPENCL) gpu_print_words_file.c -l OpenCL -o $(BIN)/gpu_print_words

gpu_print_hits: gpu_print_hits_file.c
	$(CXX) $(CFLAGS) -I$(OPENCL) gpu_print_hits_file.c -l OpenCL -o $(BIN)/gpu_print_hits

gpu_convert_frags_csv: filterFrags.c
	$(CXX) $(CFLAGS) -I$(OPENCL) filterFrags.c commonFunctions.c -l OpenCL -lm -o $(BIN)/gpu_convert_frags_csv

gpu_workflow: gpu_workflow.c
	$(CXX) $(CFLAGS) -I$(OPENCL) gpu_workflow.c -l OpenCL commonFunctions.c -o $(BIN)/gpu_workflow