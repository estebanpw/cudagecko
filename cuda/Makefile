

CC=gcc
CXX=g++ -mavx2 -fabi-version=0
NVIDIAC=nvcc
OPENCL=/usr/local/cuda-8.0/include/
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE #-DVERBOSE
BIN=../bin
NVIDIAFLAGS=-rdc=true -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -DCACHED -DAVX512CUSTOM #-DSHOWTIME #--ptxas-options=-v #-DVERBOSE

MODERNGPUFLAGS=-I ./moderngpu/src -std=c++11 --expt-extended-lambda -use_fast_math

all: gpu_cuda_workflow get_alignments

gpu_cuda_workflow: gpu_cuda_workflow.cu
	$(NVIDIAC) $(NVIDIAFLAGS) $(MODERNGPUFLAGS) gpu_cuda_workflow.cu kernels.cu -o $(BIN)/gpu_cuda_workflow

get_alignments: get_alignments.c
	$(CXX) $(CFLAGS) get_alignments.c common.c -o $(BIN)/gpu_get_alignments
