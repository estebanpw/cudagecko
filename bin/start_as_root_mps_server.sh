#!/bin/bash
# the following must be performed with root privilege
export CUDA_VISIBLE_DEVICES="0"
# WARNING: Change the "-i 0" value to the device id that you want to use
nvidia-smi -i 0 -c EXCLUSIVE_PROCESS
nvidia-cuda-mps-control -d
