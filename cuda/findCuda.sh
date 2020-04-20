#!/bin/bash

# create a 'here document' that is code we compile and use to probe the card
cat << EOF > cudaComputeVersion.cu
#include <stdio.h>
int main()
{
cudaDeviceProp prop;
cudaGetDeviceProperties(&prop,0);
int v = prop.major * 10 + prop.minor;
printf("-gencode arch=compute_%d,code=sm_%d\n",v,v);
}
EOF

# probe the card and cleanup
nvcc cudaComputeVersion.cu -o cudaComputeVersion
./cudaComputeVersion
rm cudaComputeVersion.cu
rm cudaComputeVersion
